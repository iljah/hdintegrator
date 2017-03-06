/*
Integrator program for N-sphere using gsl plain mc to do the work.

Copyright 2017 Ilja Honkonen

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.
*/


#include "algorithm"
#include "cmath"
#include "cstdlib"
#include "iomanip"
#include "ios"
#include "iostream"
#include "sstream"
#include "string"
#include "vector"

#include "gsl/gsl_monte_plain.h"


/*
Same as integrand in example1.py
*/
double integrand(double* x, size_t dim, void* /*params*/) {
	double arg_to_sqrt = 1;
	for (size_t i = 0; i < dim; i++) {
		arg_to_sqrt -= (*(x + i)) * (*(x + i));
	}
	return std::sqrt(std::max(0.0, arg_to_sqrt));
}


/*
Reads integration volume from stdin and prints the result to stdout.
*/
int main(int argc, char* argv[])
{
	if (argc != 2) {
		std::cerr << "Invalid number of arguments: " << argc - 1 << " should be 1 (number of calls)" << std::endl;
		return EXIT_FAILURE;
	}

	const double calls = std::atof(argv[1]);

	std::cout << std::setprecision(15) << std::scientific;

	gsl_rng_env_setup();
	const gsl_rng_type* const rng_t = gsl_rng_default;
	gsl_rng* rng = gsl_rng_alloc(rng_t);

	gsl_monte_function function;
	function.f = &integrand;
	function.params = nullptr;

	size_t dimensions = 0;
	bool first_integration = true;

	decltype(gsl_monte_plain_alloc(0)) state{};

	std::string line;
	while (std::getline(std::cin, line)) {

		std::istringstream iss(line);
		std::vector<double> mins, maxs;
		double item;
		while (iss >> item) {
			mins.push_back(item);
			if (iss >> item) {
				maxs.push_back(item);
			} else {
				break;
			}
		}
		if (mins.size() == 0) {
			break;
		}
		if (mins.size() != maxs.size()) {
			std::cerr << "Number of minimum and maximum extents differs" << std::endl;
			return EXIT_FAILURE;
		}

		if (first_integration) {
			first_integration = false;
			state = gsl_monte_plain_alloc(mins.size());
		} else if (dimensions != mins.size()) {
			gsl_monte_plain_free(state);
			state = gsl_monte_plain_alloc(mins.size());
		}
		dimensions = mins.size();

		function.dim = dimensions;
		double result = 0, abserr = 0;
		const auto ret_val = gsl_monte_plain_integrate(
			&function,
			mins.data(),
			maxs.data(),
			dimensions,
			size_t(std::round(calls)),
			rng,
			state,
			&result,
			&abserr
		);
		if (ret_val != 0) {
			std::cerr << "Integration failed." << std::endl;
			return EXIT_FAILURE;
		}
		std::cout << result << " " << abserr << std::endl;
	}

	if (not first_integration) {
		gsl_monte_plain_free(state);
	}
}
