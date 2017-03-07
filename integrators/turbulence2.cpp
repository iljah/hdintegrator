/*
Integrator program for turbulence using gsl mc to do the work.

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


/*
Selects which method to use.

1 == plain
2 == miser
3 == vegas
*/
#ifndef METHOD
#define METHOD 1
#endif


#include "algorithm"
#include "cmath"
#include "cstdlib"
#include "iomanip"
#include "ios"
#include "iostream"
#include "sstream"
#include "string"
#include "vector"

#include "boost/program_options.hpp"

#if METHOD == 1
#include "gsl/gsl_monte_plain.h"
#elif METHOD == 2
#include "gsl/gsl_monte_miser.h"
#elif METHOD == 3
#include "gsl/gsl_monte_vegas.h"
#endif


template<class T> T SQR(const T& t)
{
	return t*t;
}


/*
Same as integrand in example1.py
*/
double integrand(double* x, size_t dimensions, void* params) {
	const double
		vel1 = [x, params]{
			const auto vel1_i = *(static_cast<int*>(params) + 0);
			if (vel1_i < 0) {
				return 1.0;
			} else {
				return *(x + vel1_i);
			}
		}(),
		vel2 = [x, params]{
			const auto vel2_i = *(static_cast<int*>(params) + 1);
			if (vel2_i < 0) {
				return 1.0;
			} else {
				return *(x + vel2_i);
			}
		}(),
		all_sum2 = [x, dimensions](){
			double sum2 = 1;
			for (size_t i = 0; i < dimensions; i++) {
				sum2 *= std::exp(-0.5 * SQR(x[i]));
			}
			return sum2;
		}(),
		all_dx = [x, dimensions](){
			double dx
				= std::exp(-0.5 * SQR(x[0] * (x[1] - x[dimensions - 1]) + x[1] - 2*x[0] + x[dimensions - 1]))
				* std::exp(-0.5 * SQR(x[dimensions - 1] * (x[0] - x[dimensions - 2]) + x[0] - 2*x[dimensions - 1] + x[dimensions - 2]));
			for (size_t i = 1; i < dimensions - 1; i++) {
				dx *= std::exp(-0.5 * SQR(x[i] * (x[i+1] - x[i-1]) + x[i+1] - 2*x[i] + x[i-1]));
			}
			return dx;
		}();

	return vel1 * vel2 * all_sum2 * all_dx;
}


/*
Reads integration volume from stdin and prints the result to stdout.
*/
int main(int argc, char* argv[])
{
	int corr1 = 0, corr2 = 0;
	double calls = 0;

	boost::program_options::options_description
		options("Usage: program_name [options], where options are");
	options.add_options()
		("help", "Print help")
		("calls",
			boost::program_options::value<double>(&calls)->required(),
			"Number of calls to make when integrating")
		("corr1",
			boost::program_options::value<int>(&corr1),
			"Number of first correlation dimension starting from 0")
		("corr2",
			boost::program_options::value<int>(&corr2),
			"Number of second correlation dimension starting from 0");

	boost::program_options::variables_map var_map;
	try {
		boost::program_options::store(
			boost::program_options::parse_command_line(
				argc,
				argv,
				options,
				boost::program_options::command_line_style::unix_style ^ boost::program_options::command_line_style::allow_short
			),
			var_map
		);
	} catch (std::exception& e) {
		std::cerr <<  __FILE__ << "(" << __LINE__ << "): "
			<< "Couldn't parse command line options: " << e.what()
			<< std::endl;
		return EXIT_FAILURE;
	}
	boost::program_options::notify(var_map);

	if (var_map.count("help") > 0) {
		std::cout << options << std::endl;
		return EXIT_SUCCESS;
	}


	gsl_rng_env_setup();
	const gsl_rng_type* rng_t = gsl_rng_default;
	gsl_rng* rng = gsl_rng_alloc(rng_t);

	std::cout << std::setprecision(15) << std::scientific;

	gsl_monte_function function;
	function.f = &integrand;

	size_t dimensions = 0;
	bool first_integration = true;

	#if METHOD == 1
	decltype(gsl_monte_plain_alloc(0)) state{};
	#elif METHOD == 2
	decltype(gsl_monte_miser_alloc(0)) state{};
	#elif METHOD == 3
	decltype(gsl_monte_vegas_alloc(0)) state{};
	#endif

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
			#if METHOD == 1
			state = gsl_monte_plain_alloc(mins.size());
			#elif METHOD == 2
			state = gsl_monte_miser_alloc(mins.size());
			#elif METHOD == 3
			state = gsl_monte_vegas_alloc(mins.size());
			#endif
		} else if (dimensions != mins.size()) {
			#if METHOD == 1
			gsl_monte_plain_free(state);
			state = gsl_monte_plain_alloc(mins.size());
			#elif METHOD == 2
			gsl_monte_miser_free(state);
			state = gsl_monte_miser_alloc(mins.size());
			#elif METHOD == 3
			gsl_monte_vegas_free(state);
			state = gsl_monte_vegas_alloc(mins.size());
			#endif
		}
		dimensions = mins.size();

		function.dim = dimensions;

		std::array<int, 2> corrs{corr1, corr2};
		function.params = corrs.data();

		double result = 0, error = 0;
		#if METHOD == 1
		auto ret_val = gsl_monte_plain_integrate(
		#elif METHOD == 2
		auto ret_val = gsl_monte_miser_integrate(
		#elif METHOD == 3
		auto ret_val = gsl_monte_vegas_integrate(
		#endif
			&function,
			mins.data(),
			maxs.data(),
			dimensions,
			size_t(std::round(calls)),
			rng,
			state,
			&result,
			&error
		);
		if (ret_val != 0) {
			std::cerr << "Integration failed." << std::endl;
			return EXIT_FAILURE;
		}

		std::cout << result << " " << error << std::endl;
	}

	if (not first_integration) {
		#if METHOD == 1
		gsl_monte_plain_free(state);
		#elif METHOD == 2
		gsl_monte_miser_free(state);
		#elif METHOD == 3
		gsl_monte_vegas_free(state);
		#endif
	}

	return EXIT_SUCCESS;
}
