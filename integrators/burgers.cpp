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
#include "array"
#include "cmath"
#include "cstdlib"
#include "iomanip"
#include "ios"
#include "iostream"
#include "iterator"
#include "sstream"
#include "stdexcept"
#include "string"
#include "vector"

#include "boost/program_options.hpp"

#if METHOD == 1
#include "gsl_monte_plain2.h"
#include "gsl/gsl_monte_plain.h"
#elif METHOD == 2
#include "gsl_monte_miser2.h"
#include "gsl/gsl_monte_miser.h"
#elif METHOD == 3
#include "gsl_monte_vegas2.h"
#include "gsl/gsl_monte_vegas.h"
#else
#error You must choose a method when compiling (e.g. -DMETHOD=1)
#endif


template<class T> constexpr T SQR(const T& t)
{
	return t*t;
}


struct Integrand_Params
{
	// correlate in these dimensions, (nx-1)*nt - 1...nx*nt-1
	int corr1 = -1, corr2 = -1;
	size_t nx = 0, nt = 0;
};


size_t index(const size_t x_i, const size_t t_i, const size_t nx)
{
	return x_i % nx + t_i * nx;
}


/*
Integrand for single-time correlation function of burgers equation.
*/
double integrand(double* x, size_t dimensions, void* integrand_params) {
	const auto& params = *static_cast<Integrand_Params*>(integrand_params);

	if (params.corr1 >= int(dimensions) or params.corr2 >= int(dimensions)) {
		std::cerr << __FILE__ "(" << __LINE__ << ")" << std::endl;
		abort();
	}

	const bool correlate = [&](){
		if (params.corr1 < 0 or params.corr2 < 0) {
			return false;
		} else {
			return true;
		}
	}();

	const double
		vel1 = [&]{
			if (correlate) {
				return x[index(params.corr1, params.nt - 1, params.nx)];
			} else {
				return 1.0;
			}
		}(),
		vel2 = [&]{
			if (correlate) {
				return x[index(params.corr2, params.nt - 1, params.nx)];
			} else {
				return 1.0;
			}
		}(),
		arg4exp = [&](){
			double ret_val = 0;
			for (size_t x_i = 0; x_i < params.nx; x_i++) {
			for (size_t t_i = 0; t_i < params.nt - 1; t_i++) {
				ret_val += SQR(
					+ x[index(x_i, t_i + 1, params.nx)]
					- x[index(x_i, t_i    , params.nx)]
					+ 0.5 * x[index(x_i, t_i, params.nx)]
						* (
							+ x[index(x_i + 1, t_i, params.nx)]
							- x[index(x_i - 1, t_i, params.nx)]
						)
					- (
						+     x[index(x_i + 1, t_i, params.nx)]
						- 2 * x[index(x_i    , t_i, params.nx)]
						+     x[index(x_i - 1, t_i, params.nx)]
					)
				);
			}}
			return ret_val;
		}();

	return vel1 * vel2 * exp(-0.5 * arg4exp);
}


/*
Reads integration volume from stdin and prints the result to stdout.

Input format, line by line:
nr_calls v0min v0max v1min v1max ...
*/
int main(int argc, char* argv[])
{
	int corr1 = 0, corr2 = 0;
	size_t nx = 0, nt = 0;

	boost::program_options::options_description
		options("Usage: program_name [options], where options are");
	options.add_options()
		("help", "Print help")
		("corr1",
			boost::program_options::value<int>(&corr1)->required(),
			"Number of first correlation dimension starting from 0, not calculated if < 0")
		("corr2",
			boost::program_options::value<int>(&corr2)->required(),
			"Number of second correlation dimension starting from 0, not calculated if < 0")
		("nx",
			boost::program_options::value<size_t>(&nx)->required(),
			"Number of grid points in x direction, nx*nt must equal number of dimension given on stdin")
		("nt",
			boost::program_options::value<size_t>(&nt)->required(),
			"Number of grid points in t direction, nx*nt must equal number of dimension given on stdin");

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

	if (nx == 0) {
		std::cerr <<  __FILE__ << "(" << __LINE__ << "): "
			<< "Number of grid points in x direction must be > 0"
			<< std::endl;
		return EXIT_FAILURE;
	}
	if (nt == 0) {
		std::cerr <<  __FILE__ << "(" << __LINE__ << "): "
			<< "Number of grid points in t direction must be > 0"
			<< std::endl;
		return EXIT_FAILURE;
	}

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

	std::vector<int> split_dims;

	std::string line;
	while (std::getline(std::cin, line)) {

		double calls;
		std::istringstream iss(line);
		std::vector<double> mins, maxs;
		double item;
		iss >> calls;
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
		if (dimensions != nx * nt) {
			std::cerr << "Number of dimensions not equal to nx*nt" << std::endl;
			return EXIT_FAILURE;
		}

		function.dim = dimensions;

		Integrand_Params params;
		params.corr1 = corr1;
		params.corr2 = corr2;
		params.nx = nx;
		params.nt = nt;
		function.params = &params;

		std::vector<int> split_dims(dimensions);
		double result = 0, error = 0;
		#if METHOD == 1
		auto ret_val = gsl_monte_plain_integrate2(
		#elif METHOD == 2
		auto ret_val = gsl_monte_miser_integrate2(
		#elif METHOD == 3
		auto ret_val = gsl_monte_vegas_integrate2(
		#endif
			&function,
			mins.data(),
			maxs.data(),
			dimensions,
			size_t(std::round(calls)),
			rng,
			state,
			&result,
			&error,
			split_dims.data()
		);
		if (ret_val != 0) {
			std::cerr << "Integration failed." << std::endl;
			return EXIT_FAILURE;
		}

		const auto max_elem = std::max_element(split_dims.cbegin(), split_dims.cend());
		std::cout << result << " " << error << " " << std::distance(split_dims.cbegin(), max_elem) << std::endl;
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
