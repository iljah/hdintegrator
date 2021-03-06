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


size_t index(const size_t x_i, const size_t t_i, const size_t nx, const size_t nt)
{
	return x_i % nx + (t_i % nt) * nx;
}


/*
Integrand for single-time correlation function of burgers equation.
*/
double integrand(double* x, size_t dimensions, void* integrand_params) {
	const auto& params = *static_cast<Integrand_Params*>(integrand_params);
	const auto
		corr1 = params.corr1,
		corr2 = params.corr2;
	const auto
		nx = params.nx,
		nt = params.nt;

	if (corr1 >= int(dimensions) or corr2 >= int(dimensions)) {
		std::cerr << __FILE__ "(" << __LINE__ << ")" << std::endl;
		abort();
	}

	if (nx * nt != dimensions) {
		std::cerr << __FILE__ "(" << __LINE__ << ")" << std::endl;
		abort();
	}

	const bool correlate = [&](){
		if (corr1 < 0 or corr2 < 0) {
			return false;
		} else {
			return true;
		}
	}();

	// transformed version, integration range -1..1 instead of -inf..inf
	const std::vector<double> t = [x, dimensions](){
		std::vector<double> ret_val(x, x + dimensions);
		for (auto& i: ret_val) {
			i = i / (1 - i*i);
		}
		return ret_val;
	}();

	const double
		transform_factor = [&](){
			double ret_val = 1;
			for (size_t t_i = 0; t_i < nt; t_i++) {
			for (size_t x_i = 0; x_i < nx; x_i++) {
				const double x2 = SQR(x[index(x_i, t_i, nx, nt)]);
				ret_val *= (1 + x2) / SQR(1 - x2);
			}}
			return ret_val;
		}(),
		vel1 = [&]{
			if (correlate) {
				return t[index(corr1, 0, nx, nt)];
			} else {
				return 1.0;
			}
		}(),
		vel2 = [&]{
			if (correlate) {
				return t[index(corr2, 0, nx, nt)];
			} else {
				return 1.0;
			}
		}(),
		arg4exp = [&](){
			double ret_val = 0;
			for (size_t t_i = 0; t_i < nt; t_i++) {
			for (size_t x_i = 0; x_i < nx; x_i++) {
				ret_val += SQR(
					+ t[index(x_i         , t_i + 1, nx, nt)]
					+ t[index(x_i         , t_i    , nx, nt)]
					- t[index(x_i + 1     , t_i    , nx, nt)]
					- t[index(x_i + nx - 1, t_i    , nx, nt)]
					+ 0.5 * t[index(x_i, t_i, nx, nt)]
						* (
							+ t[index(x_i + 1     , t_i, nx, nt)]
							- t[index(x_i + nx - 1, t_i, nx, nt)]
						)
				);
			}}
			return ret_val;
		}();

	return transform_factor * vel1 * vel2 * exp(-0.5 * arg4exp);
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

		for (size_t i = 0; i < mins.size(); i++) {
			if (mins[i] >= maxs[i]) {
				std::cerr << "Starting coordinate of " << i+1
					<< "th dimension is not smaller than ending coordinate: "
					<< mins[i] << " >= " << maxs[i]
					<< std::endl;
				return EXIT_FAILURE;
			}
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

		Integrand_Params params{corr1, corr2, nx, nt};
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
