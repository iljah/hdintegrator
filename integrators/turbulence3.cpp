/*
Integrator program for turbulence using cubature to do the work.

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
#include "array"
#include "cmath"
#include "cstdlib"
#include "iomanip"
#include "ios"
#include "iostream"
#include "sstream"
#include "string"
#include "vector"

#include "boost/program_options.hpp"
#include "cubature.h"


template<class T> constexpr T SQR(const T& t)
{
	return t*t;
}


int func(
	unsigned ndim,
	#ifdef VECTORIZE
	unsigned npts,
	#endif
	const double *x,
	void *fdata,
	unsigned /*dim*/,
	double *fval
) {
	#ifndef VECTORIZE
	constexpr unsigned npts = 1;
	#endif

	const auto
		vel1_i = *(static_cast<int*>(fdata) + 0),
		vel2_i = *(static_cast<int*>(fdata) + 1);

	if (vel1_i < 0 or vel2_i < 0) {
		for (size_t i = 0; i < npts; i++) {
			fval[i] = 1;
		}
	} else {
		for (size_t i = 0; i < npts; i++) {
			fval[i] = x[i*ndim + vel1_i] * x[i*ndim + vel2_i];
		}
	}

	std::vector<double> all_sum2(npts, 0), all_dx(npts, 0);
	for (size_t i = 0; i < npts; i++) {
		all_sum2[i] += SQR(x[i*ndim]);
		all_sum2[i] += SQR(x[i*ndim + ndim - 1]);
		all_dx[i]
			= SQR(-0.5 * x[i*ndim] * (x[i*ndim + 1] - x[i*ndim + ndim - 1]) + x[i*ndim + 1] - 2*x[i*ndim] + x[i*ndim + ndim - 1])
			+ SQR(-0.5 * x[i*ndim + ndim - 1] * (x[i*ndim] - x[i*ndim + ndim - 2]) + x[i*ndim] - 2*x[i*ndim + ndim - 1] + x[i*ndim + ndim - 2]);

		for (size_t d = 1; d < ndim - 1; d++) {
			all_sum2[i] += SQR(x[i*ndim + d]);
			all_dx[i] += SQR(
				-0.5 * x[i*ndim + d] * (x[i*ndim + d + 1] - x[i*ndim + d - 1])
				+ x[i*ndim + d + 1]
				+ x[i*ndim + d - 1]
				- 2*x[i*ndim + d]
			);
		}
	}

	for (size_t i = 0; i < npts; i++) {
		fval[i] *= exp(-0.5 * (all_sum2[i] + all_dx[i]));
	}

	return 0;
}


int main(int argc, char* argv[])
{
	int corr1 = 0, corr2 = 0;

	boost::program_options::options_description
		options("Usage: program_name [options], where options are");
	options.add_options()
		("help", "Print help")
		("corr1",
			boost::program_options::value<int>(&corr1)->required(),
			"Number of first correlation dimension starting from 0, not calculated if < 0")
		("corr2",
			boost::program_options::value<int>(&corr2)->required(),
			"Number of second correlation dimension starting from 0, not calculated if < 0");

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

	std::array<int, 2> corrs{corr1, corr2};

	std::cout << std::setprecision(15) << std::scientific;

	size_t dimensions = 0;

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

		dimensions = mins.size();

		double result = 0, abserr = 0;

		#ifdef PCUBATURE
		#ifdef VECTORIZE
		const auto ret_val = pcubature_v(
		#else
		const auto ret_val = pcubature(
		#endif
		#else
		#ifdef VECTORIZE
		const auto ret_val = hcubature_v(
		#else
		const auto ret_val = hcubature(
		#endif
		#endif
			1,
			func,
			corrs.data(),
			dimensions,
			mins.data(),
			maxs.data(),
			size_t(std::round(calls)),
			0,
			0,
			ERROR_L1,
			&result,
			&abserr
		);

		if (ret_val != 0) {
			std::cerr << "Integration failed." << std::endl;
			return EXIT_FAILURE;
		}

		std::cout << result << " " << abserr << std::endl;
	}

	return EXIT_SUCCESS;
}
