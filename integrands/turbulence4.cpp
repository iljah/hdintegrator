/*
Integrator program for turbulence using mt19937.

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
#include "random"
#include "sstream"
#include "string"
#include "vector"

#include "boost/program_options.hpp"
#include "prettyprint.hpp"


template<class T> constexpr T SQR(const T& t)
{
	return t*t;
}

/*
Returns value of integral and dimension with largest difference in average value on both sides of middle.

Calculates integral within given volume, if corr1 or corr2 > 0
calculates correlation integral of given dimensions.

Supported Scalars are at least float, double and vector versions
of them from Agner's vector class.
*/
template<
	class Scalar
> std::pair<
	Scalar,
	size_t
> integrate(
	const std::vector<Scalar>& r_min,
	const std::vector<Scalar>& r_max,
	const int corr1,
	const int corr2,
	const size_t samples
) {
	using std::abs;
	using std::exp;

	if (r_min.size() != r_max.size()) {
		std::cerr << __FILE__ "(" << __LINE__ << ")" << std::endl;
		abort();
	}

	const auto ndims = r_min.size();

	std::vector<std::uniform_real_distribution<>> dists;
	dists.reserve(ndims);
	for (size_t i = 0; i < ndims; i++) {
		dists.emplace_back(r_min[i], r_max[i]);
	}

	Scalar average{0};

	std::vector<std::array<Scalar, 2>> quad_avgs(ndims, {0, 0});
	std::vector<std::array<size_t, 2>> quad_nr(ndims, {0, 0});

	std::mt19937 rnd;
	std::vector<Scalar> r(ndims);
	for (size_t sample = 0; sample < samples; sample++) {

		size_t i = 0;
		std::generate(
			r.begin(),
			r.end(),
			[&](){ return dists[i++](rnd); }
		);

		auto result
			= [&]()->Scalar {
				if (corr1 < 0 or corr2 < 0) {
					return 1;
				} else if (size_t(corr1) >= ndims or size_t(corr2) >= ndims) {
					std::cout << __FILE__ "(" << __LINE__ << ")" << std::endl;
					abort();
				} else {
					return r[size_t(corr1)] * r[size_t(corr2)];
				}
			}();

		Scalar
			all_sum2 = SQR(r[0]) + SQR(r[ndims - 1]),
			all_dx
				= SQR(-0.5 * r[0] * (r[1] - r[ndims - 1]) + r[1] - 2*r[0] + r[ndims - 1])
				+ SQR(-0.5 * r[ndims - 1] * (r[0] - r[ndims - 2]) + r[0] - 2*r[ndims - 1] + r[ndims - 2]);

		for (size_t d = 1; d < ndims - 1; d++) {
			all_sum2 += SQR(r[d]);
			all_dx += SQR(
				-0.5 * r[d] * (r[d + 1] - r[d - 1])
				+ r[d + 1]
				+ r[d - 1]
				- 2*r[d]
			);
		}

		result *= exp(-0.5 * (all_sum2 + all_dx));
		average += result;

		for (size_t d = 0; d < ndims; d++) {
			if (r[d] - r_min[d] < r_max[d] - r[d]) {
				quad_avgs[d][0] += result;
				quad_nr[d][0]++;
			} else {
				quad_avgs[d][1] += result;
				quad_nr[d][1]++;
			}
		}
	}

	average /= samples;

	Scalar max_diff = -1;
	size_t max_diff_d = 0;
	for (size_t d = 0; d < ndims; d++) {
		quad_avgs[d][0] /= quad_nr[d][0];
		quad_avgs[d][1] /= quad_nr[d][1];
		const Scalar diff = abs(quad_avgs[d][0] - quad_avgs[d][1]);
		if (max_diff < diff) {
			max_diff = diff;
			max_diff_d = d;
		}
	}

	return {average, max_diff_d};
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

	std::cout << std::setprecision(15) << std::scientific;

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

		const auto [result, max_diff_dim] = integrate(mins, maxs, corr1, corr2, calls);

		std::cout << result << " -1 " << max_diff_dim << std::endl;
	}

	return EXIT_SUCCESS;
}
