#! /usr/bin/env python3
'''
Integrator program for turbulence using scipy nquad to do the work.

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
'''


from argparse import ArgumentParser
from math import exp
from sys import stdin, stdout

from scipy.integrate import nquad


'''
Integrand for scipy nquad function.

Calculates:
velN * velM * exp(sum)
where sum is:
v[0]**2*(v[-1]-v[1])**2 + v[1]**2*(v[0]-v[2])**2 + ... + v[-1]**2*(v[0]-v[-2])**2
with constraint:
v[0]+v[1]+v[2]+...+v[-1] = 0
and v[0] is replaced by other velocities so sum becomes:
(v[1]+...+v[-1])**2*(v[-1]-v[1])**2 + v[1]**2...
'''
def integrand(*args):
	if args[-1] < 0 or args[-2] < 0:
		vel1, vel2 = 1, 1
	else:
		vel1, vel2 = args[args[-2]], args[args[-1]]
	dimensions = len(args) - 2
	if dimensions <= args[-2] or dimensions <= args[-1]:
		exit('Too few dimension for calculating correlation between dimensions ' + str(args[-2]) + ' and ' + str(args[-1]))

	if dimensions == 1:
		exit('At least 2 dimensions required')

	all_sum = 0.0
	for i in range(dimensions):
		all_sum += args[i]

	non_boundary_term = 0
	for i in range(1, dimensions - 1):
		non_boundary_term += (args[i] * (args[i + 1] - args[i - 1]))**2

	return \
		vel1 \
		* vel2 \
		* exp(
			-(
				(all_sum * (args[0] - args[dimensions - 1]))**2
				+ (args[       0      ] * (all_sum + args[       1      ]))**2
				+ (args[dimensions - 1] * (all_sum + args[dimensions - 2]))**2
				+ non_boundary_term
			)
		)


'''
Reads integration volume from stdin and prints the result to stdout.
'''
if __name__ == '__main__':

	parser = ArgumentParser()
	parser.add_argument(
		'--corr1',
		required = True,
		metavar = 'C1',
		type = int,
		help = 'Calculate correlation between velocities at (zero-based) indices C1 and C2, not calculated if < 0'
	)
	parser.add_argument('--corr2', required = True, metavar = 'C2', type = int, help = 'See above')
	args = parser.parse_args()

	while True:
		instr = stdin.readline()
		if instr == '':
			break
		instr = instr.strip().split()
		extents = []
		for i in range(1, len(instr), 2):
			extents.append((float(instr[i - 1]), float(instr[i])))
		result, error1 = nquad(integrand, extents, args = [args.corr1, args.corr2])
		normalization, error2 = nquad(integrand, extents, args = [-1, -1])
		stdout.write(
			'{:.15e} {:.15e}\n'.format(
				result / normalization,
				((error1 / result)**2 + (error2 / normalization)**2)**0.5
			)
		)
		stdout.flush()
