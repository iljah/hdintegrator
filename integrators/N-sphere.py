#! /usr/bin/env python3
'''
Integrator program for N-sphere using scipy nquad to do the work.

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


from random import randint
from sys import stdin, stdout

from scipy.integrate import nquad


'''
N-sphere.

Returns value of scalar function f at position r = (r1, r2, ..., rN)
representing an N-sphere of radius 1 centered at origin:
f(r1, r2, ..., rN) = (1 - r1**2 - r2**2 - ... - rN**2)**0.5.
For example in two dimensions (1 input + 1 output), in other words
for a circle: f(r1) = (1 - r1**2)**0.5.

Returns 0 if r1**2 + ... + rN**2 > 1.

r1 = args[0], r2 = args[1], etc.

If used for integrating ([-1,1]) result will be half area of unit circle,
if used for ([-1,1],[-1,1]) result will be half volume of unit sphere, etc...
'''
def integrand(*args):
	arg_for_sqrt = 1.0
	for r in args:
		arg_for_sqrt -= r**2
	return max(0, arg_for_sqrt)**0.5


'''
Reads integration volume from stdin and prints the result to stdout.

Input format, line by line:
nr_calls v0min v0max v1min v1max ...

Ignores nr_calls.
'''
if __name__ == '__main__':

	while True:
		instr = stdin.readline()
		if instr == '':
			break
		instr = instr.strip().split()
		extents = []
		for i in range(2, len(instr), 2):
			extents.append((float(instr[i - 1]), float(instr[i])))
		result, error = nquad(integrand, extents)
		split_dim = 0
		if len(extents) > 2:
			split_dim = randint(0, int(len(extents)/2) - 1)
		stdout.write('{:.15e} {:.15e} {:d}\n'.format(result, error, split_dim))
		stdout.flush()
