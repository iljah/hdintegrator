#! /usr/bin/env python3
'''
High-dimensional mathematical function integrator.

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

import argparse
from shlex import split
from subprocess import Popen, PIPE


if __name__ == '__main__':

	parser = argparse.ArgumentParser(
		description = 'Integrate a mathematical function',
		formatter_class = argparse.ArgumentDefaultsHelpFormatter
	)

	parser.add_argument(
		'--integrator',
		required = True,
		help = 'Path to integrator program to use, relative to current working directory'
	)
	parser.add_argument(
		'--dimensions',
		required = True,
		type = int,
		help = 'Number of dimensions to use'
	)
	parser.add_argument(
		'--min-extent',
		type = float,
		default = 0,
		help = 'Minimum extent of integration in every dimension'
	)
	parser.add_argument(
		'--max-extent',
		type = float,
		default = 1,
		help = 'Maximum extent of integration in every dimension'
	)
	parser.add_argument(
		'--args',
		help = 'Arguments to pass to integrator program, given as one string (e.g. quoted) which are passed on to integrator after splitting with shlex.split'
	)

	args = parser.parse_args()

	if args.dimensions < 1:
		exit('Number of dimensions must be at least 1')

	extents = [(args.min_extent, args.max_extent) for i in range(args.dimensions)]

	arg_list = [args.integrator]
	if args.args != None:
		arg_list += split(args.args)
	integrator = Popen(arg_list, stdin = PIPE, stdout = PIPE, universal_newlines = True, bufsize = 1)

	for extent in extents:
		integrator.stdin.write('{:.15e} {:.15e} '.format(extent[0], extent[1]))
	integrator.stdin.write('\n')
	integrator.stdin.flush()

	print(integrator.stdout.readline(), end = '')

	integrator.kill()
	integrator.communicate()
