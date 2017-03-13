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
from sys import stdout

from cell import cell
from ndgrid import ndgrid


if __name__ == '__main__':

	parser = argparse.ArgumentParser(
		description = 'Integrate a mathematical function',
		formatter_class = argparse.ArgumentDefaultsHelpFormatter
	)

	parser.add_argument(
		'--verbose',
		action = 'store_true',
		help = 'Print diagnostic information during integration'
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
	parser.add_argument(
		'--calls',
		type = float,
		default = 1e6,
		help = 'Request this number of calls to integrand from integrator'
	)
	parser.add_argument(
		'--calls-factor',
		metavar = 'A',
		type = float,
		default = 2,
		help = 'Increase number of calls by factor F when checking for convergence'
	)
	parser.add_argument(
		'--convergence-factor',
		metavar = 'O',
		type = float,
		default = 1.05,
		help = 'Consider result converged when using A times more calls gives a result within factor O.'
	)
	parser.add_argument(
		'--convergence-diff',
		metavar = 'D',
		type = float,
		default = 1e-3,
		help = 'Consider result converged when using A times more calls gives a result within difference D.'
	)

	args = parser.parse_args()

	if args.dimensions < 1:
		exit('Number of dimensions must be at least 1')

	# prepare integrator program
	arg_list = [args.integrator]
	if args.args != None:
		arg_list += split(args.args)
	integrator = Popen(arg_list, stdin = PIPE, stdout = PIPE, universal_newlines = True, bufsize = 1)

	# initialize grid for integration
	c = cell()
	c.data['id'] = 1
	c.data['converged'] = False
	c.data['value'] = None
	c.data['error'] = None
	for i in range(args.dimensions):
		c.set_extent(i, args.min_extent, args.max_extent)
	grid = ndgrid(c)

	done = False
	while not done:

		done = True
		for c in grid.get_cells():

			if c.data['converged']:
				continue
			else:
				done = False

			dimensions = sorted(c.get_dimensions())

			if c.data['value'] == None:

				extents = [c.get_extent(dim) for dim in dimensions]
				if args.verbose:
					print('Initial processing of cell', c.data['id'], 'with', len(extents), 'dimensions:', end = ' ')
				integrator.stdin.write('{:.15e} '.format(args.calls))
				for extent in extents:
					if args.verbose:
						print(extent[0], extent[1], end = ' ')
					integrator.stdin.write('{:.15e} {:.15e} '.format(extent[0], extent[1]))
				integrator.stdin.write('\n')
				integrator.stdin.flush()
				stdout.flush()

				value, error = integrator.stdout.readline().strip().split()
				c.data['value'], c.data['error'] = float(value), float(error)

			# check convergence
			integrator.stdin.write('{:.15e} '.format(args.calls * args.convergence_factor))
			for extent in extents:
				integrator.stdin.write('{:.15e} {:.15e} '.format(extent[0], extent[1]))
			integrator.stdin.write('\n')
			integrator.stdin.flush()

			new_value, new_error = integrator.stdout.readline().strip().split()
			new_value, new_error = float(new_value), float(new_error)

			try:
				convg_fact = max(abs(c.data['value']), abs(new_value)) / min(abs(c.data['value']), abs(new_value))
			except:
				convg_fact = 0.0
			convg_diff = abs(c.data['value'] - new_value)
			c.data['value'] = new_value
			c.data['error'] = new_error

			if convg_fact < args.convergence_factor or convg_diff < args.convergence_diff:
				if args.verbose:
					print(', converged')
				c.data['converged'] = True
			else:
				if args.verbose:
					print(", didn't converge, splitting into cells ", end = '')
				c.data['value'] = c.data['error'] = None
				cells_to_split = [c]
				new_cs_to_split = []
				for dim in dimensions:
					if args.verbose:
						print('  (' + str(dim) + ')', end = ' ')
					for c_to_split in cells_to_split:
						old_id = c_to_split.data['id']
						for new_cell in grid.split(c_to_split, dim):
							new_cs_to_split.append(new_cell)
						new_cs_to_split[-2].data['id'] = old_id * 2
						new_cs_to_split[-1].data['id'] = old_id * 2 + 1
						if args.verbose:
							print(old_id * 2, old_id * 2 + 1, end = ' ')
					cells_to_split = new_cs_to_split
					new_cs_to_split = []
				if args.verbose:
					print()


value, error = 0.0, 0.0
for c in grid.get_cells():
	value += c.data['value']
	error += c.data['error']
print(value, error)
