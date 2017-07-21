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
from os.path import exists
from random import choice, randint
import shlex
from subprocess import Popen, PIPE
from sys import stdout
from time import sleep

from cell import cell
from mpi4py import MPI
from ndgrid import ndgrid


'''
Splits given cell splits times in given dimensions.

ID of child cell is parent id * 2 + (0 or 1).
'''
def split(cell, splits, dimensions, grid):
	# TODO logging
	cells_to_split = [cell]
	new_cells_to_split = []
	for dim in dimensions:
		for i in range(splits):
			for c_to_split in cells_to_split:
				old_id = c_to_split.data['id']
				for new_cell in grid.split(c_to_split, dim):
					new_cells_to_split.append(new_cell)
				new_cells_to_split[-2].data['id'] = old_id * 2
				new_cells_to_split[-1].data['id'] = old_id * 2 + 1
			cells_to_split = new_cells_to_split
			new_cells_to_split = []


class Work_Item:
	def __init__(self):
		self.volume = None
		self.cell_id = None
		self.converged = None
		self.value = None
		self.error = None
		self.split_dim = None


class Work_Tracker:
	def __init__(self):
		self.item = None
		self.processing = None


if __name__ == '__main__':

	comm = MPI.COMM_WORLD
	rank = comm.Get_rank()
	if comm.size < 2:
		if rank == 0:
			print('At least 2 processes required')
		exit(1)

	if rank == 0:
		print('Starting with', comm.size, 'processes')
		stdout.flush()

	parser = argparse.ArgumentParser(
		description = 'Integrate a mathematical function',
		formatter_class = argparse.ArgumentDefaultsHelpFormatter,
		add_help = False
	)

	parser.add_argument(
		'-h', '--help',
		action = 'store_true',
		help = 'Show this help message and exit'
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
		help = 'Arguments to pass to integrator program, given as one string (e.g. --args "-a b -c d") which are passed on to integrator after splitting with shlex.split'
	)
	parser.add_argument(
		'--prerefine',
		type = int,
		default = 0,
		metavar = 'S',
		help = 'Split S times a random grid cell in random dimension before integrating'
	)
	parser.add_argument(
		'--calls',
		type = float,
		default = 1e6,
		help = 'Request this number of calls to integrand from integrator'
	)
	parser.add_argument(
		'--calls-factor',
		metavar = 'F',
		type = float,
		default = 2,
		help = 'Increase number of calls by factor F when checking for convergence'
	)
	parser.add_argument(
		'--convergence-factor',
		metavar = 'O',
		type = float,
		default = 1.01,
		help = 'Consider result converged when using F times more calls gives a result within factor O.'
	)
	parser.add_argument(
		'--convergence-diff',
		metavar = 'D',
		type = float,
		default = 1e-3,
		help = 'Consider result converged when using F times more calls gives a result within difference D.'
	)
	parser.add_argument(
		'--min-value',
		metavar = 'M',
		type = float,
		default = 1e-3,
		help = 'Consider result converged when using F times more calls gives an absolute result less than M.'
	)

	args = parser.parse_args()

	if rank == 0:
		print('Arguments parsed')
		stdout.flush()

	if args.help:
		if rank == 0:
			parser.print_help()
			stdout.flush()
		exit()

	if args.dimensions < 1:
		if rank == 0:
			print('Number of dimensions must be at least 1')
		exit(1)

	if not exists(args.integrator):
		print('Integrator', args.integrator, "doesn't exist")
		exit(1)

	# prepare integrator program
	integrator = None
	if rank > 0:
		arg_list = [args.integrator]
		if args.args != None:
			arg_list += shlex.split(args.args)
		integrator = Popen(arg_list, stdin = PIPE, stdout = PIPE, universal_newlines = True, bufsize = 1)
		if args.verbose:
			print('Integrator initialized by rank', rank)
		stdout.flush()

	dimensions = list(range(args.dimensions))

	# initialize grid for integration
	grid = None
	cells_to_process = None
	if rank == 0:
		c = cell()
		c.data['id'] = 1
		c.data['processing'] = False
		c.data['converged'] = False
		c.data['value'] = None
		c.data['error'] = None
		for i in dimensions:
			c.set_extent(i, args.min_extent, args.max_extent)
		grid = ndgrid(c)
		cells_to_process = [c]

		for i in range(args.prerefine):
			split(choice(grid.get_cells()), 1, [randint(0, len(dimensions) - 1)], grid)
		if args.verbose:
			print('Grid initialized by rank', rank, 'with', len(grid.get_cells()), 'cells')
			stdout.flush()

	if rank > 0:

		while True:
			if args.verbose:
				print('Rank', rank, 'waiting for work')
				stdout.flush()

			work_item = comm.recv(source = 0, tag = 1)
			if work_item.cell_id == None:
				if args.verbose:
					print('Rank', rank, 'exiting')
					stdout.flush()
				exit()

			if args.verbose:
				print('Rank', rank, 'processing cell', work_item.cell_id)
				stdout.flush()

			integrator.stdin.write('{:.15e} '.format(args.calls))
			for extent in work_item.volume:
				integrator.stdin.write('{:.15e} {:.15e} '.format(extent[0], extent[1]))
			integrator.stdin.write('\n')
			integrator.stdin.flush()
			stdout.flush()

			value, error, split_dim = integrator.stdout.readline().strip().split()
			work_item.value, work_item.error, work_item.split_dim = float(value), float(error), int(split_dim)

			# check convergence
			integrator.stdin.write('{:.15e} '.format(args.calls * args.convergence_factor))
			for extent in work_item.volume:
				integrator.stdin.write('{:.15e} {:.15e} '.format(extent[0], extent[1]))
			integrator.stdin.write('\n')
			integrator.stdin.flush()

			new_value, new_error, new_split_dim = integrator.stdout.readline().strip().split()
			new_value, new_error, new_split_dim = float(new_value), float(new_error), int(new_split_dim)

			try:
				convg_fact = max(abs(work_item.value), abs(new_value)) / min(abs(work_item.value), abs(new_value))
			except:
				convg_fact = 0.0
			convg_diff = abs(work_item.value - new_value)
			work_item.value = new_value
			work_item.error = new_error

			if \
				convg_fact < args.convergence_factor \
				or convg_diff < args.convergence_diff \
				or abs(new_value) < args.min_value \
			:
				if args.verbose:
					print('Rank', rank, 'converged')
					stdout.flush()
				work_item.converged = True
			else:
				if args.verbose:
					print('Rank', rank, "didn't converge, returning split dimension", new_split_dim)
					stdout.flush()
				work_item.value = work_item.error = None
				work_item.split_dim = new_split_dim

			if args.verbose:
				print('Rank', rank, 'returning work')
				stdout.flush()
			comm.send(obj = work_item, dest = 0, tag = 1)

	else: # rank == 0

		work_trackers = [Work_Tracker() for i in range(comm.size - 1)]
		for work_tracker in work_trackers:
			work_tracker.processing = False
			work_tracker.item = Work_Item()

		if args.verbose:
			print('Number of work item slots:', len(work_trackers))
			stdout.flush()

		while True:
			sleep(0.1)

			work_left = 0
			for c in grid.get_cells():
				if not c.data['converged']:
					work_left += 1
			for proc in range(len(work_trackers)):

				# idle worker
				if not work_trackers[proc].processing:

					# find cell to process
					for c in grid.get_cells():
						if c.data['converged'] or c.data['processing']:
							continue

						# found
						work_trackers[proc].processing = True
						all_idle = False
						c.data['processing'] = True
						work_trackers[proc].item.converged = False
						work_trackers[proc].item.cell_id = c.data['id']
						work_trackers[proc].item.volume = [c.get_extent(dim) for dim in dimensions]
						if args.verbose:
							print('Sending cell', c.data['id'], 'for processing to rank', proc + 1)
							stdout.flush()
						comm.send(obj = work_trackers[proc].item, dest = proc + 1, tag = 1)
						break

				else:

					# check if result ready
					if comm.Iprobe(source = proc + 1, tag = 1):
						work_left -= 1
						work_trackers[proc].processing = False
						work_trackers[proc].item = comm.recv(source = proc + 1, tag = 1)
						cell_id = work_trackers[proc].item.cell_id
						if args.verbose:
							print('Received result for cell', cell_id, 'from process', proc + 1)
							stdout.flush()
						found = False
						for c in grid.get_cells():
							if c.data['id'] == cell_id:
								found = True
						if not found:
							print('Cell', cell_id, 'not in grid')
							stdout.flush()
							exit(1)
						for c in grid.get_cells():
							if c.data['id'] == cell_id:
								c.data['processing'] = False
								c.data['converged'] = work_trackers[proc].item.converged
								if work_trackers[proc].item.value == None and work_trackers[proc].item.converged:
									print('Got none value for converged result')
									stdout.flush()
									exit(1)
								c.data['value'] = work_trackers[proc].item.value
								c.data['error'] = work_trackers[proc].item.error
								split_dim = work_trackers[proc].item.split_dim
								if not c.data['converged']:
									if args.verbose:
										print("Cell didn't converge, splitting along dimension", split_dim)
										stdout.flush()
									split(c, 1, [split_dim], grid)
									work_left += 2
								break

			#print('Work left', work_left)
			#stdout.flush()
			if work_left <= 0:
				stdout.flush()
				break

		for i in range(1, comm.size):
			comm.send(obj = Work_Item(), dest = i, tag = 1)

	if rank == 0:
		value, error = 0.0, 0.0
		for c in grid.get_cells():
			value += c.data['value']
			error += c.data['error']
		print(value, error)
