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
from datetime import datetime, timedelta
from math import isnan
from os import rename
from os.path import dirname, exists, join, realpath
from pickle import dump, load
from random import choice, randint
import shlex
from subprocess import Popen, PIPE
from sys import path, stdout
from time import sleep

path.append(join(dirname(realpath(__file__)), 'submodules/ndgrid/source'))
try:
	from cell import cell
	from ndgrid import ndgrid
except Exception as e:
	print("Couldn't import required submodule(s): ", e)
	print("Did you remember to clone them along with hdintegrator by using:")
	print('git clone --recursive https://github.com/iljah/hdintegrator.git')
	exit(1)

try:
	from mpi4py import MPI
except Exception as e:
	exit("Couldn't import mpi4py: " + str(e))


'''
Splits a cell.

\param cell Cell to split
\param splits Number of times to split given cell, it's children, etc.
\param dimensions List of dimensions in which to split
\param grid Grid in which to split given cell.

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


'''
Used for transferring work between rank 0 and other ranks.

\var volume List of pairs indicating minimum and maximum extent of integration volume in each dimension
\var cell_id Unique id of a grid cell
\var converged Whether result of integration converged
\var value Value of integral
\var error Estimate of absolute error for calculated integral
\var split_dim Suggested dimension for splitting the volume in case result didn't converge
'''
class Work_Item:
	def __init__(self):
		self.volume = None
		self.cell_id = None
		self.converged = None
		self.value = None
		self.error = None
		self.split_dim = None

	def __str__(self):
		ret_val = 'Id: ' + str(self.cell_id) + ', Vol: '
		for extent in self.volume:
			ret_val += '[' + str(extent[0]) + ', ' + str(extent[1]) + '], '
		return ret_val

	__repr__ = __str__


'''
Used by rank 0 to keep track of worker ranks.
'''
class Work_Tracker:
	def __init__(self):
		self.item = None
		self.processing = None
		self.start_time = None


'''
Returns basic info about the calculated solution.

\param grid Integration grid.

\return Tuple with current integral's value, error, NaN volume, total volume, number of converged cells and total number of grid cells.
'''
def get_info(grid):
	converged_cells = 0
	# sum up final result
	total_vol, nan_vol = 0.0, 0.0
	value, error = 0.0, 0.0

	cells = grid.get_cells()
	for c in cells:
		if c.data['converged']:
			converged_cells += 1

		vol = 1
		extents = c.get_extents()
		for extent in extents:
			vol *= extents[extent][1] - extents[extent][0]
		total_vol += vol

		if c.data['value'] != None:
			if isnan(c.data['value']):
				nan_vol += vol
			else:
				value += c.data['value']

		if c.data['error'] != None and not isnan(c.data['error']):
			error += c.data['error']

	return value, error, nan_vol, total_vol, converged_cells, len(cells)


'''
Prepares an integrand with Popen.

\param args Result from parse_args() of argparse.ArgumentParser in __main__.

\return Value returned by Popen.
'''
def prepare_integrand(args):
	arg_list = [args.integrand]
	if args.args != None:
		arg_list += shlex.split(args.args)
	integrand = Popen(arg_list, stdin = PIPE, stdout = PIPE, universal_newlines = True, bufsize = 1)
	if args.verbose:
		print('Integrand initialized by rank', rank)
	stdout.flush()
	return integrand


if __name__ == '__main__':

	comm = MPI.COMM_WORLD
	rank = comm.Get_rank()

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
		'--integrand',
		default = '',
		help = 'Path to integrand program to use, relative to current working directory'
	)
	parser.add_argument(
		'--dimensions',
		type = int,
		default = 0,
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
		help = 'Arguments to pass to integrand program, given as one string (e.g. --args "-a b -c d") which are passed on to integrand after splitting with shlex.split'
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
		help = 'Request this number of calls to integrand'
	)
	parser.add_argument(
		'--timer',
		type = int,
		default = 9999,
		metavar = 'T',
		help = 'Consider workers that do not return a result within T seconds as failed'
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
		help = 'Consider result converged when using F times more calls gives a result within factor O'
	)
	parser.add_argument(
		'--convergence-diff',
		metavar = 'D',
		type = float,
		default = 1e-3,
		help = 'Consider result converged when using F times more calls gives a result within difference D'
	)
	parser.add_argument(
		'--min-value',
		metavar = 'M',
		type = float,
		default = 1e-3,
		help = 'Consider result converged when using F times more calls gives an absolute result less than M'
	)
	parser.add_argument(
		'--restart',
		metavar = 'R',
		default = '',
		help = 'If R exists continue integration from result in R, write result to R every I seconds (do not continue from result files of untrusted sources)'
	)
	parser.add_argument(
		'--restart-interval',
		metavar = 'I',
		type = int,
		default = -1,
		help = 'If I > 0 write result to file R every I seconds during integration'
	)
	parser.add_argument(
		'--inspect',
		default = '',
		help = 'If not empty, print information about given restart file and exit'
	)

	args = parser.parse_args()

	if args.help:
		if rank == 0:
			parser.print_help()
			stdout.flush()
		exit()

	if args.inspect != '':
		if rank == 0:
			if not exists(args.inspect):
				print('Restart file', args.inspect, "doesn't exist")
				exit(1)
			with open(args.inspect, 'rb') as restartfile:
				grid = load(restartfile)
			value, error, nan_vol, total_vol, converged, nr_cells = get_info(grid)
			print('Value:', value, 'error:', error, 'NaN volume/total:', nan_vol / total_vol, ',', converged, '/', nr_cells, 'converged cells')
		exit()

	if comm.size < 2:
		if rank == 0:
			print('At least 2 processes required')
		exit(1)

	if rank == 0 and args.verbose:
		print('Starting with', comm.size, 'processes')
		stdout.flush()

	if args.dimensions < 1:
		if rank == 0:
			print('Number of dimensions must be at least 1')
		exit(1)

	if not exists(args.integrand):
		print('Integrand', args.integrand, "doesn't exist")
		exit(1)


	dimensions = list(range(args.dimensions))

	if rank == 0:

		# prepare grid for integration
		grid = None
		restart = False

		if args.restart != '' and exists(args.restart):
			restart = True

		if restart:
			if args.verbose:
				print('Restarting from', args.restart, end = '...  ')
			with open(args.restart, 'rb') as restartfile:
				grid = load(restartfile)

			converged = 0
			for c in grid.get_cells():
				if c.data['converged']:
					converged += 1
				c.data['processing'] = False
			print(converged, '/', len(grid.get_cells()), 'converged')

		else:
			c = cell()
			c.data['id'] = 1
			c.data['processing'] = False
			c.data['converged'] = False
			c.data['value'] = None
			c.data['error'] = None
			for i in dimensions:
				c.set_extent(i, args.min_extent, args.max_extent)
			grid = ndgrid(c)

			for i in range(args.prerefine):
				split(choice(grid.get_cells()), 1, [randint(0, len(dimensions) - 1)], grid)
			if args.verbose:
				print('Grid initialized by rank', rank, 'with', len(grid.get_cells()), 'cells')
				stdout.flush()

		# tracker for every rank > 0
		work_trackers = [Work_Tracker() for i in range(comm.size - 1)]
		for work_tracker in work_trackers:
			work_tracker.processing = False
			work_tracker.item = Work_Item()

		if args.verbose:
			print('Number of work item slots:', len(work_trackers))
			stdout.flush()

		next_restart = datetime.now() + timedelta(seconds = args.restart_interval)
		while True:
			sleep(0.1)

			# write restart if needed
			now = datetime.now()
			if args.restart_interval > 0 and next_restart <= now:
				if args.verbose:
					print('Writing restart file at', now.isoformat().split('.')[0])
				next_restart += timedelta(seconds = args.restart_interval)

				if exists(args.restart):
					rename(args.restart, args.restart + '-' + now.isoformat().split('.')[0])
				with open(args.restart, 'wb') as restartfile:
					dump(grid, restartfile)


			work_left = 0
			processing = 0
			for c in grid.get_cells():
				if not c.data['converged']:
					work_left += 1
				if c.data['processing']:
					processing += 1
			if args.verbose:
				print(work_left, 'work left,', processing, 'processing')

			nr_failed = 0
			for proc in range(len(work_trackers)):

				# failed worker
				if work_trackers[proc].processing == None:
					nr_failed += 1
					continue

				# idle worker
				if not work_trackers[proc].processing:

					# find cell to process
					for c in grid.get_cells():
						if c.data['converged'] or c.data['processing']:
							if c.data['processing']:
								processing += 1
							continue

						# found
						found = True
						work_trackers[proc].processing = True
						c.data['processing'] = True
						work_trackers[proc].item.converged = False
						work_trackers[proc].item.cell_id = c.data['id']
						work_trackers[proc].item.volume = [c.get_extent(dim) for dim in dimensions]
						if args.verbose:
							print('Sending cell', c.data['id'], 'for processing to rank', proc + 1)
							stdout.flush()
						comm.send(obj = work_trackers[proc].item, dest = proc + 1, tag = 1)
						work_trackers[proc].start_time = datetime.now()
						break

				else:

					# if result ready
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
								break
						if not found:
							print('Cell', cell_id, 'not in grid')
							stdout.flush()
							exit(1)
						for c in grid.get_cells():
							if c.data['id'] == cell_id:
								c.data['processing'] = False
								c.data['converged'] = work_trackers[proc].item.converged
								if work_trackers[proc].item.value == None and work_trackers[proc].item.converged:
									print('Worker', proc + 1, 'failed')
									stdout.flush()
									work_trackers[proc].processing = None
									break
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

					# if result not ready
					else:
						processing_time = (datetime.now() - work_trackers[proc].start_time).seconds
						if processing_time > args.timer:
							print('Marking rank', proc + 1, 'as failed due to exceeded processing time, work item', work_trackers[proc].item)
							work_trackers[proc].processing = None
							break


			if work_left <= 0:
				stdout.flush()
				break

			if nr_failed >= comm.size - 1:
				print('All workers failed, exiting...')
				stdout.flush()
				break

		# tell others to quit
		for i in range(1, comm.size):
			comm.send(obj = Work_Item(), dest = i, tag = 1)

		value, error, nan_vol, total_vol, converged, nr_cells = get_info(grid)
		print(value, error, nan_vol / total_vol)


	else: # if rank == 0

		integrand = prepare_integrand(args)

		# work loop
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

			work_item.value = float('NaN')
			work_item.error = float('NaN')
			work_item.converged = False

			failed = False
			to_stdin = '{:.16e} '.format(args.calls)
			for extent in work_item.volume:
				ext_str = '{:.16e} {:.16e} '.format(extent[0], extent[1])
				first, second = ext_str.split()
				if first == second or float(first) >= float(second):
					failed = True
					break
				to_stdin += ext_str
			if failed:
				print('Rank', rank, 'invalid extent, returning NaN')
				comm.send(obj = work_item, dest = 0, tag = 1)
				continue

			try:
				integrand.stdin.write(to_stdin + '\n')
				integrand.stdin.flush()
			except Exception as e:
				print('Rank', rank, 'request to integrand failed with input', to_stdin, ', error:', e)
				integrand = prepare_integrand(args)
				comm.send(obj = work_item, dest = 0, tag = 1)
				continue

			stdout.flush()

			try:
				answer = integrand.stdout.readline()
				value, error, split_dim = answer.strip().split()
				work_item.value, work_item.error, work_item.split_dim = float(value), float(error), int(split_dim)
			except Exception as e:
				print('Rank', rank, 'call to integrand failed with result:', answer, ', returning NaN, input string:', to_stdin, ', exception:', e)
				comm.send(obj = work_item, dest = 0, tag = 1)
				continue

			# check convergence
			failed = False
			time_start = datetime.now()
			to_stdin = '{:.16e} '.format(args.calls * args.calls_factor)
			for extent in work_item.volume:
				ext_str = '{:.16e} {:.16e} '.format(extent[0], extent[1])
				first, second = ext_str.split()
				if first == second or float(first) >= float(second):
					failed = True
					break
				to_stdin += ext_str
			if failed:
				print('Rank', rank, 'invalid extent, returning NaN')
				comm.send(obj = work_item, dest = 0, tag = 1)
				continue

			try:
				integrand.stdin.write(to_stdin + '\n')
				integrand.stdin.flush()
			except Exception as e:
				print('Rank', rank, 'request to integrand failed with input', to_stdin, ', error:', e)
				integrand = prepare_integrand(args)
				comm.send(obj = work_item, dest = 0, tag = 1)
				continue

			try:
				answer = integrand.stdout.readline()
				new_value, new_error, new_split_dim = answer.strip().split()
				new_value, new_error, new_split_dim = float(new_value), float(new_error), int(new_split_dim)
			except Exception as e:
				print('Rank', rank, 'call to integrand failed with result:', answer, ', returning NaN, input string:', to_stdin, ', exception:', e)
				comm.send(obj = work_item, dest = 0, tag = 1)
				continue

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
