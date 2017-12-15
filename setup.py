#! /usr/bin/env python3

from setuptools import setup

setup(
	name = 'hdintegrator',
	version = '1.0.1',
	classifiers = [
		'Development Status :: 5 - Production/Stable',
		'Environment :: Console',
		'Intended Audience :: Science/Research',
		'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
		'Natural Language :: English',
		'Operating System :: Unix',
		'Programming Language :: Python :: 3',
		'Topic :: Scientific/Engineering :: Mathematics'
	],
	python_requires = '>=3',
	install_requires = ['networkx', 'mpi4py', 'scipy'],
	scripts = [
		'hdintegrator.py',
		'submodules/ndgrid/source/cell.py',
		'submodules/ndgrid/source/ndgrid.py',
		'integrands/N-sphere.py'
	],
	author = 'Ilja Honkonen',
	author_email = 'ilja.honkonen@fmi.fi',
	url = 'https://github.com/iljah/hdintegrator',
	description = 'High-Dimensional Integrator, see github.com/iljah/hdintegrator for full version'
)
