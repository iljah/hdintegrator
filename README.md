# HDIntegrator

High-Dimensional Integrator (HDIntegrator) is a program for parallel numerical
integration of functions in arbitrary number of dimensions, see [paper.md](paper.md) for a
general overview.


# Installation

Installation and testing is detailed in the [INSTALL.md](INSTALL.md) file but at least Python 3
and an implementation of Message Passing Interface is required along with mpi4py
and the NetworkX python package. Separate programs are used for evaluating
integrals numerically and are provided in the integrands directory. These can
have their own prerequisites as detailed in [INSTALL](INSTALL.md).


# Usage

To get a list of mandatory and optional arguments supported by the program run
[hdintegrator.py](hdintegrator.py) with --help argument: `./hdintegrator.py --help`
(or `hdintegrator.py --help` if you installed it through `pip`).

In general the program must be run with at least two MPI processes and must
include a serial program to use as an integrand as well as the number of
dimensions to integrate in:

    mpiexec -n 5 ./hdintegrator.py --integrator integrands/N-sphere.py --dimensions 2

The result consists of one line with the integral's value, absolute error and
fraction of volume relative to total volume in which the integrand failed to
return a value.


# Input and output formats

To define your own integral you must write a program that will be called by
hdintegrator.py for evaluating that intergral within a given volume and give that
program as an argument to hdintegrator.py. For an example see e.g. the program in
[integrands/N-sphere.py](integrands/N-sphere.py) which evaluates the integral for an N-dimensional sphere.
Communication between hdintegrator.py and integrands is handled via standard input
and output in ASCII format. Each line given to the integrand by hdintegrator.py
consists of floating point numbers separated by spaces:

    C a0 a1 b0 b1 c0 c1 ...

where C is the number of samples to use for evaluating the integral and a0, a1,
b0, etc. represent the minimum and maximum extent of the integration volume
respectively. Note that the number of dimensions might change from one line to
another, although that is not the case currently. Each line received from the
integrand by hdintegrator.py must consist of three floating point numbers separated
by spaces in ASCII format:

    V E S

where V is the value of the integral, E is an estimate of the absolute
integration error and S is the suggested dimension starting from 0 in which to
split the integration volume in order to minimize subsequent integration errors.


# Support

To seek support or report an issue in HDIntegrator please create a new issue at
https://github.com/iljah/hdintegrator/issues


# Contributing

To contribute to HDIntegrator please create a new pull request at
https://github.com/iljah/hdintegrator/pulls
