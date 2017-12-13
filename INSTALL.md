# Installation


## Downloading the program

Easiest way to obtain HDIntegrator is with [pip](https://pip.pypa.io):

    `pip install hdintegrator`

after setting up the environment as described in section *Local installation*. Due
to technical reasons only integrands implemented in Python will be available.
The above command will install necessary prerequisites for the main program, and
the Python integrands, under your home directory.

To use integrands implemented in C++, for higher performance or different algorithms,
you'll have to clone the repository:

    git clone --recursive https://github.com/iljah/hdintegrator.git

after which the the main program will be in `hdintegrator/hdintegrator.py` and
all the integrands in `hdintegrator/integrands/`. See the Integrands section
for details.


## Prerequisites

HDIntegrator requires **Python 3**, **NetworkX** and **mpi4py**, while mpi4py requires an
implementation of the Message Passing Interface standard such as **OpenMPI**.


## System-wide installation

Prerequisites are available for example in Fedora and Ubuntu repositories and
can be installed by following the respective distributions' documentation on
installing additional software. Note that there might be separate packages for
different versions of Python and/or the MPI implementation so be sure to install
packages corresponding to the same combinations of Python/MPI. For example the
correct NetworkX package is called python3-networkx in both Ubuntu and Fedora
but the mpi4py package is called python3-mpi4py in Ubuntu and in Fedora either
python3-mpi4py-mpich or python3-mpi4py-openmpi depending on which OpenMPI
implementation you are using.


## Local installation

You can install the required Python libraries under your home directory using
[virtualenv](https://virtualenv.pypa.io) while OpenMPI can be obtained from
[https://www.open-mpi.org](https://www.open-mpi.org). You should install OpenMPI
before installing HDIntegrator.

Create a new python environment with virtualenv:

    virtualenv --python=python3 ~/hdintegrator_env

switch to that environment:

    . ~/hdintegrator_env/bin/activate

and install HDIntegrator and its prerequisites:

    pip install hdintegrator

make sure that the OpenMPI wrapper compiler is in your PATH before running `pip`
so that the `mpiexec` binary from OpenMPI will be found. The main program and
Python integrands will be located in `~/hdintegrator_env/bin` and can be used
for example as:

    mpiexec -n 2 ~/hdintegrator_env/bin/hdintegrator.py --integrand ~/hdintegrator_env/bin/N-sphere.py --dimensions 1

You can leave the virtualenv environment with `deactivate`.


## Integrands

Integrands located in `integrands` directory in the git repository each have
their own requirements and prerequisites.

Integrands implemented in Python require the **SciPy** package which is called
python3-scipy in Fedora and is installed automatically if you used `pip` to
install HDIntegrator as described in section *Local installation*.

Integrands implemented in C++ require a **C++14 compiler** and can be
compiled automatically with **GNU Make** by using the Makefile located in this
directory. These integrands also require [GSL](https://www.gnu.org/software/gsl/)
and the *program options* part of [boost](http://www.boost.org/) libraries
which can be installed in Fedora by installing *gsl-devel* and *boost-devel*
and in Ubuntu by installing *libgsl-dev* and *libboost-program-options-dev*.
Compilation parameters can be customized by invoking `make` with the desired
parameter and value separated by `=`:

    make CXXFLAGS=-std=c++14

See top of `Makefile` for list of parameters used.


# Testing

The first two tests below can be executed by running `make test`. You can
specify a different executable for running the integrator in parallel by giving
the MPIEXEC argument to make, for example `make test MPIEXEC=/path/to/mpirun`.

Integrating the x >= 0 and y >= 0 quarter of a 2d unit sphere with SciPy-based
python program:

    mpiexec -n 2 ./hdintegrator.py --integrand integrands/N-sphere.py --dimensions 1

should give a value of 0.785, error less than 0.001 and a NaN volume of 0, for
example:

    0.7853981633974481 8.833911380179416e-11 0.0


The same for 3d unit sphere:

    mpiexec -n 2 ./hdintegrator.py --integrand integrands/N-sphere.py --dimensions 2

should give a value of 0.52, error less than 0.001 and NaN volume 0:

    0.5235987763835492 1.453148457120079e-08 0.0


Integrating a 15d unit sphere with GSL-based C++ program:

    mpiexec -n 2 ./hdintegrator.py --integrand integrands/N-sphere --dimensions 15 --prerefine 200

should give a value of 3.6e-6, error < 1e-6 and NaN volume of 0:

    3.6103624726239274e-06 2.263344140165195e-07 0.0
