# Installation


## Downloading the program

You can obtain HDIntegrator and its submodules by cloning the repository with:

    git clone --recursive https://github.com/iljah/hdintegrator.git


## Prerequisites

HDIntegrator requires **Python 3**, **NetworkX** and **mpi4py**, while mpi4py requires an
implementation of the Message Passing Interface standard such as **OpenMPI**.


### Integrands

Integrands located in integrands directory each have their own requirements and
prerequisites.

Integrands implemented in Python require the **SciPy** package which is called
python3-scipy in Fedora.

Integrands implemented in C++ require a **C++14 compiler** and can be
compiled automatically with **GNU Make** by using the Makefile located in this
directory. These integrands also require one of the following libraries:
**GNU GSL**, **cubature** available at [github.com/stevengj/cubature](https://github.com/stevengj/cubature).
In Fedora the required GSL package is gsl-devel while in Ubuntu it is libgsl-dev.


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
before installing mpi4py.

Create a new python environment with virtualenv:

    virtualenv --python=python3 ~/hdintegrator_env

Switch to that environment:

    . ~/hdintegrator_env/bin/activate

and install prerequisites (also for the [N-sphere.py](integrands/N-sphere.py) integrand):

    pip install networkx mpi4py scipy

make sure that the OpenMPI wrapper compiler is in your PATH before running pip.


# Testing

Integrating the x >= 0 and y >= 0 quarter of a 2d unit sphere with SciPy-based
python program:

    mpiexec -n 2 ./hdintegrator.py --integrator integrands/N-sphere.py --dimensions 1

should give a value of 0.785, error less than 0.001 and a NaN volume of 0, for
example:

    0.7853981633974481 8.833911380179416e-11 0.0


The same for 3d unit sphere:

    mpiexec -n 2 ./hdintegrator.py --integrator integrands/N-sphere.py --dimensions 2

should give a value of 0.52, error less than 0.001 and NaN volume 0:

    0.5235987763835492 1.453148457120079e-08 0.0


Integrating a 15d unit sphere with GSL-based C++ program:

    mpiexec -n 2 ./hdintegrator.py --integrator integrands/N-sphere --dimensions 15 --prerefine 200

should give a value of 3.6e-6, error < 1e-6 and NaN volume of 0:

    3.6103624726239274e-06 2.263344140165195e-07 0.0
