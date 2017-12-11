# customize the following to suite your environment
# by calling e.g. make PYTHON=python3 CXX=CC, etc.
PYTHON ?= python
MPIEXEC ?= mpiexec
DIFF ?= diff
TOUCH ?= touch
CXX ?= c++
CPPFLAGS ?=
CXXFLAGS ?= -std=c++14 -O3 -march=native -W -Wall -Wextra -Wpedantic
LDFLAGS ?=
BOOST_CPPFLAGS ?=
BOOST_CXXFLAGS ?=
BOOST_LDFLAGS ?= -lboost_program_options
GSL_CPPFLAGS ?=
GSL_CXXFLAGS ?=
GSL_LDFLAGS ?= -lgsl -lgslcblas

PROGRAMS=integrands/failing \
	integrands/maybe_failing \
	integrands/hanging \
	integrands/maybe_hanging \
	integrands/N-sphere \
	integrands/burgers_plain \
	integrands/burgers_miser \
	integrands/burgers_vegas

all: $(PROGRAMS)

BOOST_FLAGS = $(BOOST_CPPFLAGS) $(BOOST_CXXFLAGS) $(BOOST_LDFLAGS)
GSL_FLAGS = $(GSL_CPPFLAGS) $(GSL_CXXFLAGS) $(GSL_LDFLAGS)
COMP = $(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) $< -o $@

integrands/failing: integrands/failing.cpp Makefile
	$(COMP)

integrands/maybe_failing: integrands/maybe_failing.cpp Makefile
	$(COMP)

integrands/hanging: integrands/hanging.cpp Makefile
	$(COMP)

integrands/maybe_hanging: integrands/maybe_hanging.cpp Makefile
	$(COMP)

integrands/N-sphere: integrands/N-sphere.cpp Makefile
	$(COMP) integrands/gsl/plain2.c -I integrands/gsl $(GSL_FLAGS)

integrands/burgers_plain: integrands/burgers.cpp integrands/gsl/plain2.c Makefile
	$(COMP) integrands/gsl/plain2.c -DMETHOD=1 -I integrands/gsl $(GSL_FLAGS) $(BOOST_FLAGS)

integrands/burgers_miser: integrands/burgers.cpp integrands/gsl/miser2.c Makefile
	$(COMP) integrands/gsl/miser2.c -DMETHOD=2 -I integrands/gsl $(GSL_FLAGS) $(BOOST_FLAGS)

integrands/burgers_vegas: integrands/burgers.cpp Makefile
	$(COMP) integrands/gsl/vegas2.c -DMETHOD=3 -I integrands/gsl $(GSL_FLAGS) $(BOOST_FLAGS)

c: clean
clean:
	rm -f $(PROGRAMS) tests/*out tests/*ok

t: test
test: tests/2d_ok tests/3d_ok

tests/2d_ok: hdintegrator.py integrands/N-sphere.py Makefile
	@printf 'TEST N-sphere.py 2d... ' && $(MPIEXEC) -n 2 ./hdintegrator.py --integrand integrands/N-sphere.py --dimensions 1 | $(PYTHON) -c "from sys import stdin; val,err,vol=stdin.read().split(); print('{:.12e} {:.4e} {:.12e}'.format(float(val),float(err),float(vol)))" > tests/2d_out
	@$(DIFF) -q tests/2d_ref tests/2d_out && $(TOUCH) tests/2d_ok && echo PASSED

tests/3d_ok: hdintegrator.py integrands/N-sphere.py Makefile
	@printf 'TEST N-sphere.py 3d... ' && $(MPIEXEC) -n 2 ./hdintegrator.py --integrand integrands/N-sphere.py --dimensions 2 | $(PYTHON) -c "from sys import stdin; val,err,vol=stdin.read().split(); print('{:.12e} {:.4e} {:.12e}'.format(float(val),float(err),float(vol)))" > tests/3d_out
	@$(DIFF) -q tests/3d_ref tests/3d_out && $(TOUCH) tests/3d_ok && echo PASSED
