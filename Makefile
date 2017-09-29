all: \
	integrators/failing \
	integrators/hanging \
	integrators/N-sphere \
	integrators/burgers_plain \
	integrators/burgers_miser \
	integrators/burgers_vegas \
	integrators/turbulence3_hcubature \
	integrators/turbulence3_pcubature \
	integrators/turbulence3_hcubature_v \
	integrators/turbulence3_pcubature_v \
	integrators/turbulence4

COMP = g++ -std=c++17 -O3 -march=native -W -Wall -Wextra -Wpedantic $< -o $@

integrators/failing: integrators/failing.cpp Makefile
	$(COMP)

integrators/hanging: integrators/hanging.cpp Makefile
	$(COMP)

integrators/N-sphere: integrators/N-sphere.cpp Makefile
	$(COMP) integrators/gsl/plain2.c -I integrators/gsl -lgsl -lgslcblas

integrators/burgers_plain: integrators/burgers.cpp integrators/gsl/plain2.c Makefile
	$(COMP) integrators/gsl/plain2.c -DMETHOD=1 -I integrators/gsl -lgsl -lgslcblas -lboost_program_options

integrators/burgers_miser: integrators/burgers.cpp integrators/gsl/miser2.c Makefile
	$(COMP) integrators/gsl/miser2.c -DMETHOD=2 -I integrators/gsl -lgsl -lgslcblas -lboost_program_options

integrators/burgers_vegas: integrators/burgers.cpp Makefile
	$(COMP) integrators/gsl/vegas2.c -DMETHOD=3 -I integrators/gsl -lgsl -lgslcblas -lboost_program_options

integrators/turbulence3_hcubature: integrators/turbulence3.cpp Makefile
	$(COMP) -lhcubature -lboost_program_options

integrators/turbulence3_pcubature: integrators/turbulence3.cpp Makefile
	$(COMP) -DPCUBATURE -lpcubature -lboost_program_options

integrators/turbulence3_hcubature_v: integrators/turbulence3.cpp Makefile
	$(COMP) -DVECTORIZE -lhcubature -lboost_program_options

integrators/turbulence3_pcubature_v: integrators/turbulence3.cpp Makefile
	$(COMP) -DVECTORIZE -DPCUBATURE -lpcubature -lboost_program_options

integrators/turbulence4: integrators/turbulence4.cpp Makefile
	$(COMP) -I $(HOME)/ohjelmat/cxx-prettyprint -lboost_program_options

c: clean
clean:
	rm -f integrators/example2 integrators/burgers_plain integrators/burgers_miser integrators/burgers_vegas integrators/turbulence3_hcubature integrators/turbulence3_pcubature integrators/turbulence3_hcubature_v integrators/turbulence3_pcubature_v integrators/turbulence4
