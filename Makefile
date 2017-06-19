all: \
	integrators/example2 \
	integrators/turbulence2_plain \
	integrators/turbulence2_miser \
	integrators/turbulence2_vegas \
	integrators/turbulence3_hcubature \
	integrators/turbulence3_pcubature \
	integrators/turbulence3_hcubature_v \
	integrators/turbulence3_pcubature_v \
	integrators/turbulence4

COMP = g++ -std=c++17 -O3 -march=native -W -Wall -Wextra -Wpedantic $< -o $@

integrators/example2: integrators/example2.cpp Makefile
	$(COMP) -lgsl -lgslcblas

integrators/turbulence2_plain: integrators/turbulence2.cpp integrators/gsl/plain2.c Makefile
	$(COMP) integrators/gsl/plain2.c -DMETHOD=1 -I integrators/gsl -lgsl -lgslcblas -lboost_program_options

integrators/turbulence2_miser: integrators/turbulence2.cpp integrators/gsl/miser2.c Makefile
	$(COMP) integrators/gsl/miser2.c -DMETHOD=2 -I integrators/gsl -lgsl -lgslcblas -lboost_program_options

integrators/turbulence2_vegas: integrators/turbulence2.cpp Makefile
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
	rm -f integrators/example2 integrators/turbulence2_plain integrators/turbulence2_miser integrators/turbulence2_vegas integrators/turbulence3_hcubature integrators/turbulence3_pcubature integrators/turbulence3_hcubature_v integrators/turbulence3_pcubature_v integrators/turbulence4
