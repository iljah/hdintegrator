all: \
	integrators/example2 \
	integrators/turbulence2_plain \
	integrators/turbulence2_miser \
	integrators/turbulence2_vegas \
	integrators/turbulence3_hcubature \
	integrators/turbulence3_pcubature \
	integrators/turbulence3_hcubature_v \
	integrators/turbulence3_pcubature_v

COMP = g++ -std=c++11 -O3 -march=native -W -Wall -Wextra -Wpedantic $< -o $@

integrators/example2: integrators/example2.cpp Makefile
	$(COMP) -lgsl -lgslcblas

integrators/turbulence2_plain: integrators/turbulence2.cpp Makefile
	$(COMP) -DMETHOD=1 -lgsl -lgslcblas -lboost_program_options

integrators/turbulence2_miser: integrators/turbulence2.cpp Makefile
	$(COMP) -DMETHOD=2 -lgsl -lgslcblas -lboost_program_options

integrators/turbulence2_vegas: integrators/turbulence2.cpp Makefile
	$(COMP) -DMETHOD=3 -lgsl -lgslcblas -lboost_program_options

integrators/turbulence3_hcubature: integrators/turbulence3.cpp Makefile
	$(COMP) -lhcubature -lboost_program_options

integrators/turbulence3_pcubature: integrators/turbulence3.cpp Makefile
	$(COMP) -DPCUBATURE -lpcubature -lboost_program_options

integrators/turbulence3_hcubature_v: integrators/turbulence3.cpp Makefile
	$(COMP) -DVECTORIZE -lhcubature -lboost_program_options

integrators/turbulence3_pcubature_v: integrators/turbulence3.cpp Makefile
	$(COMP) -DVECTORIZE -DPCUBATURE -lpcubature -lboost_program_options

c: clean
clean:
	rm -f integrators/example2 integrators/turbulence2_plain integrators/turbulence2_miser integrators/turbulence2_vegas integrators/turbulence3_hcubature integrators/turbulence3_pcubature integrators/turbulence3_hcubature_v integrators/turbulence3_pcubature_v
