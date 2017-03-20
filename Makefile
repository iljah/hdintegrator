all: integrators/example2 integrators/turbulence2_plain integrators/turbulence2_miser integrators/turbulence2_vegas

COMP = g++ -std=c++11 -O3 -march=native -W -Wall -Wextra -Wpedantic $< -o $@ -lgsl -lgslcblas

integrators/example2: integrators/example2.cpp Makefile
	$(COMP)

integrators/turbulence2_plain: integrators/turbulence2.cpp Makefile
	$(COMP) -DMETHOD=1 -lboost_program_options

integrators/turbulence2_miser: integrators/turbulence2.cpp Makefile
	$(COMP) -DMETHOD=2 -lboost_program_options

integrators/turbulence2_vegas: integrators/turbulence2.cpp Makefile
	$(COMP) -DMETHOD=3 -lboost_program_options

c: clean
clean:
	rm -f integrators/example2 integrators/turbulence2_plain integrators/turbulence2_miser integrators/turbulence2_vegas
