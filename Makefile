all: integrators/example2 integrators/turbulence2_plain integrators/turbulence2_miser integrators/turbulence2_vegas

integrators/example2: integrators/example2.cpp Makefile
	g++ -std=c++11 -O3 -march=native -W -Wall -Wpedantic $< -o $@ -lgsl -lgslcblas

integrators/turbulence2_plain: integrators/turbulence2.cpp Makefile
	g++ -DMETHOD=1 -std=c++11 -O3 -march=native -W -Wall -Wpedantic $< -o $@ -lgsl -lgslcblas -lboost_program_options

integrators/turbulence2_miser: integrators/turbulence2.cpp Makefile
	g++ -DMETHOD=2 -std=c++11 -O3 -march=native -W -Wall -Wpedantic $< -o $@ -lgsl -lgslcblas -lboost_program_options

integrators/turbulence2_vegas: integrators/turbulence2.cpp Makefile
	g++ -DMETHOD=3 -std=c++11 -O3 -march=native -W -Wall -Wpedantic $< -o $@ -lgsl -lgslcblas -lboost_program_options

c: clean
clean:
	rm -f integrators/example2 integrators/turbulence2_plain integrators/turbulence2_miser integrators/turbulence2_vegas
