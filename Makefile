all: integrators/example2

integrators/example2: integrators/example2.cpp Makefile
	g++ -std=c++11 -O3 -march=native -W -Wall -Wpedantic integrators/example2.cpp -o integrators/example2 -lgsl -lgslcblas

c: clean
clean:
	rm -f integrators/example2
