CC = gcc
CXX = g++
CPPFLAGS = -std=c++11 

all: pyskimage mcmc pytest
	@echo Make Complete

pytest: bin/pybasc.o bin/skimage.o bayesys/bayesys3.o bayesys/random.o bayesys/hilbert.o bayesys/app.o bin/options.o
	$(CXX) -g -std=c++11 bin/pybasc.o bin/skimage.o bayesys/bayesys3.o bayesys/random.o bayesys/hilbert.o bin/options.o bayesys/app.o -opytest

mcmc: bin/mcmc.o bin/skimage.o bayesys/bayesys3.o bayesys/random.o bayesys/hilbert.o bayesys/app.o bin/options.o
	$(CXX) -g -std=c++11 bin/mcmc.o bin/skimage.o bayesys/bayesys3.o bayesys/random.o bayesys/hilbert.o bin/options.o bayesys/app.o -omcmc

bin/%.o: source/%.cpp binfolder
	$(CXX) $(CPPFLAGS) $< -c -o$@

bayesys/%.o: bayesys/%.c
	$(CC) $< -c -o$@

pyskimage:
	CC='$(CC)' CXX='$(CXX)' CXX_STANDARD='-std=c++11' python3 setup.py build

binfolder:
	mkdir -p bin

clean:
	rm -rf build
	rm -f bin/*
	rm -f bayesys/*.o
	rm -f mcmc
	rm -f model.pyc
	rm -f *.so
