CC = gcc
CXX = g++

all: pyskimage mcmc pytest
	@echo Make Complete

pytest: bin/pybasc.o bin/skimage.o bayesys/bayesys3.o bayesys/random.o bayesys/hilbert.o bayesys/app.o bin/options.o
	$(CXX) -g -std=c++11 bin/pybasc.o bin/skimage.o bayesys/bayesys3.o bayesys/random.o bayesys/hilbert.o bin/options.o bayesys/app.o -opytest

mcmc: bin/mcmc.o bin/skimage.o bayesys/bayesys3.o bayesys/random.o bayesys/hilbert.o bayesys/app.o bin/options.o
	$(CXX) -g -std=c++11 bin/mcmc.o bin/skimage.o bayesys/bayesys3.o bayesys/random.o bayesys/hilbert.o bin/options.o bayesys/app.o -omcmc

bin/mcmc.o: source/mcmc.cpp
	$(CXX) -g -std=c++11 -c source/mcmc.cpp -obin/mcmc.o

bin/skimage.o: source/skimage.cpp
	$(CXX) -g -std=c++11 -c source/skimage.cpp -obin/skimage.o

bin/options.o: source/options.cpp
	$(CXX) -g -std=c++11 -c source/options.cpp -obin/options.o

bin/pybasc.o: source/pybasc.cpp
	$(CXX) -g -std=c++11 -c -DSTANDALONE source/pybasc.cpp -obin/pybasc.o

pyskimage:
	CC='$(CC)' CXX='$(CC)' CXX_STANDARD='-std=c++11' python3 setup.py build

binfolder:
	mkdir -p bin

clean:
	rm -rf build
	rm -f bin/*
	rm -f bayesys/*.o
	rm -f mcmc
	rm -f model.pyc
	rm -f *.so
