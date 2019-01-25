all: pyskimage mcmc pytest
	@echo Make Complete

pytest: bin/pybasc.o bin/skimage.o bayesys/bayesys3.o bayesys/random.o bayesys/hilbert.o bayesys/app.o bin/options.o
	g++ -g -std=c++11 bin/pybasc.o bin/skimage.o bayesys/bayesys3.o bayesys/random.o bayesys/hilbert.o bin/options.o bayesys/app.o -opytest

mcmc: bin/mcmc.o bin/skimage.o bayesys/bayesys3.o bayesys/random.o bayesys/hilbert.o bayesys/app.o bin/options.o
	g++ -g -std=c++11 bin/mcmc.o bin/skimage.o bayesys/bayesys3.o bayesys/random.o bayesys/hilbert.o bin/options.o bayesys/app.o -omcmc

bin/mcmc.o: source/mcmc.cpp binfolder
	g++ -g -std=c++11 -c source/mcmc.cpp -obin/mcmc.o

bin/skimage.o: source/skimage.cpp binfolder
	g++ -g -std=c++11 -c source/skimage.cpp -obin/skimage.o

bin/options.o: source/options.cpp binfolder
	g++ -g -std=c++11 -c source/options.cpp -obin/options.o

bin/pybasc.o: source/pybasc.cpp binfolder
	g++ -g -std=c++11 -c -DSTANDALONE source/pybasc.cpp -obin/pybasc.o

pyskimage:
	CXXFLAGS='-std=c++11' python setup.py build

binfolder:
	mkdir bin

clean:
	rm -rf build
	rm -f bin/*
	rm -f bayesys/*.o
	rm -f mcmc
	rm -f model.pyc
	rm -f *.so
