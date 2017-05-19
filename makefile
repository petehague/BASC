all: pyskimage mcmc
	@echo Make Complete

mcmc: bin/mcmc.o bin/skimage.o bayesys/bayesys3.o bayesys/random.o bayesys/hilbert.o bayesys/app.o
	g++ -std=c++11 bin/mcmc.o bin/skimage.o bayesys/bayesys3.o bayesys/random.o bayesys/hilbert.o bayesys/app.o -omcmc

bin/mcmc.o: source/mcmc.cpp
	g++ -std=c++11 -c source/mcmc.cpp -obin/mcmc.o

bin/skimage.o: source/skimage.cpp
	g++ -std=c++11 -c source/skimage.cpp -obin/skimage.o

pyskimage:
	./setup.py build

clean:
	rm -rf build
	rm -f bin/*
	rm -f bayesys/*.o
	rm -f mcmc
	rm model.pyc
