# README #

BaSC is short for Bayesian Source Characterisation. It is an MCMC process that performs source detection and characterisation on dirty maps, taking into account the properties of the beam in a more rigorous way that CLEAN does. 

It is based on a method developed by Steve Gull, and uses the BayeSys MCMC driver.

### Installation ###

BaSC requires a c compiler, Python 3, and the numpy and astropy libraries. On Linux, type `make` to compile the extension, or if using Mac, type 'make mac'. Currently Windows is not supported - you will have to run it in an Ubuntu terminal on a Windows machine. 

### Usage ###

To try out BaSC, run the example script `example.py` in the BaSC folder. This should locate a single source in the centre of the map, and return the models from the burned in chain to `chain.txt`

To run BaSC on your own data, use

```
basc.py <dirty map file> <dirty psf file> <primary beam flux file>
```

In general, import basc into your Python program and use it as shown in `example.py`

### Contact ###

For more information, please email Peter Hague at prh44 AT cam.ac.uk
