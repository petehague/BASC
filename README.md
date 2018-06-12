# README #

Basc is short for Bayesian Source Characterisation. It is an MCMC process that performs source detection and characterisation on dirty maps, taking into account the properties of the beam in a more rigorous way that CLEAN does. 

It is based on a method developed by Steve Gull, and uses the BayeSys MCMC driver.

### Usage ###

Basc requires a c compiler and Python 3. On Linux/MacOS, type `make` to compile the library and BayeSys and then test the program script with

```
basc.py <dirty map file> <dirty psf file> <primary beam flux file>
```

In general, import basc into your Python program and use it as shown in example.py

### Contact ###

For more information, please email Peter Hague at prh44 AT cam.ac.uk
