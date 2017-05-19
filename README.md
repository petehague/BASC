# README #

Basc is short for Bayesian Source Characterisation. It is an MCMC process that performs source detection and characterisation on dirty maps, taking into account the properties of the beam in a more rigorous way that CLEAN does. 

It is based on a method developed by Steve Gull, and uses the BayeSys MCMC driver.


### Usage ###

Type make to compile the C++ portion (this will automatically compile the relevant parts of Bayesis as well) and then invoke the Python script with

basc.py <dirty map file> <dirty psf file> <primary beam flux file>

The Python script requires astropy and numpy in order to function. In addition, there is a script use_pymc.py which for comparison applies the same likelihood model using PyMC (which doesn't have all the functionality of BayeSys unfortunately)

### Ongoing Work ###

The code will be more customisable soon, and will also include more facilities for post-processing.

### Contact ###

For more information, please email Peter Hague at prh44 AT cam.ac.uk