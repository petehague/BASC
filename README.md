# README #

Basc is short for Bayesian Source Characterisation. It is an MCMC process that performs source detection and characterisation on dirty maps, taking into account the properties of the beam in a more rigorous way that CLEAN does. 

It is based on a method developed by Steve Gull, and uses the BayeSys MCMC driver.


### Usage ###

Type make to compile the C++ portion (this will automatically compile the relevant parts of Bayesis as well) and then invoke the Python script with

basc.py <dirty map file> <dirty psf file> <primary beam flux file>

The Python script requires astropy and numpy in order to function. In addition, there is a script use_pymc.py which for comparison applies the same likelihood model using PyMC (which doesn't have all the functionality of BayeSys unfortunately)

### Python Module ###

There is a python module, basc2.py, which is being developed alongside the command line utility. In the long term this will
be the only way of using basc. It is still experimental at this time.

### Contact ###

For more information, please email Peter Hague at prh44 AT cam.ac.uk