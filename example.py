#!/usr/bin/env python

import basc
from astropy.table import Table

basc.readConfig("config.txt")

# Generate a view objecct to work with
newView = basc.view()

# Load in the fits files
newView.loadMap("ex_image.fits")
newView.loadBeam("ex_psf.fits")
newView.loadPBCor("ex_flux.fits")

# Run the MCMC process
newView.run()
result = newView.getChain()

# Add in any processing of results here

result.write("chain.txt", format="ascii")
print("Models written to chain.txt")

print("Sources detected:")
print(newView.clusters())
