#!/usr/bin/env python

import basc
from astropy.table import Table

basc.readConfig("config.txt")

basc.loadMap("ex_image.fits")
basc.loadBeam("ex_psf.fits")
basc.loadPBCor("ex_flux.fits")

basc.run()

result = basc.getChain()

# Add in any processing of results here

result.write("chain.txt", format="ascii")
print("Models written to chain.txt")
