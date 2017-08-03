#!/usr/bin/env python

import basc
from astropy.io import fits

dmap = fits.open("ex_image.fits")
dbeam = fits.open("ex_psf.fits")
flux = fits.open("ex_flux.fits")

basc.run(dmap, dmap, flux)

print("Models written to chain.txt")
