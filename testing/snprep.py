#!/usr/bin/env python

from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm
from numpy.fft import fft2, ifft2
from astropy.io import fits
import os

nmaps = 5 

factor = 10
scale = 0.4 

inset = 512
nchan = 1

fluxval = 0.1

realmap = np.zeros(shape=(1,nchan,inset,inset), dtype=float)
dmap = fits.open("NGC1808_cont.dirty5.fits")

np.random.seed(1234)

pointlog = open("pointlog.txt", "w")
pointlog.write("index r theta\n")
for i in range(nmaps):
    pointmap = np.zeros(shape=(1,nchan,inset,inset), dtype=float)
    for f in range(0,nchan):
        pointmap[0,f,int(inset/2),int(inset/2)] += fluxval/(2**i) 

    # WARNING! WILL PRODUCE NEGATIVE CRPIX VALUES IF REFERENCE PIXEL ISNT INSIDE INSET AREA!
    newheader = dmap[0].header.copy()
    newheader['CRPIX1'] -= 0.5*(newheader['NAXIS1']-inset)
    newheader['CRPIX2'] -= 0.5*(newheader['NAXIS2']-inset)
    newheader['CDELT1'] *= scale
    newheader['CDELT2'] *= scale
    newheader['NAXIS1'] = inset
    newheader['NAXIS2'] = inset
    newheader['NAXIS3'] = nchan
    newheader['CDELT3'] = 15625000 #15 MHz
    newheader['CRPIX3'] = min(nchan,64)

    os.system("rm -f snmap_{}.fits".format(i))
    fits.writeto("snmap_{}.fits".format(i), pointmap, newheader)

