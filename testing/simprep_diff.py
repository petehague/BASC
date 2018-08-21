#!/usr/bin/env python

from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm
from numpy.fft import fft2, ifft2
from astropy.io import fits
import os
import sys

nmaps = 100

factor = 10
scale = 0.4 

inset = 512
nchan = 1

code = sys.argv[1]

fluxval = 0.00625
contrast = 10*float(code)

realmap = np.zeros(shape=(1,nchan,inset,inset), dtype=float)
dmap = fits.open("NGC1808_cont.dirty5.fits")

np.random.seed(1234)

pointlog = open("pointlog_diff{}.txt".format(code), "w")
pointlog.write("index r theta dx dy\n")
for i in range(nmaps):
    pointmap = np.zeros(shape=(1,nchan,inset,inset), dtype=float)
    r = float(np.random.randint(factor))*factor
    theta = np.random.rand()*(2.0*np.pi)
    x = r*np.cos(theta)
    y = r*np.sin(theta)
    pointlog.write("{} {} {} {} {}\n".format(i,r,theta,x,y))
    for f in range(0,nchan):
        pointmap[0,f,int(inset/2),int(inset/2)] += fluxval*contrast 
        pointmap[0,f,int(inset/2 + x),int(inset/2 + y)] += fluxval

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

    os.system("rm -f diffmap{}_{}.fits".format(code,i))
    fits.writeto("diffmap{}_{}.fits".format(code,i), pointmap, newheader)
pointlog.close()
