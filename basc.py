#!/usr/bin/env python

import os
from astropy.io import fits
import sys
import re
from numpy.fft import fft2,ifft2
import numpy as np
from astropy.table import Table

def imageout(arr, name):
    sh = arr.shape
    xmax = sh[0]
    ymax = sh[1]

    output = open(name, "w")
    for y in range(ymax):
        for x in range(xmax):
            output.write("{} ".format(float(arr[x,y])))
        output.write("\n")
    output.close()


if len(sys.argv)<4:
    print("sf <dirty map> <dirty psf> <dirty primary beam>")
    sys.exit(1)

image=fits.open(sys.argv[1])
mapsize = image[0].header['NAXIS1']

filelist = []
for arg in sys.argv[1:]:
    infile = arg
    outfile = re.split("\.",infile)[:-1]
    outfile.append("txt")
    outfile = '.'.join(outfile)
    if os.path.exists(outfile):
        filelist.append(outfile)
        continue

    print("Extracting "+infile)
    image = fits.open(infile)
    output = open(outfile, "w")
    nx = image[0].header['NAXIS1']
    ny = image[0].header['NAXIS2']
    for y in range(ny):
      for x in range(nx):
        output.write(" {}".format(image[0].data[0,0,y,x]))
      output.write("\n")
    output.close()
    filelist.append(outfile)

os.system(os.path.dirname(os.path.abspath(__file__))+"/mcmc {} {} {} {}".format(filelist[0],filelist[1],filelist[2],mapsize))

dirtyMap = fits.open(sys.argv[1])
dirtyBeam = fits.open(sys.argv[2])
chain = Table.read("chain.txt", format="ascii")
afactor = 1.0/float(len(chain))

xsize = dirtyMap[0].header['NAXIS1']
ysize = dirtyMap[0].header['NAXIS2']

atoms = np.zeros(shape=(ysize,xsize), dtype=float)

for line in chain:
    x = int(line[0]*xsize)
    y = int(line[1]*ysize)
    f = line[2]
    atoms[y,x] += f*afactor

atomFT = fft2(atoms)
dmap = dirtyMap[0].data[0,0]
dmapFT = fft2(np.nan_to_num(dmap))

totalFT = dmapFT - atomFT

result = ifft2(totalFT)

imageout(result.real, "resid.txt")
imageout(atoms, "atoms.txt")
imageout(atomFT, "atomsft.txt")
imageout(dmapFT, "dmapft.txt")

sqsum = 0.0
for y in range(ysize):
    for x in range(xsize):
        sqsum = sqsum + (dmap[x,y]-result.real[x,y])**2
print("RMS residual = {}".format(np.sqrt(sqsum/(xsize*ysize))))

sys.exit(0)
