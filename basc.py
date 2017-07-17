#!/usr/bin/env python

import os
from astropy.io import fits
import sys
import re
from numpy.fft import fft2, ifft2
import numpy as np
from astropy.table import Table


def imageout(arr, name):
    sh = arr.shape
    xmax = sh[0]
    ymax = sh[1]

    output = open(name, "w")
    for y in range(ymax):
        for x in range(xmax):
            output.write("{} ".format(float(arr[x, y])))
        output.write("\n")
    output.close()


if len(sys.argv) < 4:
    print("basc.py <dirty map> <dirty psf> <dirty primary beam>")
    sys.exit(1)


dirtyMap = fits.open(sys.argv[1])
dirtyBeam = fits.open(sys.argv[2])
mapsize = dirtyMap[0].header['NAXIS1']
mapdepth = dirtyMap[0].header['NAXIS3']
beamsize = dirtyBeam[0].header['NAXIS1']

filelist = []
for arg in sys.argv[1:]:
    infile = arg
    outfile = re.split("\.", infile)[:-1]
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
    nf = image[0].header['NAXIS3']
    for y in range(ny):
        for x in range(nx):
            for f in range(nf):
                output.write(" {}".format(image[0].data[0, f, y, x]))
        output.write("\n")
    output.close()
    filelist.append(outfile)

os.system(os.path.dirname(os.path.abspath(__file__)) +
          "/mcmc {} {} {} {} {}".format(filelist[0],
                                        filelist[1],
                                        filelist[2],
                                        mapsize,
                                        mapdepth))

chain = Table.read("chain.txt", format="ascii")
nmodels = open("info.txt", "r")
afactor = 1.0/float(nmodels.readline())

xsize = int(dirtyMap[0].header['NAXIS1'])
ysize = int(dirtyMap[0].header['NAXIS2'])
bxsize = int(dirtyBeam[0].header['NAXIS1'])
bysize = int(dirtyBeam[0].header['NAXIS2'])

inset = 0
if (xsize == bxsize):
    xsize = int(xsize/2)
    ysize = int(ysize/2)
    inset = bxsize/4

atoms = np.zeros(shape=(bysize, bxsize), dtype=float)
for line in chain:
    x = int(round(line[0]) + inset)
    y = int(round(line[1]) + inset)
    f = line[2]
    atoms[y, x] += f

atoms = np.multiply(atoms, afactor)
atomFT = np.fft.fftshift(fft2(atoms))

beam = np.zeros(shape=(bysize, bxsize), dtype=float)
xoff = xsize - bxsize/2
yoff = ysize - bysize/2
for x in range(0, bxsize):
    for y in range(0, bysize):
        u = int(y-yoff)
        v = int(x-xoff)
        if (u > 0 and u < bxsize and v > 0 and v < bysize):
            beam[y, x] = dirtyBeam[0].data[0, 0, v, u]/float(bxsize*bysize)

beamFT = np.fft.fftshift(fft2(beam))
result = np.absolute(ifft2(np.fft.ifftshift(np.multiply(atomFT, beamFT))))

imageout(beam, "zbeam.txt")
imageout(np.absolute(atomFT), "atomft.txt")
imageout(np.absolute(beamFT), "beamft.txt")
imageout(np.absolute(atomFT*beamFT), "totalft.txt")

imageout(result, "resid.txt")
imageout(atoms, "atoms.txt")
imageout(np.absolute(ifft2(np.fft.ifftshift(atomFT))), "atoms2.txt")

clresult = np.zeros(shape=(ysize, xsize), dtype=float)
for x in range(0, xsize):
    for y in range(0, ysize):
        clresult[y, x] = result[int(y+ysize/2), int(x+xsize/2)]

sqsum = 0.0
dmap = dirtyMap[0].data[0, 0]
for y in range(ysize):
    for x in range(xsize):
        sqsum = sqsum + (dmap[y, x]-clresult[y, x])**2
print("RMS residual = {}".format(np.sqrt(sqsum/(xsize*ysize))))

sqsum = 0.0
for y in range(ysize):
    for x in range(xsize):
        sqsum += dmap[y, x]**2
print("Base residual = {}".format(np.sqrt(sqsum/(xsize*ysize))))

sys.exit(0)
