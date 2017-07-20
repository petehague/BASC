#!/usr/bin/env python

import os
from astropy.io import fits
import sys
import re
from numpy.fft import fft2, ifft2
import numpy as np
from astropy.table import Table
import matplotlib.pyplot as plt

# Global Variables
chain = Table()
nmodels = 0
mapsize = 0
mapdepth = 0


def run(dmap, dbeam, dflux):
    global chain, nmodels, mapsize, mapdepth

    dirtyMap = fits.open(dmap)
    mapsize = dirtyMap[0].header['NAXIS1']
    mapdepth = dirtyMap[0].header['NAXIS3']

    filelist = []
    for infile in [dmap, dbeam, dflux]:
        outfile = re.split("\.", infile)[:-1]
        outfile.append("txt")
        outfile = '.'.join(outfile)
        if os.path.exists(outfile):
            filelist.append(outfile)
            continue

        print("Extracting "+infile)
        image = fits.open(infile)
        output = open(outfile, "w")
        for y in range(mapsize):
            for x in range(mapsize):
                for f in range(mapdepth):
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


if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("basc.py <dirty map> <dirty psf> <dirty primary beam>")
        sys.exit(1)

    run(sys.argv[1], sys.argv[2], sys.argv[3])

    afactor = 1.0/float(nmodels.readline())

    insize = int(mapsize/2)
    inset = mapsize/4

    atoms = np.zeros(shape=(mapsize, mapsize), dtype=float)
    for line in chain:
        x = int(round(line[0]) + inset)
        y = int(round(line[1]) + inset)
        f = line[2]
        atoms[y, x] += f

    atoms = np.multiply(atoms, afactor)
    atomFT = np.fft.fftshift(fft2(atoms))

    dirtyMap = fits.open(sys.argv[1])
    dirtyBeam = fits.open(sys.argv[2])

    beam = np.zeros(shape=(mapsize, mapsize), dtype=float)
    offset = mapsize/2
    beamfactor = float(mapsize*mapsize)
    for x in range(0, mapsize):
        for y in range(0, mapsize):
            u = int(y-offset)
            v = int(x-offset)
            if (u > 0 and u < mapsize and v > 0 and v < mapsize):
                beam[y, x] = dirtyBeam[0].data[0, 0, v, u]/beamfactor

    beamFT = np.fft.fftshift(fft2(beam))
    result = np.absolute(ifft2(np.fft.ifftshift(np.multiply(atomFT, beamFT))))

    plt.imsave("zbeam.png", beam)
    plt.imsave("atomft.png", np.absolute(atomFT))
    plt.imsave("beamft.png", np.absolute(beamFT))
    plt.imsave("totalft.png", np.absolute(atomFT*beamFT))

    plt.imsave("resid.png", result)
    plt.imsave("atoms.png", atoms)
    plt.imsave("atoms2.png", np.absolute(ifft2(np.fft.ifftshift(atomFT))))

    clresult = np.zeros(shape=(insize, insize), dtype=float)
    for x in range(0, insize):
        for y in range(0, insize):
            clresult[y, x] = result[int(y+inset), int(x+inset)]

    sqsum = 0.0
    dmap = dirtyMap[0].data[0, 0]
    for y in range(insize):
        for x in range(insize):
            sqsum = sqsum + (dmap[y, x]-clresult[y, x])**2
    print("RMS residual: {}".format(np.sqrt(sqsum/(insize*insize))))
