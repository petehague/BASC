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

# Results
atoms = 0


# Perform MCMC run
def run(dmap, dbeam, dflux):
    global chain, nmodels, mapsize, mapdepth

    dirtyMap = fits.open(dmap)
    mapsize = dirtyMap[0].header['NAXIS1']
    mapdepth = dirtyMap[0].header['NAXIS3']

    output = open("meta.txt", "w")
    output.write("{:f} {:f} {:d} {:f} {:f} {:d}\n".format(dirtyMap[0].header['CRVAL1'],
                                              dirtyMap[0].header['CDELT1'],
                                              int(dirtyMap[0].header['CRPIX1']),
                                              dirtyMap[0].header['CRVAL2'],
                                              dirtyMap[0].header['CDELT2'],
                                              int(dirtyMap[0].header['CRPIX2'])))
    output.close()

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
                    mapdata = image[0].data[0, f, y, x]
                    if np.isfinite(mapdata):
                        output.write(" {}".format(mapdata))
                    else:
                        output.write(" 0")
            output.write("\n")
        output.close()
        filelist.append(outfile)
        image.close()

    os.system(os.path.dirname(os.path.abspath(__file__)) +
              "/mcmc {} {} {} meta.txt {} {}".format(filelist[0],
                                            filelist[1],
                                            filelist[2],
                                            mapsize,
                                            mapdepth))
    chain = Table.read("chain.txt", format="ascii")
    infofile = open("info.txt", "r")
    nmodels = infofile.readline()
    infofile.close()


# Retrieve fourier transform of beam
def beamft(filename):
    dirtyBeam = fits.open(filename)

    beam = np.zeros(shape=(mapsize, mapsize), dtype=float)
    offset = mapsize/2
    beamfactor = float(mapsize*mapsize)
    for x in range(0, mapsize):
        for y in range(0, mapsize):
            u = int(y-offset)
            v = int(x-offset)
            if (u > 0 and u < mapsize and v > 0 and v < mapsize):
                beam[y, x] = dirtyBeam[0].data[0, 0, v, u]/beamfactor

    return np.fft.fftshift(fft2(beam))


# Get a map of flux weighted atom positions
def getAtoms():
    global atoms

    atoms = np.zeros(shape=(mapsize/2, mapsize/2), dtype=float)

    afactor = 1.0/float(nmodels)

    for line in chain:
        x = int(round(line[0]))
        y = int(round(line[1]))
        f = line[2]
        atoms[y, x] += f

    atoms = np.multiply(atoms, afactor)
    return atoms


# Version of the above, but pads map to make it the same size as dirty map
def getAtomsPad():
    global atoms

    atoms = np.zeros(shape=(mapsize, mapsize), dtype=float)

    afactor = 1.0/float(nmodels)
    inset = mapsize/4

    for line in chain:
        x = int(round(line[0]) + inset)
        y = int(round(line[1]) + inset)
        f = line[2]
        atoms[y, x] += f

    atoms = np.multiply(atoms, afactor)
    return atoms


# Gets the fourier transform of the atoms, padded to dirty map size
def atomft():
    return np.fft.fftshift(fft2(getAtomsPad()))


# Computes the RMS residual of the convolved atom map/beam with dirty map
def resid(filename, beamfile):
    insize = int(mapsize/2)
    inset = int(mapsize/4)

    result = ifft2(np.fft.ifftshift(np.multiply(atomft(), beamft(beamfile))))
    result = np.absolute(result)
    clresult = np.zeros(shape=(insize, insize), dtype=float)
    for x in range(0, insize):
        for y in range(0, insize):
            clresult[y, x] = result[int(y+inset), int(x+inset)]

    sqsum = 0.0
    dirtyMap = fits.open(filename)
    dmap = dirtyMap[0].data[0, 0]
    for y in range(insize):
        for x in range(insize):
            sqsum = sqsum + (dmap[y, x]-clresult[y, x])**2
    return np.sqrt(sqsum/float(insize/insize))


if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("basc.py <dirty map> <dirty psf> <dirty primary beam>")
        sys.exit(1)

    run(sys.argv[1], sys.argv[2], sys.argv[3])

    atomFT = atomft()
    beamFT = beamft(sys.argv[2])

    plt.imsave("atomft.png", np.absolute(atomFT))
    plt.imsave("beamft.png", np.absolute(beamFT))
    plt.imsave("totalft.png", np.absolute(atomFT*beamFT))
    plt.imsave("atoms.png", atoms)

    print("RMS residual: {}".format(resid(sys.argv[1], sys.argv[2])))
