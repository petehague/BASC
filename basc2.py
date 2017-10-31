#!/usr/bin/env python
import sys
import os
import numpy as np

#sys.path.append(os.path.dirname(os.path.realpath(__file__))+"/build/lib.macosx-10.12-intel-2.7/")

import bascmod
from astropy.io import fits
from astropy.table import Table

# TODO: This needs to be repaced with a faster interface
def fitsraster(image, x, y):
    result = []
    count = 0
    for v in range(y):
        for u in range(x):
            result.append(image[0,0,v,u])
            count += 1
    return result

def mapraster(rawmap):
    mapsize = int(np.sqrt(len(rawmap)))
    result = np.zeros(shape=(mapsize,mapsize))
    for y in range(mapsize):
        for x in range(mapsize):
            result[x,y] = rawmap[y*mapsize + x]
    return result

def loadFitsFile(filename, index):
    source = fits.open(filename)
    mx = source[0].header['NAXIS1']
    my = source[0].header['NAXIS2']
    bascmod.map(fitsraster(source[0].data, mx, my), mx, my, index)
    crpix1 = source[0].header['CRPIX1']
    crval1 = source[0].header['CRVAL1']
    cdelt1 = source[0].header['CDELT1']
    crpix2 = source[0].header['CRPIX2']
    crval2 = source[0].header['CRVAL2']
    cdelt2 = source[0].header['CDELT2']


def loadMap(filename):
    loadFitsFile(filename, 0)

def loadBeam(filename):
    loadFitsFile(filename, 1)

def loadPBCor(filename):
    loadFitsFile(filename, 2)

def map(index):
    return mapraster(bascmod.getmap(index))

def run():
    bascmod.run()

def showall():
    bascmod.show(0)
    bascmod.show(1)
    bascmod.show(2)

def getChain():
    x = bascmod.chain(0)
    y = bascmod.chain(1)
    f = bascmod.chain(2)
    k = bascmod.chain(3)
    result = Table([x, y, f, k], names = ('x', 'y', 'F', 'k'))
    return result

def getSlice(k):
    result = getChain()
    models = result.group_by('k')
    mask = []
    for index in range(1,len(models.groups.keys)):
        n = models.groups.indices[index] - models.groups.indices[index-1]
        if n==k:
            mask += np.arange(models.groups.indices[index-1],models.groups.indices[index]).tolist()
    # result = Table([x[mask],y[mask],f[mask],k[mask)]], names = ('x', 'y', 'F', 'k'))
    return result[mask]


if __name__ == "__main__":
    print("Load Image")
    loadMap("ex_image.fits")
    #bascmod.show(0)
    print("Load PSF")
    loadBeam("ex_psf.fits")
    #bascmod.show(1)
    print("Load Primary Beam")
    loadPBCor("ex_flux.fits")
    #bascmod.show(2)
    print("Run BASC")
    bascmod.run()
