#!/usr/bin/env Python
import sys
import os

#sys.path.append(os.path.dirname(os.path.realpath(__file__))+"/build/lib.macosx-10.12-intel-2.7/")

import bascmod
from astropy.io import fits

# TODO: This needs to be repaced with a faster interface
def fitsraster(image, x, y):
    result = []
    count = 0
    for v in range(y):
        for u in range(x):
            result.append(image[0,0,v,u])
            count += 1
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
