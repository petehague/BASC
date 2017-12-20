#!/usr/bin/env python
import sys
import os
import numpy as np
import re

sys.path.append(os.path.dirname(os.path.realpath(__file__))+"/build/lib.macosx-10.13-intel-2.7/")

import bascmod
from astropy.io import fits
from astropy.table import Table

bascmod.init()

'''
Utilities that don't invoke bascmod
'''

def readConfig(filename):
    optionsfile = open(filename, "r")
    for line in optionsfile:
        tokens = re.split("=", line)
        key = tokens[0].strip()
        value = tokens[1].strip()
        bascmod.option(key, value)

def setOption(key, value):
    bascmod.option(key,"{}".format(value))

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

'''
The class that encapsulates all bascmod calls
'''

class view():
    def __init__(self):
        self.mx = 0
        self.my = 0
        self.crpix1 = 0
        self.crval1 = 0
        self.cdelt1 = 0
        self.crpix2 = 0
        self.crval2 = 0
        self.cdelt2 = 0
        self.cindex= bascmod.new()


    def loadFitsFile(self, filename, index):
        source = fits.open(filename)
        mx = source[0].header['NAXIS1']
        my = source[0].header['NAXIS2']
        mapdata = fitsraster(source[0].data, mx, my)
        bascmod.map(self.cindex, mapdata , mx, my, index)
        self.crpix1 = source[0].header['CRPIX1']
        self.crval1 = source[0].header['CRVAL1']
        self.cdelt1 = source[0].header['CDELT1']
        self.crpix2 = source[0].header['CRPIX2']
        self.crval2 = source[0].header['CRVAL2']
        self.cdelt2 = source[0].header['CDELT2']


    def loadMap(self,filename):
        self.loadFitsFile(filename, 0)

    def loadBeam(self,filename):
        self.loadFitsFile(filename, 1)

    def loadPBCor(self,filename):
        self.loadFitsFile(filename, 2)

    def blankPBCor(self,mx, my):
        bascmod.map(self.cindex, np.ones(shape=(mx * my)).tolist(), mx, my, 2)

    def map(self,index):
        return mapraster(bascmod.getmap(self.cindex, index))

    def run(self):
        bascmod.run(self.cindex)

    def showall(self):
        bascmod.show(self.cindex,0)
        bascmod.show(self.cindex,1)
        bascmod.show(self.cindex,2)

    def getChain(self):
        x = bascmod.chain(self.cindex,0)
        y = bascmod.chain(self.cindex,1)
        f = bascmod.chain(self.cindex,2)
        k = bascmod.chain(self.cindex,3)
        result = Table([x, y, f, k], names=('x', 'y', 'F', 'k'))
        return result

    def getSlice(self,k):
        result = self.getChain()
        models = result.group_by('k')
        mask = []
        for index in range(1, len(models.groups.keys)):
            n = models.groups.indices[index] - models.groups.indices[index - 1]
            if n == k:
                mask += np.arange(models.groups.indices[index - 1], models.groups.indices[index]).tolist()
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
