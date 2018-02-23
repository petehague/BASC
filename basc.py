#!/usr/bin/env python
import sys
import os
import numpy as np
import re

#sys.path.append(os.path.dirname(os.path.realpath(__file__))+"/build/lib.macosx-10.13-intel-2.7/")

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

def mapprof(maparr, limits=[0,0]):
    if limits==[0,0]:
        maxval = np.max(maparr)
        minval = np.min(maparr)
    else:
        minval = limits[0]
        maxval = limits[1]
    fbins = np.zeros(100)
    for f in np.ndarray.flatten(maparr):
        index = int(np.floor(100*((f-minval)/(maxval-minval))))
        if index>=0 and index<100:
            fbins[index]+=1
    return np.linspace(minval, maxval, num=100),fbins

def arrshift(maparr):
    mx,my = maparr.shape
    newmap = np.zeros(maparr.shape)
    for y in range(my):
        for x in range(mx):
            newmap[x,y] = maparr[int((x+mx/2)%mx),int((y+my/2)%my)]
    return newmap

def logimage(maparr):
    minval = np.min(maparr)-1e-6
    maparr -= minval
    return np.log10(maparr)

def cutin(maparr):
    nx,ny = maparr.shape
    nx = int(nx*2)
    ny = int(ny*2)
    offx = int(nx/4)
    offy = int(ny/4)
    result = np.zeros(shape=(nx,ny))
    for y in range(int(ny/2)):
        for x in range(int(nx/2)):
            result[x+offx,y+offy] = maparr[x,y]
    return result

def cutout(maparr):
    nx,ny = maparr.shape
    nx = int(nx/2)
    ny = int(ny/2)
    offx = int(nx/2)
    offy = int(ny/2)
    result = np.zeros(shape=(nx,ny))
    for y in range(ny):
        for x in range(nx):
            result[x,y] = maparr[x+offx,y+offy]
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
        self.resid = []


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

    def setNoise(self, noise):
        bascmod.noise(self.cindex, noise)

    def map(self,index):
        return mapraster(bascmod.getmap(self.cindex, index))

    def run(self):
        return bascmod.run(self.cindex)

    def showall(self):
        bascmod.show(self.cindex,0)
        bascmod.show(self.cindex,1)
        bascmod.show(self.cindex,2)

    def getEvidence(self):
        return bascmod.evidence(self.cindex)

    def getChain(self):
        x = bascmod.chain(self.cindex,0)
        y = bascmod.chain(self.cindex,1)
        f = bascmod.chain(self.cindex,2)
        k = bascmod.chain(self.cindex,3)
        L = bascmod.chain(self.cindex,4)
        result = Table([x, y, f, k, L], names=('x', 'y', 'F', 'k', 'L'))
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

    def getRMS(self):
        dmap = mapraster(bascmod.getmap(self.cindex, 0)) #Try implementing some caching
        dbeam = mapraster(bascmod.getmap(self.cindex, 1))
        result = self.getChain()

        offsetx, offsety = dbeam.shape
        xygrid = np.zeros(shape=dmap.shape)
        oldk = -1
        nmodels = 0
        for line in result:
            x = int(np.floor(line['x']))
            y = int(np.floor(line['y']))
            xygrid[x,y] += line['F']
            if line['k']!=oldk:
                nmodels += 1
            oldk = line['k']
        xygrid /= nmodels

        ftbeam = np.fft.fft2(arrshift(dbeam))
        ftpoints = np.fft.fft2(arrshift(cutin(xygrid)))
        convmap = ftbeam*ftpoints

        propmap = arrshift(np.fft.ifft2(convmap).real)
        propmap = np.rot90(np.fliplr(propmap)) 

        self.resid = dmap-cutout(propmap)
        rms = np.std(self.resid)

        return rms

    def getResid(self):
        return self.resid

if __name__ == "__main__":
    newView = view()
    print("Load Image")
    newView.loadMap("ex_image.fits")
    print("Load PSF")
    newView.loadBeam("ex_psf.fits")
    print("Load Primary Beam")
    newView.loadPBCor("ex_flux.fits")
    newView.showall()
    print("Run BASC")
    newView.run()
