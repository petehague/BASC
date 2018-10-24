#!/usr/bin/env python
import sys
import os
import numpy as np
import re
import bascmod
from astropy.io import fits
from astropy.table import Table
import clustering

bascmod.init()

bascopts = {}

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
        bascopts[key] = value


def setOption(key, value):
    bascmod.option(key,"{}".format(value))
    bascopts[key] = value;

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
        self.propmap = []
        self.fluxmap = []
        self.dmap = []
        self.dbeam = []
        self.pbcor = []
        self.mapname = ""


    def loadFitsFile(self, filename, index):
        source = fits.open(filename)
        self.mx = source[0].header['NAXIS1']
        self.my = source[0].header['NAXIS2']
        if index==0:
            self.dmap = source[0].data[0,0]
        if index==1:
            self.dbeam = source[0].data[0,0]
        if index==2:
            self.pbcor = source[0].data[0,0]
        mapdata = fitsraster(source[0].data, self.mx, self.my)
        bascmod.map(self.cindex, mapdata , self.mx, self.my, index)
        self.crpix1 = source[0].header['CRPIX1']
        self.crval1 = source[0].header['CRVAL1']
        self.cdelt1 = source[0].header['CDELT1']
        self.crpix2 = source[0].header['CRPIX2']
        self.crval2 = source[0].header['CRVAL2']
        self.cdelt2 = source[0].header['CDELT2']


    def loadMap(self,filename):
        self.loadFitsFile(filename, 0)
        self.mapname = filename

    def loadBeam(self,filename):
        self.loadFitsFile(filename, 1)

    def loadPBCor(self,filename):
        self.loadFitsFile(filename, 2)

    def blankPBCor(self,mx, my):
        bascmod.map(self.cindex, np.ones(shape=(mx * my)).tolist(), mx, my, 2)

    def setNoise(self, noise):
        bascmod.noise(self.cindex, noise)

    def setFlux(self, flux):
        bascmod.flux(self.cindex, flux)

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
        if 'objmode' in bascopts:
            if bascopts['objmode']==1:
                pa = bascmod.chain(self.cindex,5)
                major = bascmod.chain(self.cindex,6)
                minor = bascmod.chain(self.cindex,7)
                result = Table([x,y,f,k,L,pa,major,minor],
                               names=('x','y','F','k','L','pa','maj','min'))
                return result
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
        result = self.getChain()

        xygrid = np.zeros(shape=self.dbeam.shape)
        oldk = -1
        nmodels = 0
        offx, offy = self.dbeam.shape
        ncells = offx*offy
        offx = int(offx/4)-1
        offy = int(offy/4)-1
        for line in result:
            x = int(line['x'])+offx
            y = int(line['y'])+offy
            xygrid[x,y] += line['F']
            if line['k']!=oldk:
                nmodels += 1
            oldk = line['k']
        xygrid /= nmodels
        xygrid /= ncells
        xygrid = np.rot90(xygrid)
        xygrid = np.fliplr(xygrid)
        self.fluxmap = cutout(xygrid)

        ftbeam = np.fft.fft2(arrshift(self.dbeam))
        ftpoints = np.fft.fft2(xygrid)
        convmap = ftbeam*ftpoints
        self.propmap = cutout(np.fft.fft2(convmap).real)

        self.resid = cutout(self.dmap)-self.propmap
        rms = np.std(self.resid)

        return rms

    def getResid(self):
        return self.resid

    def saveResid(self, filename):
        source = fits.open(self.mapname)
        newimage = np.ndarray(shape=(1,1, self.mx, self.my),dtype=float,buffer=cutin(self.resid))
        if os.path.exists(filename):
            os.remove(filename)
        fits.writeto(filename,newimage,source[0].header)       

    def saveResult(self, filename):
        source = fits.open(self.mapname)
        newimage = np.ndarray(shape=(1,1, self.mx, self.my),dtype=float,buffer=cutin(self.fluxmap))
        if os.path.exists(filename):
            os.remove(filename)
        fits.writeto(filename,newimage,source[0].header)       

    def saveProp(self, filename):
        source = fits.open(self.mapname)
        newimage = np.ndarray(shape=(1,1, self.mx, self.my),dtype=float,buffer=cutin(self.propmap))
        if os.path.exists(filename):
            os.remove(filename)
        fits.writeto(filename,newimage,source[0].header)       

    def clusters(self, min_samples=10, eps=2):
        result = self.getChain()
        lastk = -1
        maxk = 0
        curk = 0
        for line in result:
            if line['k']==lastk:
                curk += 1
            else:
                if curk>maxk:
                    maxk = curk
                curk = 0
            lask = line['k']
        xydata = np.zeros(shape=(len(result),2))
        xydata[:,0]=result['x'].data
        xydata[:,1]=result['y'].data
        fluxdata = result['F'].data
        atom, xy, flux, noise, labels, centers, widths = clustering.find_center(xydata, fluxdata, maxk, min_samples, eps)
        F = []
        for fluxlist in flux:
            F.append(np.mean(fluxlist))
        clout = Table([centers[:,0],centers[:,1],widths[:,0],widths[:,1], F],names=("x", "y", "dx", "dy", "F"))
        return clout,len(noise)
        

if __name__ == "__main__":
    if len(sys.argv)<4:
        print("basc.py <dirty image> <dirty beam> <primary beam flux>")
    else:
       newView = view()
       print("Load Image "+sys.argv[1])
       newView.loadMap(sys.argv[1])
       print("Load PSF "+sys.argv[2])
       newView.loadBeam(sys.argv[2])
       print("Load Primary Beam "+sys.argv[3])
       newView.loadPBCor(sys.argv[3])
       newView.showall()
       print("Run BASC")
       newView.run()
