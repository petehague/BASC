#!/usr/bin/env python

from astropy.io import fits
import sys
from pymc import deterministic, Uniform, Exponential, Poisson, potential
import numpy as np

sys.path.append(os.getcwd()+"/build/lib.macosx-10.12-intel-2.7/")

import skimage as sk

def fitsraster(image, x, y):
    result = []
    count = 0
    for v in range(y):
        for u in range(x):
            result.append(image[0,0,v,u])
            count += 1
    return result

ff = fits.open(sys.argv[1])
ix = ff[0].header['NAXIS1']
iy = ff[0].header['NAXIS2']
imagedata = fitsraster(ff[0].data,ix,iy)

ff = fits.open(sys.argv[2])
bx = ff[0].header['NAXIS1']
by = ff[0].header['NAXIS2']
beamdata = fitsraster(ff[0].data,bx,by)

sk.render(imagedata,ix,iy)
#sk.render(beamdata,bx,by)

sk.makeMap(imagedata,ix,iy)
sk.makeBeam(beamdata,bx,by)

sigma = sk.mapNoise()
sigmasq = sigma*sigma

sigmasq = 1e-3

base = 0

fluxscale = sk.fluxScale()

print("Noise = {} Jy/Beam".format(sigma))
print("Base Flux = {} Jy/Beam".format(fluxscale))

#PyMC stuff

natoms = 1

xpos = []
ypos = []
flux = []

for i in range(0,natoms):
    xpos.append(Uniform('xpos{}'.format(i), lower=0, upper=511))
    ypos.append(Uniform('ypos{}'.format(i), lower=0, upper=511))
    flux.append(Exponential('flux{}'.format(i), beta=sigma))

@deterministic
def chisq(x=xpos,y=ypos,f=flux):
    atoms = []
    for i in range(0,natoms):
        atoms.append((x[i],y[i],fluxscale*np.exp(f[i])))
    return base+(1./sigmasq)*sk.deconvolve(atoms)

@potential
def logfitness(c=chisq):
    return -0.5*c

#Diagnostic main function
if __name__ == "__main__":
    output = open("like.txt","w")
    for y in range(0,512):
        for x in range(0,512):
            atoms = [(x,y,1)]
            result = sk.deconvolve(atoms)
            output.write("{} ".format(result))
        output.write("\n")
    output.close()
