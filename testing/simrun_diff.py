import os
import sys
from astropy.io import fits
import re

#sys.path.append("/home/prh44/stoa/ALMA/")
#import almautil

acfg = "2.6"
mapsize = 2048
cellsize = '0.01arcsec'
nmaps = 100

codelist = [1,2,3,4]


for code in codelist:
	for i in range(nmaps):
	    source = fits.open("simmap_{}.fits".format(i))
	    dobs = source[0].header['DATE-OBS'][0:10]
	    dobs = "/".join(re.split("-",dobs)) 
	    source.close()
	    print("Processing diffmap{}_{}".format(code,i))

	    os.system("rm -rf diffmap{}_{}".format(code,i))

	    simalma(project="diff{}_{}".format(code,i),
		    skymodel="diffmap{}_{}.fits".format(code,i),
		    antennalist="alma.cycle{}.cfg".format(acfg),
		    totaltime="698s",
		    imsize=[mapsize,mapsize],
		    cell=cellsize,
		    overwrite=True,
		    niter=0,
		    dryrun=False)

	    exportfits(imagename="diff{2}_{0}/diff{2}_{0}.alma.cycle{1}.noisy.image".format(i,acfg,code), fitsimage="diff{2}_{0}/diff{2}_{0}.alma.cycle{1}.noisy.image.fits".format(i,acfg,code))
	    exportfits(imagename="diff{2}_{0}/diff{2}_{0}.alma.cycle{1}.noisy.psf".format(i,acfg,code), fitsimage="diff{2}_{0}/diff{2}_{0}.alma.cycle{1}.noisy.psf.fits".format(i,acfg,code))
	    exportfits(imagename="diff{2}_{0}/diff{2}_{0}.alma.cycle{1}.noisy.flux".format(i,acfg,code), fitsimage="diff{2}_{0}/diff{2}_{0}.alma.cycle{1}.noisy.flux.fits".format(i,acfg,code))
