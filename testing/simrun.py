import os
import sys
import re

acfg = "2.6"
mapsize = 2048
cellsize = '0.01arcsec'
nmaps = 1

for i in range(nmaps):
    print("Processing simmap{}".format(i))

    os.system("rm -rf simmap{}".format(i))

    simalma(project="new{}".format(i),
            skymodel="newmap_{}.fits".format(i),
            antennalist="alma.cycle{}.cfg".format(acfg),
            totaltime="698s",
            imsize=[mapsize,mapsize],
            cell=cellsize,
            overwrite=True,
            niter=0,
            dryrun=False)

    simalma(project="newc{}".format(i),
            skymodel="newmap_{}.fits".format(i),
            antennalist="alma.cycle{}.cfg".format(acfg),
            totaltime="698s",
            imsize=[mapsize,mapsize],
            cell=cellsize,
            overwrite=True,
            niter=1000,
            dryrun=False)

    exportfits(imagename="newc{0}/newc{0}.alma.cycle{1}.noisy.image".format(i,acfg), fitsimage="new{0}/new{0}.alma.cycle{1}.clean.image.fits".format(i,acfg))
    exportfits(imagename="new{0}/new{0}.alma.cycle{1}.noisy.image".format(i,acfg), fitsimage="new{0}/new{0}.alma.cycle{1}.noisy.image.fits".format(i,acfg))
    exportfits(imagename="new{0}/new{0}.alma.cycle{1}.noisy.psf".format(i,acfg), fitsimage="new{0}/new{0}.alma.cycle{1}.noisy.psf.fits".format(i,acfg))
    exportfits(imagename="new{0}/new{0}.alma.cycle{1}.noisy.flux".format(i,acfg), fitsimage="new{0}/new{0}.alma.cycle{1}.noisy.flux.fits".format(i,acfg))
    exportuvfits(vis='new{0}/new{0}.alma.cycle{1}.noisy.ms'.format(i, acfg),
                 fitsfile='new{0}/new{0}.uv.fits'.format(i),
                 multisource=False)
