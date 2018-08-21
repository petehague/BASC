import os

acfg = "2.6"
imsize = 2048
cellsize = '0.01arcsec'
nmaps = 5

for i in range(nmaps):
    print("Processing snmap{}".format(i))

    os.system("rm -rf sn{}".format(i))

    simalma(project="sn{}".format(i),
            skymodel="snmap_{}.fits".format(i),
            antennalist="alma.cycle2.6.cfg",
            totaltime="698s",
            imsize=[imsize,imsize],
            cell=cellsize,
            overwrite=True,
            niter=0,
            dryrun=False)

    exportfits(imagename="sn{0}/sn{0}.alma.cycle{1}.noisy.image".format(i,acfg), fitsimage="sn{0}/sn{0}.alma.cycle{1}.noisy.image.fits".format(i,acfg))
    exportfits(imagename="sn{0}/sn{0}.alma.cycle{1}.noisy.psf".format(i,acfg), fitsimage="sn{0}/sn{0}.alma.cycle{1}.noisy.psf.fits".format(i,acfg))
    exportfits(imagename="sn{0}/sn{0}.alma.cycle{1}.noisy.flux".format(i,acfg), fitsimage="sn{0}/sn{0}.alma.cycle{1}.noisy.flux.fits".format(i,acfg))

