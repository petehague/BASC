#!/usr/bin/env python

from astropy.io import fits
import sys
import re

for arg in sys.argv[1:]:
    infile = arg
    outfile = re.split("\.",infile)[:-1]
    outfile.append("txt")
    outfile = '.'.join(outfile)

    print(infile,outfile)
    image = fits.open(infile)
    output = open(outfile, "w")
    nx = image[0].header['NAXIS1']
    ny = image[0].header['NAXIS2']
    for y in range(ny):
      for x in range(nx):
        output.write(" {}".format(image[0].data[0,0,y,x]))
      output.write("\n")
    output.close()
