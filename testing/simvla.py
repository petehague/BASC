## Generate the visibilities for the aliasing experiments.

import os
import sys
import glob
from scipy import constants
import numpy as np

a = np.loadtxt('pointlog.txt', delimiter = ' ').T
index = a[0]
r = a[1]//10 * 3 
theta = a[2]
loc_x_pixel = np.round(r*np.cos(theta))
loc_y_pixel = np.round(r*np.sin(theta))


pixel_size = 0.03
miu = 4.0e+09

for i in range(len(r)):##2 sources
	print ('processing file ', i)
	project_name = 'vla_sim_' + str(i)
	os.system('rm -rf Point_source.cl')
	os.system('rm -rf ' + project_name)
	os.system('rm -rf ' + project_name + '.dirty.image.fits')
	os.system('rm -rf ' + project_name + '.dirty.flux.fits')
	os.system('rm -rf ' + project_name + '.dirty.psf.fits')
	os.system('rm -rf ' + project_name + '.clean.image.fits')
	os.system('rm -rf ' + project_name + '.clean.flux.fits')
	loc_x = loc_x_pixel[i] * pixel_size /15 #0.002 arcsec/pixel, 15 pixel
	loc_y = loc_y_pixel[i] * pixel_size 
	if loc_x >= 0 and loc_y >= 0:
		cl.addcomponent(dir='J2000 0h0m' + str(loc_x) + 's 0d0m' + str(loc_y) + 's', flux=0.05, fluxunit='Jy', freq='4.0 GHz', spectrumtype='Constant', shape='point')
	elif loc_x >= 0 and loc_y <= 0:
		cl.addcomponent(dir='J2000 0h0m' + str(loc_x) + 's -0d0m' + str(-loc_y) + 's', flux=0.05, fluxunit='Jy', freq='4.0 GHz', spectrumtype='Constant', shape='point')
	elif loc_x <= 0 and loc_y >= 0:
		cl.addcomponent(dir='J2000 -0h0m' + str(-loc_x) + 's 0d0m' + str(loc_y) + 's', flux=0.05, fluxunit='Jy', freq='4.0 GHz', spectrumtype='Constant', shape='point')
	elif loc_x <= 0 and loc_y <= 0:
		cl.addcomponent(dir='J2000 -0h0m' + str(-loc_x) + 's -0d0m' + str(-loc_y) + 's', flux=0.05, fluxunit='Jy', freq='4.0 GHz', spectrumtype='Constant', shape='point')
	cl.addcomponent(dir='J2000 0h0m0s 0d0m0s', flux=0.05, fluxunit='Jy', freq='4.0 GHz', spectrumtype='Constant', shape='point')
	cl.rename('Point_source.cl')
	cl.done()

	simobserve(project = project_name,
		complist = 'Point_source.cl',
		direction = "J2000 0h0m0s 0d0m0s",
		obsmode = 'int',
		antennalist = 'vla.a.cfg',
		totaltime = '1h',
		thermalnoise = 'tsys-atm')
	os.chdir(project_name + '/')
	## make .csv vis file ####
	sm.openfromms(project_name + '.vla.a.ms') 
	sm.setnoise(mode='simplenoise', simplenoise='2Jy')   
	sm.corrupt()
	sm.done()
	sm.close()  
	wavelength = constants.c/miu
	np.set_printoptions(threshold='nan') #This prevents numpy truncating the print out of large arrays.
	tb.open(project_name + '.vla.a.ms', nomodify=F) #Opens your Measurement Set
	uvw = tb.getcol('UVW') #Puts UVW data information into array uvw.
	uvw = uvw.transpose() #transposes uvw array into uvw_col such that there are 3 columns representing u v and w values.
	uvw = uvw/wavelength
	mydata = tb.getcol('DATA') #Puts visibility data information into array mydata.
	data = mydata[0].transpose()
	wt = np.ones([len(data),1])
	vis = np.concatenate((uvw,data.real,data.imag,wt),axis=1)
	np.savetxt(project_name + '.csv',vis, delimiter = ',') # saves uvw_col to a text file named 'namefile'.
	tb.close()#Closes your Measurement Set

	clean(vis = project_name + '.vla.a.ms',
	 	imagename = project_name +'.dirty',
		imsize = 500,
		niter = 0,
		cell = pixel_size,
		usescratch = F,
		interactive = False,
		weighting = 'natural')

	clean(vis = project_name + '.vla.a.ms',
	 	imagename = project_name +'.clean',
		imsize = 500,
		niter = 1000,
		cell = pixel_size,
		usescratch = F,
		interactive = False,
		weighting = 'natural')


	exportfits(imagename = project_name +'.clean.model', fitsimage = '../' + project_name +'.clean.model.fits')
	exportfits(imagename = project_name +'.dirty.image', fitsimage = '../' + project_name +'.dirty.image.fits')
	exportfits(imagename = project_name +'.dirty.flux', fitsimage = '../' + project_name +'.dirty.flux.fits')
	exportfits(imagename = project_name +'.dirty.psf', fitsimage = '../' + project_name +'.dirty.psf.fits')
	exportfits(imagename = project_name +'.clean.image', fitsimage = '../' + project_name +'.clean.image.fits')
	exportfits(imagename = project_name +'.clean.flux', fitsimage = '../' + project_name +'.clean.flux.fits')
	os.chdir('../')
	
