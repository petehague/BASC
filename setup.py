#!/usr/bin/env python

from distutils.core import setup, Extension
from shutil import copy
import glob

skmodule = Extension('skimage',
                     sources=['source/pyskimage.cpp', 'source/skimage.cpp'],
                     extra_compile_args=['-std=c++11'])
setup(ext_modules=[skmodule])

bascmodule = Extension('bascmod',
                       sources=['source/pybasc.cpp', 'source/skimage.cpp',
                                'bayesys/bayesys3.c', 'bayesys/random.c',
                                'bayesys/hilbert.c', 'bayesys/app.c',
                                'source/options.cpp'])

setup(ext_modules=[bascmodule])

filelist = glob.glob("build/lib*")
# copy(filelist[0]+"/bascmod.so", ".")
# copy(filelist[0]+"/skimage.so", ".")
