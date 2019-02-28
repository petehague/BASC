#!/usr/bin/env python

from distutils.core import setup, Extension
from shutil import copy
import glob
import os
import re

skmodule = Extension('skimage',
                     sources=['source/pyskimage.cpp', 'source/skimage.cpp'],
                     language='c++',
                     extra_compile_args=['-std=c++11'])
setup(ext_modules=[skmodule])

bascmodule = Extension('bascmod',
                       sources=['source/pybasc.cpp', 'source/skimage.cpp',
                                'bayesys/bayesys3.c', 'bayesys/random.c',
                                'bayesys/hilbert.c', 'bayesys/app.c',
                                'source/options.cpp'],
                       language='c++')

setup(ext_modules=[bascmodule])

filelist = glob.glob("build/lib*")
if "PYTHONPATH" in os.environ:
    pypath = os.environ["PYTHONPATH"]
    newpypath = "/".join(re.split("/",os.path.realpath(__file__))[:-1])
    newpypath += "/build/"
    newpypath += re.split("/",filelist[0])[-1]
    print("To set your python path: PYTHONPATH=$PYTHONPATH:"+newpypath+"\n")
  
libraries = glob.glob(filelist[0]+"/*.so")
for lib in libraries:
    newlib = re.split("/", lib)[-1]
    newlib = newlib[0:7]+".so"
    print(lib, newlib)
    copy(lib,newlib)
