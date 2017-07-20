#!/usr/bin/env python

from distutils.core import setup, Extension

skmodule = Extension('skimage',
                     sources=['source/pyskimage.cpp', 'source/skimage.cpp'],
                     extra_compile_args=['-std=c++11'])
setup(ext_modules=[skmodule])
