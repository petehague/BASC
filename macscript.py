#!/usr/bin/env python3

import os
import sys
import re

'''

    Big old hack to get around the problems caused by Anaconda and its implementation of gcc
    (I shouldn't be the one to have to sort this out should I?)

'''

oldpath = os.environ['PATH']

pathnames = re.split(":",oldpath)

newpath = ""
for pn in pathnames:
  if 'naconda' not in pn:
    newpath += pn + ":"

newpath = newpath[:-1]

os.environ['PATH']=newpath
os.system("make")
os.environ['PATH']=oldpath
