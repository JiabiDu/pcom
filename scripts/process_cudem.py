#!/usr/bin/env python3
from pylib import *
close("all")

fnames=[i for i in os.listdir('.') if i.endswith('.tif')]
for m,fname in enumerate(fnames):
    print(fname,'ncei19_CB_{}.npz'.format(m+1))
    convert_dem_format(fname,sname='ncei19_CB_{}.npz'.format(m+1),fmt=1)
