#!/usr/bin/env python3
import os
fnames=open('filelist','r').readlines()
for fname in fnames:
    if os.path.exists(fname): continue
    cmd= f'wget --no-check-certificate https://coast.noaa.gov/htdata/raster2/elevation/NCEI_ninth_Topobathy_2014_8483/chesapeake_bay/{fname}'
    print(cmd); os.system(cmd)
