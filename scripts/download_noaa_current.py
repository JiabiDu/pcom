#!/usr/bin/env python3
from pylib import *
close("all")

stations=['g08010', #available since 2015-11, bin 1-40
          'g06010', #available since 2009-2, bin 1-20
         ]
for station in stations:
  for year in arange(2009,2022):
    mday=[31,28,31,30,31,30,31,31,30,31,30,31]
    if year%4==0: mday[1]=29
    for mon in arange(1,13):
      begin_date='{:4d}{:02d}{:02d}'.format(year,mon,1)
      end_date='{:4d}{:02d}{:02d}'.format(year,mon,mday[mon-1])
      for bin in arange(1,41):       
        url='https://api.tidesandcurrents.noaa.gov/api/prod/datagetter?begin_date={}&end_date={}&station={}&product=currents&time_zone=gmt&units=metric&format=csv&bin={}'.format(begin_date,end_date,station,bin)
        fname='data/{}_{}_{}_bin{:02d}.csv'.format(station,begin_date,end_date,bin)
        if os.path.isfile(fname): continue
        try: #when there is no data at the given station, error may occurs. So use "try"
          urlsave(url,fname)
          print('.. save into {}'.format(fname))
        except:
          print('... bad request for {} (likely data do not exist)'.format(product_name))
          pass

