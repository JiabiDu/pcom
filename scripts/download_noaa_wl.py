#!/usr/bin/env python3
from pylib import *
close("all")
from functions import *
#--------------------------------------------
read_stainfo=False
xlims=[-98,-87]
ylims=[25,31]
product_names=['hourly_height','water_temperature','wind']
years=arange(2000,2023)

#-------------------------------------------
#%% get station info
if read_stainfo or not os.path.exists('npz/stainfo.npz'):
    f=open('stainfo_noaa.txt','r')
    lon,lat,station,station_name=[],[],[],[]
    ist=0
    for i in f.readlines():
      i=i.split('"')
      if len(i)<8: continue
      station.append(i[1])
      station_name.append(i[3])
      lat.append(i[5])
      lon.append(i[7])
      ist+=1
      print(ist)
    S=zdata()
    S.lon,S.lat,S.station,S.station_name=array(lon).astype('float'), array(lat).astype('float'),array(station),array(station_name)
    savez('npz/stainfo.npz',S)
S=loadz('npz/stainfo.npz')
lon,lat,station=S.lon,S.lat,S.station
fp=(lon>xlims[0])*(lon<xlims[1])*(lat>ylims[0])*(lat<ylims[1])
stations=station[fp]
print('# of noaa stations to download: {}'.format(sum(fp))) 

#------------------------------------------    
#%% get and process noaa data
get_noaa_tide_current(stations=stations,years=years,varnames=product_names,load_again=False)
process_noaa_tide_current(stations=stations,years=years,varnames=product_names,read_again=False)

