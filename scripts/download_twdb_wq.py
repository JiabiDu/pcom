#!/usr/bin/env python3
"""
This file will automatically download the latest archive of salinity from TWDB
and save into csv format. 
See https://www.waterdatafortexas.org/api#introduction for the api
"""
from pylib import *
import requests
import json
import csv
import os
from datetime import datetime
#-----------------------------------------------------------------
reload=True #reddownalod the data
save_csv=False #if False, data will be only saved into npz
vars=[]
sname='npz/twdb_wq.npz'
vars=['seawater_salinity'] #if vars is empty, all available parameters will be downloaded
sname='npz/twdb_salinity.npz'

#-----------------------------------------------------------------
now = datetime.now()
now = str(now)[0:10] #in format yyyy-mm-dd

if not os.path.exists('data'): os.mkdir('data')
if not os.path.exists('npz'): os.mkdir('npz')

#-- get station information
payload = {"all_stations": True}
r = requests.get("https://waterdatafortexas.org/coastal/api/stations", params=payload)
station_info = r.json()

#- get station list
stations,lons,lats=[],[],[]
lon,lat={},{}
for rc in station_info:
    stations.append(rc['station_code'])
    lons.append(rc['coordinates'][0])
    lats.append(rc['coordinates'][1])
    lon[stations[-1]]=rc['coordinates'][0]
    lat[stations[-1]]=rc['coordinates'][1]

#-- download data
ist=0
stationi=[]
datai=[]
timei=[]
vari=[]
for m,station in enumerate(stations):
    #- get parameters at this station
    r = requests.get(f'https://waterdatafortexas.org/coastal/api/stations/{station}/parameters')
    vars_info=r.json()
    vars2=[i['code'] for i in vars_info]
    if len(vars)==0: vars=vars2
    
    #- download data and save inoto csv format
    for var in vars: 
        if not var in vars2: continue #in case the specified variable is not found
        payload = {"start_date": "1900-01-01", "end_date": now, "binning": "hour"}
        r = requests.get("https://waterdatafortexas.org/coastal/api/stations/{}/data/{}".format(station,var), params=payload)
        json_data=json.loads(r.text)
        if save_csv:
          fname="data/{}_{}.csv".format(station,var)
          if os.path.isfile(fname) and not reload: print(fname,'exist'); continue 
          if len(json_data)!=0:
            keys = json_data[0].keys()
            a_file = open(fname, "w", newline='')
            dict_writer = csv.DictWriter(a_file, keys)
            dict_writer.writeheader()
            dict_writer.writerows(json_data)
            a_file.close()
            print('#{} got {} for station {} from 1900-01-01 to {}'.format(m+1,var,station,now))
    
        #- read data and save into npz fomrat, need station, parameter, depth
        rn=len(json_data) 
        if rn==0: continue
        datai.extend([i['value'] for i in json_data])
        #timei.extend(quickdatenum([i['datetime_utc'] for i in json_data])) #run into error 
        stationi.extend(tile(station,rn))
        vari.extend(tile(var,rn))
        for jsoni in json_data:
            timei.append(datenum(jsoni['datetime_utc']))
        print('got {} for station #{}/{} {} rn={}'.format(var,m+1,len(stations),station,rn))
S=zdata()
S.time=array(timei)
S.data=array(datai).astype('float')
S.station=array(stationi)
S.var=array(vari)
S.lon=lon
S.lat=lat
S.lons=lons
S.lats=lats
S.stations=stations
savez(sname,S)
