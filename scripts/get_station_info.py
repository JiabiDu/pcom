#!/usr/bin/env python3
from pylib import *
import json

#%%
#get noaa tide and current stations, lon, lat, has_tide, has_current, has_meteo, active
active={}
historic={}
stations=[]
lon={}
lat={}
name={}
# get historical tide stations
lines=open('html/noaa_tide_historic.html','r').readlines()
substr='.html?id='

tstation=[]
for line in lines:
    res = [i for i in range(len(line)) if line.startswith(substr, i)]
    for i in res:
        tstation.append(line[i+18:i+25])
tstation=unique(tstation)

for i in tstation:
    historic[i]=True
    
# get active tide stations
lines=open('html/noaa_tide_active.html','r').readlines()

tstation=[]
for line in lines:
    res = [i for i in range(len(line)) if line.startswith(substr, i)]
    for i in res:
        tstation.append(line[i+18:i+25])
tstation=unique(tstation)

for i in tstation:
    active[i]=True

stations = array([i for i in historic])
for station in stations:
    if not station in active.keys():
        active[station]=False
for m,station in enumerate(stations):
    url='https://api.tidesandcurrents.noaa.gov/mdapi/prod/webapi/stations/{}.json?expand=details,products&units=metric'.format(station)
    name[station]=station
    lon[station]=nan
    lat[station]=nan
    try: 
        urlsave(url,'tmp.json')
    except:
        print('no data found for station ',station)
        continue
    with open('tmp.json') as f: jdata = json.load(f)
    tlat,tlon,tname=jdata['stations'][0]['lat'],jdata['stations'][0]['lng'],jdata['stations'][0]['name']
    print('{}/{} {}:{}, {}N, {}W'.format(m+1,len(stations),station,tname,tlat,-tlon))
    name[station]=tname
    lon[station]=tlon
    lat[station]=tlat
lons=array([lon[i] for i in stations])
lats=array([lat[i] for i in stations])
S=zdata()
S.lon,S.lat,S.name,S.lons,S.lats,S.stations,S.active,S.historic=lon,lat,name,lons,lats,stations,active,historic
savez('noaa_tide_station_info.npz',S)
os.remove('tmp.json')

#%%
#get noaa tide and current stations, lon, lat, has_tide, has_current, has_meteo, active
active={}
historic={}
stations=[]
lon={}
lat={}
name={}
# get historical tide stations
lines=open('html/noaa_meteo_historic.html','r').readlines()
substr='.html?id='

tstation=[]
for line in lines:
    res = [i for i in range(len(line)) if line.startswith(substr, i)]
    for i in res:
        tstation.append(line[i+18:i+25])
tstation=unique(tstation)

for i in tstation:
    historic[i]=True
    active[i]=False

# get active tide stations
lines=open('html/noaa_meteo_active.html','r').readlines()

tstation=[]
for line in lines:
    res = [i for i in range(len(line)) if line.startswith(substr, i)]
    for i in res:
        tstation.append(line[i+18:i+25])
tstation=unique(tstation)

for i in tstation:
    active[i]=True

for station in active:
    if not station in historic.keys():
        historic[station]=False

stations=array([i for i in active])

for m,station in enumerate(stations):
    url='https://api.tidesandcurrents.noaa.gov/mdapi/prod/webapi/stations/{}.json?expand=details,products&units=metric'.format(station)
    name[station]=station
    lon[station]=nan
    lat[station]=nan
    try: 
        urlsave(url,'tmp.json')
    except:
        print('no data found for station ',station)
        continue
    with open('tmp.json') as f: jdata = json.load(f)
    tlat,tlon,tname=jdata['stations'][0]['lat'],jdata['stations'][0]['lng'],jdata['stations'][0]['name']
    print('{}/{} {}:{}, {}N, {}W'.format(m+1,len(stations),station,tname,tlat,-tlon))
    name[station]=tname
    lon[station]=tlon
    lat[station]=tlat
lons=array([lon[i] for i in stations])
lats=array([lat[i] for i in stations])
S=zdata()
S.lon,S.lat,S.name,S.lons,S.lats,S.stations,S.active,S.historic=lon,lat,name,lons,lats,stations,active,historic
savez('noaa_meteo_station_info.npz',S)
os.remove('tmp.json')
#%% get current information
active={}
historic={}
stations=[]
lon={}
lat={}
name={}

data=array([array(i.strip().split(','),dtype=object) for i in open('html/coops-historiccurrentstations.csv','r').readlines()])
data=data[1:]
tstation=data[:,0]
tname=data[:,1]
tlat=data[:,2].astype('float')
tlon=data[:,3].astype('float')

for station,namei,loni,lati in zip(tstation,tname,tlon,tlat):
    lon[station]=float(loni)
    lat[station]=float(lati)
    name[station]=namei.strip()
    historic[station]=True
    active[station]=False
stations=tstation

data=array([array(i.strip().split(','),dtype=object) for i in open('html/coops-activecurrentstations.csv','r').readlines()])
data=data[1:]
tstation=data[:,0]
tname=data[:,1]
tlat=data[:,2]
tlon=data[:,3]
for station,namei,loni,lati in zip(tstation,tname,tlon,tlat):
    lon[station]=float(loni)
    lat[station]=float(lati)
    name[station]=namei.strip()
    active[station]=True
    
for station in active:
    if not station in historic.keys():
        historic[station]=False
        
stations=array([i for i in active])       
lons=array([lon[i] for i in stations]).astype('float')
lats=array([lat[i] for i in stations]).astype('float')
S=zdata()
S.lon,S.lat,S.name,S.lons,S.lats,S.stations,S.active,S.historic=lon,lat,name,lons,lats,stations,active,historic
savez('noaa_current_station_info.npz',S)
