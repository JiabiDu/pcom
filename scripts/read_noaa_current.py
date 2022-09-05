#!/usr/bin/env python3
from pylib import *
close("all")

stations=['g06010','g08010']
year=2018
for station in stations:
    fnames=[i for i in os.listdir('data/') if i.startswith(station+'_{}'.format(year))]
    time,speed,direction,binn=[],[],[],[]
    for fname in fnames:
        data=array([i.strip().split(',') for i in open('data/'+fname,'r').readlines()])
        data=data[1:] #remove the head line
        if len(data)<=10: continue #no data, can use len(data)<=1
        print(fname)
        time.extend(quickdatenum(data[:,0]))
        speed.extend(data[:,1].astype('float'))
        direction.extend(data[:,2].astype('float'))
        binn.extend(data[:,3].astype('int'))
    S=zdata()
    S.time,S.speed,S.direction,S.bin=array(time),array(speed),array(direction),array(binn)
    savez(station+'_{}.npz'.format(year),S)
        