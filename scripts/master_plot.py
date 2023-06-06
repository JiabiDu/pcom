#!/usr/bin/env python3
from pylib import *
from pcom import *
from fgom_local import *
close("all")
p=zdata()
p.run='run17e2'
#p.run='run18a2'
p.brun='run17e'
#p.brun=None
p.reftime=datenum(2018,1,1)

p.salt_plot    = 1 
p.current_plot = 0
p.elev_plot    = 0
p.temp_plot    = 0
p.salt_obs='../Observations/twdb_wq/npz/twdb_salinity.npz'
p.elev_obs='../Observations/noaa_wl_temp/npz/hourly_height_2000_2022.npz'
p.elev_obs_info='../Observations/noaa_wl_temp/npz/stainfo.npz'
p.temp_obs='../Observations/noaa_wl_temp/npz/water_temperature_2000_2022.npz'
p.gd=read_schism_hgrid(f'hgrid.ll')
os.makedirs(p.run+'/figs',exist_ok=True)
ioff()
#%%
if p.salt_plot:
    M=loadz(p.run+'/salt.npz')
    O=loadz(p.salt_obs)
    B=None
    if p.brun!=None: B=loadz(p.brun+'/salt.npz')
    stations=[i for i in M.bp.station if not i.startswith('87')]; stations=unique(stations)
    lons=[mean(M.bp.x[M.bp.station==i]) for i in stations]
    stations=array(stations)[argsort(lons)]
    plot_salinity(stations,M,B,O,row=6,col=4,p=p)
    bay_station={'Galveston':[('BOLI','MIDG','FISH','TRIN'),(3,3,1.5,1.5)],
                 'Laguna_Corpus':[['LMA2','BIRD','SPCG','INPT'],(1,1,1,1)],
                 'Saint':[('CHKN','GEA-S','GEA-B','DELT','MOSQ','SOPA'),(1,1,1,1,1,1)],
                 'Matagorda':[('TPAL','6985','6990','6980','EMB'),(1,1,1,1,1)]}         
    for fmt in [0,1]: 
        for bay in ['Galveston']:
            plot_salinity_one_layer(M,B,O,fmt=fmt,stations=bay_station[bay][0],deps=bay_station[bay][1],bay=bay,p=p)

#%%
if p.elev_plot:
    run,brun,reftime=p.run,p.brun,p.reftime
    M=loadz(p.run+'/elevation.npz')
    B=None
    if p.brun!=None: B=loadz(p.brun+'/elevation.npz')
    O=loadz(p.elev_obs)
    # get station information
    E=loadz(p.elev_obs_info)
    p.lon={i:j for i,j in zip(E.station,E.lon)}
    p.lat={i:j for i,j in zip(E.station,E.lat)}
    p.station_name={i:j for i,j in zip(E.station,E.station_name)}
    stations=[i for i in M.bp.station if i.startswith('87')]
    plot_elevation(stations,M,B,O,row=6,col=4,margin=[0.05,0.02,0.06,0.05],dxy=[0.03,0.00],fmt=0,xlims=[datenum(2018,6,1),datenum(2018,7,1)],p=p) #full signal
    plot_elevation(stations,M,B,O,row=6,col=4,margin=[0.05,0.02,0.06,0.05],dxy=[0.03,0.00],fmt=2,xlims=[datenum(2018,6,1),datenum(2018,7,1)],ylims=[-0.5,0.5],p=p) #tidal signal
    plot_elevation(stations,M,B,O,row=6,col=4,margin=[0.05,0.02,0.06,0.05],dxy=[0.03,0.00],fmt=1,p=p) #subtidal

#%% 
if p.current_plot:
    O=loadz('../Observations/noaa_current/noaa_current_gom_2018.npz')
    M=loadz(p.run+'/current.npz')
    M.speed=sqrt(M.horizontalVelX**2+M.horizontalVelY**2)
    B=None
    if p.brun!=None: 
        B=loadz(p.brun+'/current.npz')
        B.speed=sqrt(B.horizontalVelX**2+B.horizontalVelY**2)
    gd=read_schism_hgrid(f'hgrid.ll')
    plot_obs_station(M,gd=gd,var='current',p=p)
    plot_current_direction_hist(M,O,p=p)
    plot_current_speed_hist(M,O,p=p)
    plot_current(M,O,B,p=p)
    plot_current(M,O,B,p=p,zoom=30) #zoom in view the first obs-available 30 days
#%% 
if p.temp_plot:
    O=loadz(p.temp_obs)
    M=loadz(p.run+'/temp.npz')
    B=None
    if p.brun!=None and fexist(p.brun+'/temp.npz'): B=loadz(p.brun+'/temp.npz')
    plot_temp(M.bp.station,M,B,O,row=6,col=4,margin=[0.05,0.02,0.06,0.05],dxy=[0.03,0.00],p=p)
