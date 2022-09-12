#!/usr/bin/env python3

from pylib import *
import json
from scipy.stats import gaussian_kde
from scipy.interpolate import griddata
#================================================================
# AUTHOR NOTES
# This package is used for Coastal Ocean Modelling
# It has several components: data doaloading, data processing, plotting 
# contact jiabi.du@gmail.com
#===============================================================

#%% schism related
def gen_bpfile2(lons,lats,stations,fname='station.bp',cmt='station.bp',hgrid='hgrid.gr3',vgrid='vgrid.in'):
    '''
    gen bpfiles for extracting the vertical profiles at certain stations
    the script will find those stations within the model domain, and get the 
    vertical layer information at nearest grid. 
    ToDO: enable hgrid and vgrid loaded from npz
    '''
    if not os.path.exists(hgrid): sys.exit('hgrid not found')
    if not os.path.exists(vgrid): sys.exit('vgrid not found')
    gd=read_schism_hgrid(hgrid)
    vg=read_schism_vgrid(vgrid)
    # find nearest pts, ensure the distance is within some range
    stations,lons,lats=array(stations),array(lons),array(lats)
    
    #-- select only station within model domain
    ie,ip,acor=gd.compute_acor(c_[lons,lats])
    fp=ie!=-1 #points inside model modmain
    stations,lons,lats=stations[fp],lons[fp],lats[fp]
    z=vg.compute_zcor(gd.dp)
    sind=near_pts(c_[lons,lats],c_[gd.x,gd.y])
    zi2=z[sind]
    
    newstations,newlons,newlats,newdeps=[],[],[],[]
    for station,lon,lat,zi in zip(stations,lons,lats,zi2):
        for iz in unique(zi):
            newstations.append(station)
            newlons.append(lon)
            newlats.append(lat)
            newdeps.append(iz*-1)
    gen_bpfile(newlons,newlats,newstations,newdeps,fname=fname,cmt=cmt)

def gen_bpfile(lons,lats,stations,deps=0.0,fname='station.bp',cmt='station.bp'):
    '''
    To generate bpfile based on given lon,lat information
    ----------
    lons,lats,stations: location and names of stations
    deps : can be list or one scalar value

    '''
    f=open(fname,'w')
    f.write(cmt+'\n')
    f.write('{}\n'.format(len(lons)))
    ist=0
    if type(deps) in [int,float]: deps=ones(len(lons))*deps
    for station,lon,lat,dep in zip(stations,lons,lats,deps):
        ist+=1
        f.write('{} {:.10f} {:.10f} {:.4f} !{}\n'.format(ist,lon,lat,dep,station))
    f.close()

#%% mode-observaiton comparison related
def pair_data(x1,y1,x2,y2,hw=0.2/24):
    ''' 
    pair data within a certain of time window
    used to compare model and observation data
    put sparse data first (i.e., x1 and y1)
    '''
    py1,py2=[],[]
    for ix,iy in zip(x1,y1): 
        fp=abs(ix-x2)<hw
        if sum(fp)==0: continue
        if sum(fp)>1: print('multiple value will be avearged')
        py1.append(iy)
        py2.append(nanmean(y2[fp]))
    return array(py1),array(py2)

def kde_plot(x,y,ms=1):
    # Calculate the point density
    xy = np.vstack([x,y])
    z = gaussian_kde(xy)(xy)

    # Sort the points by density, so that the densest points are plotted last
    idx = z.argsort()
    x, y, z = x[idx], y[idx], z[idx]

    scatter(x, y, c=z, s=ms)

#%% Data downloading
def get_usgs_flow(stations=None,StartT='1980-1-1',EndT='2022-1-1',sname=None,reRead=False,sdir=None, reDownload=False):
    y1=num2date(datenum(StartT)).year; y2=num2date(datenum(EndT)-1/24).year
    if sdir is None: sdir=f'usgs_{y1}_{y2}'
    if not os.path.exists(sdir): os.mkdir(sdir)
    if not os.path.exists('figs'): os.mkdir('figs')
    if sname is None: sname=f'npz/usgs_flow_{y1}_{y2}'
    if os.path.exists(sname+'.npz') and not reDownload: return
    for m,station in enumerate(stations):
        #get links
        urls=['https://nwis.waterdata.usgs.gov/usa/nwis/uv/?cb_00060=on&format=rdb&site_no='+station+'&period=&begin_date='+StartT+'&end_date='+EndT,  #15-min flow
              'https://waterdata.usgs.gov/nwis/dv?cb_00060=on&format=rdb&site_no='+station+'&period=&begin_date='+StartT+'&end_date='+EndT] #daily flow
        tags=['15min','daily']
    
        #download usgs flow data; first, try 15min data; if fails, then, try daily data
        for url,tag in zip(urls,tags):
            fname='{}/{}_{}_{}_{}.txt'.format(sdir,station,y1,y2,tag) if y1!=y2 else '{}/{}_{}_{}.txt'.format(sdir,station,y1,tag)
            if fexist(fname): print('exist '+fname); break
            print('{}/{} download usgs flow: {} for {}-{} {} to {}'.format(m+1,len(stations),station,y1,y2,tag,fname))
            urlsave(url,fname)
            if station=='08072000': urlsave(url.replace('00060','00065'),fname) #download the water level for Lake Houston
            if os.path.getsize(fname)<3000: print('No data available');os.remove(fname)
            if fexist(fname): break  #if 15 min is available, hourly data will not be downloaded
    
    if os.path.exists(sname+'.npz') and not reRead:
        print('already exist ' +sname + '.npz')
    else:
        fnames=array(['{}/{}'.format(sdir,i) for i in os.listdir(sdir) if i.endswith('.txt') and i[0:8] in stations])
        station=[]; mtime=[]; flow=[]; river=[]
        for m,fname in enumerate(fnames):
            if os.path.getsize(fname)<3000: continue
            print('reading {}: {}/{}'.format(fname,m+1,len(fnames)))
            if os.path.exists(fname.replace('txt','npz')): #corresponding npz file, in case error occur for some specific stations
                print('loading existing '+fname.replace('txt','npz'))
                Z=loadz(fname.replace('txt','npz'))
                stationi,mti,flowi=Z.stationi,Z.mti,Z.flowi
            else:
                if 'daily' in fname:
                    stationi,mti,flowi=array([array(i.split('\t')[1:4]) for i in open(fname,'r').readlines() if i.startswith('USGS')]).T
                    #compute time
                    mti=datenum(mti)
                else:
                    stationi,mti,tz,flowi=array([array(i.split('\t')[1:5]) for i in open(fname,'r').readlines() if i.startswith('USGS')]).T
                    #compute time
                    mti=datenum(mti)
                St=zdata()
                St.stationi=stationi
                St.mti=mti
                St.flowi=flowi
                savez(fname.replace('txt','npz'),St)
            #save variables
            station.extend(stationi)
            mtime.extend(mti)
            flow.extend(flowi)
        #save data
        flow=array(flow); fpn=~((flow=='')|(flow=='Ssn')|(flow=='Eqp')|(flow=='***')|(flow=='Mnt'))
        S=zdata()
        S.station=array(station)[fpn]
        S.time=array(mtime)[fpn]
        S.flow=array(flow)[fpn].astype('float')*0.028316847
          
        #calculate flow for San Jacinto River based on wl-flow curve
        fpn=S.station=='08072000'
        if sum(fpn)>0:
            print('Get San Jacinto River discharge based on water level of Lake Houston')
            wl=S.flow[fpn]/0.028316847
            wtime=S.time[fpn]
            #-- get flow usin the wl-q curve
            fp=wtime>=datenum(2009,10,1)
            wl[fp]=wl[fp]+2.77
            tflow=86.76*(wl-44.5)/0.3
            fp=wl>44.8
            tflow[fp] = 5418.31*(wl[fp]-44.5)**2 -1366.29*(wl[fp]-44.5);
            tflow = tflow*0.028316847
            tflow[tflow<0]=0
            #-- update the flow
            print('Update the flow for San Jacinto River')
            S.flow[fpn]=tflow
        
        # save the data into npz format
        savez(sname,S)
def get_usgs_temp(stations=None,StartT='1980-1-1',EndT='2022-1-1',sname=None,reRead=False,sdir=None):
    y1=num2date(datenum(StartT)).year; y2=num2date(datenum(EndT)-1/24).year
    if sdir is None: sdir=f'usgs_temp_{y1}_{y2}'
    if not os.path.exists(sdir): os.mkdir(sdir)
    if not os.path.exists('figs'): os.mkdir('figs')
    if sname is None: sname=f'npz/usgs_temp_{y1}_{y2}'
    
    for m,station in enumerate(stations):
        #get links
        urls=['https://nwis.waterdata.usgs.gov/usa/nwis/uv/?cb_00010=on&format=rdb&site_no='+station+'&period=&begin_date='+StartT+'&end_date='+EndT,  #15-min flow
              'https://waterdata.usgs.gov/nwis/dv?cb_00010=on&format=rdb&site_no='+station+'&period=&begin_date='+StartT+'&end_date='+EndT] #daily flow
        tags=['15min','daily']
    
        #download usgs temp data; first, try 15min data; if fails, then, try daily data
        for url,tag in zip(urls,tags):
            fname='{}/{}_{}_{}_{}.txt'.format(sdir,station,y1,y2,tag) if y1!=y2 else '{}/{}_{}_{}.txt'.format(sdir,station,y1,tag)
            if fexist(fname): break
            print('download usgs temp: {} for {}-{} {}'.format(station,y1,y2,tag))
            urlsave(url,fname)
            if os.path.getsize(fname)<3000: print('No data available');os.remove(fname)
            if fexist(fname): break  #if 15 min is available, hourly data will not be downloaded
    
    if os.path.exists(sname+'.npz') and not reRead:
        print('already exist ' +sname + '.npz')
    else:
        fnames=array(['{}/{}'.format(sdir,i) for i in os.listdir(sdir) if i.endswith('.txt') and i[0:8] in stations])
        station=[]; mtime=[]; temp=[]; river=[]
        for m,fname in enumerate(fnames):
            if os.path.getsize(fname)<3000: continue
            print('reading {}: {}/{}'.format(fname,m+1,len(fnames)))
            if os.path.exists(fname.replace('txt','npz')): #corresponding npz file, in case error occur for some specific stations
                print('loading existing '+fname.replace('txt','npz'))
                Z=loadz(fname.replace('txt','npz'))
                stationi,mti,tempi=Z.station,Z.time,Z.temp
            else:
                if 'daily' in fname:
                    stationi,mti,tempi=array([array(i.split('\t')[1:4]) for i in open(fname,'r').readlines() if i.startswith('USGS')]).T
                    #compute time
                    mti=datenum(mti)
                else:
                    stationi,mti,tz,tempi=array([array(i.split('\t')[1:5]) for i in open(fname,'r').readlines() if i.startswith('USGS')]).T
                    #compute time
                    mti=datenum(mti)
                atemp=array(tempi); fpn=~((atemp=='')|(atemp=='Ssn')|(atemp=='Eqp')|(atemp=='***')|(atemp=='Mnt'))
                St=zdata()
                St.station=array(stationi)[fpn]
                St.time=array(mti)[fpn]
                St.temp=array(atemp)[fpn].astype('float')
                savez(fname.replace('txt','npz'),St)
                stationi,mti,tempi=St.station,St.time, St.temp #temp alareayd in float format
            #save variables
            station.extend(stationi)
            mtime.extend(mti)
            temp.extend(tempi)
            print('record ',len(station))
        #save data
        S=zdata()
        S.station=array(station)
        S.time=array(mtime)
        S.temp=array(temp).astype('float')
        # save the data into npz format
        print('save data into '+sname+'.npz')
        savez(sname,S)
        

def fill_usgs_flow_v2(S=None,begin_time=datenum(2007,1,1),end_time=datenum(2009,1,1),sname=None,recal=False):
    '''
    Get gap gree usgs flow data, use ltm to interp missing one
    Parameters
    ----------
    S : Zdata format structured data

    Returns
    -------
    newdaily

    '''
    # if file already exist
    if os.path.isfile(sname+'.npz') and not recal:
        print('already exist ',sname+'.npz')
        Z=loadz(sname+'.npz')
        return Z
    
    S.flow[S.flow<0]=0
    stations=unique(S.station)
    # get daily mean first
    dailydata={}
    ostation=[]
    otime=[]
    oflow=[]
    times=arange(begin_time,end_time)
    for station1 in stations:
        print('get daily mean for '+station1)
        fp=S.station==station1
        time1=S.time[fp]
        flow1=S.flow[fp]
        ltm=nanmean(flow1)
        time1,flow1=get_daily_mean(time1,flow1)
        dailydata[station1]=[time1,flow1]
        fp=(time1>=times.min())*(time1<=times.max())
        if sum(fp)==len(times):
            ostation.extend(tile(station1,len(times)))
            otime.extend(times)
            oflow.extend(flow1[fp])
        if sum(fp)!=len(times):
            flows=interp(times,time1,flow1,left=ltm,right=ltm)
            ostation.extend(tile(station1,len(times)))
            otime.extend(times)
            oflow.extend(flows)
    Z=zdata()
    Z.station=array(ostation)
    Z.flow=array(oflow)
    Z.time=array(otime)
    if not sname is None: 
        savez(sname,Z)
    return Z

def fill_usgs_flow(S=None,begin_time=datenum(2007,1,1),end_time=datenum(2009,1,1),sname=None,recal=False):
    '''
    Parameters
    ----------
    S : Zdata format structured data

    Returns
    -------
    newdaily

    '''
    # if file already exist
    if os.path.isfile(sname+'.npz') and not recal:
        print('already exist ',sname+'.npz')
        Z=loadz(sname+'.npz')
        return Z
    
    S.flow[S.flow<0]=0
    stations=unique(S.station)
    # get daily mean first
    dailydata={}
    times=arange(begin_time,end_time)
    for station1 in stations:
        print('get daily mean for '+station1)
        fp=S.station==station1
        time1=S.time[fp]
        flow1=S.flow[fp]
        time1,flow1=get_daily_mean(time1,flow1) #daily mean
        dailydata[station1]=[time1,flow1]

    #- get gap-free daily data
    newdaily={}
    ostation=[]
    otime=[]
    oflow=[]
    times=arange(begin_time,end_time)
    for station1 in stations:
        time1,flow1=dailydata[station1]
        ltm=nanmean(flow1)
        #find the correlation between stations
        B=[]
        BR=[]
        for station2 in stations:
            if station2==station1: continue
            time2,flow2=dailydata[station2]
            best=linear_lag_analysis(time1,flow1,time2,flow2,daily=True,shiftings=arange(-5,5)) #return shifting, ratio, linear coefficient a and b
            best.station=station2
            B.append(best)
            BR.append(best.R)
            #print('get lag correlation for '+station1+' vs '+station2 +': {:.2f}'.format(best.R))
        #print('best R: {:.2f}'.format(R))
        
        #- for the missing values, interp from best related stations
        time1=list(time1) #convert to list to make it extentable
        flow1=list(flow1) #convert to list to make it extentable
        for ist in argsort(BR)[::-1]:
            print('stations and R: ',station1,B[ist].station,BR[ist])
            missing_time=find_missing_time(time1, times)
            if BR[ist]<0.5: break
            ratio=B[ist].pair1.mean()/B[ist].pair2.mean()
            #slope, intercept, r, p, std_err = stats.linregress(B[ist].pair2, B[ist].pair1)
            time2,flow2=dailydata[B[ist].station]
            for imtime in missing_time:
                if imtime in time2:
                    fp=time2==imtime
                    time1.append(imtime)
                    flow1.append(ratio*flow2[fp])
            
            #check after filing up
            print('missing record after filling: {}'.format(len(find_missing_time(time1, times))))
            if len(find_missing_time(time1, times))==0: break
        time1=array(time1)
        flow1=array(flow1).astype('float')
        ind=argsort(time1)
        time1=time1[ind]
        flow1=flow1[ind]
        fp=(time1>=times.min())*(time1<=times.max())
        if sum(fp)-len(times)==0: flow1=flow1[fp]
        if len(times)-sum(fp)>0 and len(times)-sum(fp)<30: flow1=interp(times,time1,flow1) #interp 
        if len(times)-sum(fp)>=30: #fill with the long-term mean
            time1=list(time1) #convert to list to make it extentable
            flow1=list(flow1) #convert to list to make it extentable
            missing_time=find_missing_time(time1, times)
            time1.extend(missing_time)
            flow1.extend(tile(ltm,len(missing_time)))
            time1=array(time1)
            flow1=array(flow1).astype('float')
            ind=argsort(time1)
            time1=time1[ind]
            flow1=flow1[ind]
            fp=(time1>=times.min())*(time1<=times.max())
            flow1=flow1[fp]
        if len(times)!=len(flow1): 
            print('len(times),len(flow1)',len(times),len(flow1))
            sys.exit('dimension of flow does not match with the time')
        newdaily[station1]=[times,flow1]
        ostation.extend(tile(station1,len(times)))
        otime.extend(times)
        oflow.extend(flow1)
    Z=zdata()
    Z.station=array(ostation)
    Z.flow=array(oflow)
    Z.time=array(otime)
    if not sname is None: 
        savez(sname,Z)
    return Z

# some of my own public functions
def get_noaa_tide_current(stations=['8637689'],years=arange(2007,2022),varnames=['hourly_height'],sdir='data',load_again=False):
    print('== get data from noaa tide & current ==' )
    if not os.path.exists(sdir): os.mkdir(sdir)
    #url0='https://api.tidesandcurrents.noaa.gov/api/prod/datagetter?units=metric&time_zone=gmt&application=NCCOOS&format=csv&interval=h&datum=msl'
    for station in stations:
        # get lon lat from json txt output; optional
        url='https://api.tidesandcurrents.noaa.gov/mdapi/prod/webapi/stations/{}.json?expand=details,products&units=metric'.format(station)
        try: 
            urlsave(url,'tmp.json')
        except:
            print('no data found for station ',station)
            continue
        with open('tmp.json') as f: jdata = json.load(f)
        lat,lon,name=jdata['stations'][0]['lat'],jdata['stations'][0]['lng'],jdata['stations'][0]['name']
        os.remove('tmp.json')
        print('{}:{}, {}N, {}W'.format(station,name,lon,-lat))
        
        # get data
        print('get noaa data')
        for year in years:
            for product_name in varnames:
                if product_name=='currents':
                    url0='https://api.tidesandcurrents.noaa.gov/api/prod/datagetter?units=metric&time_zone=gmt&application=NCCOOS&format=csv&interval=h&datum=msl&vel_type=default'
                else:
                    url0='https://api.tidesandcurrents.noaa.gov/api/prod/datagetter?units=metric&time_zone=gmt&application=NCCOOS&format=csv&interval=h&datum=msl'
                url='{}&begin_date={}0101&end_date={}1231&station={}&product={}'.format(url0,year,year,station,product_name)
                #print(url)
                fname='{}/{}_{}_{}.csv'.format(sdir,product_name,station,year)
                if os.path.isfile(fname) and not load_again: continue
                try: #when there is no data at the given station, error may occurs. So use "try"
                    urlsave(url,fname)
                    print('.. save into {}'.format(fname))
                except:
                    print('... bad request for {} (likely data do not exist)'.format(product_name))
                    pass
            
    ''' Some Notes
    1. interval=h is necesary when download period between one month and a year
    2. for each api option see: https://api.tidesandcurrents.noaa.gov/api/prod/
    3. api helper in https://tidesandcurrents.noaa.gov/api-helper/url-generator.html
    '''

def process_noaa_tide_current(stations=['8637689'],years=arange(2007,2022),varnames=['hourly_height'],sname_pre='npz/',sdir='data/',read_again=False):
    print('== read noaa tide and current data ==')
    if not sdir.endswith('/'): sdir=sdir+'/'
    for varname in varnames: #['wind','hourly_height','water_temp','air_temp','conductivity']:
        if len(stations)==1: sname='{}{}_{}_{}_{}'.format(sname_pre,station,varname,years[0],years[-1])
        if len(stations)>1: sname='{}{}_{}_{}'.format(sname_pre,varname,years[0],years[-1])
        if os.path.isfile(sname+'.npz') and not read_again: print('file exist '+sname+'.npz'); continue
        time=[] #initialize the list, list is extentable
        data=[]
        station=[]
        if varname=='wind': data2=[] #for wind direction
        for stid in stations:
            for year in years:
                fname='{}{}_{}_{}.csv'.format(sdir,varname,stid,year)
                print('reading '+fname)
                if not os.path.isfile(fname) or os.path.getsize(fname)<1e3:continue #if not exist or file size less than 1kb
                tdata=array([i.split(',') for i in open(fname,'r').readlines() if i.startswith(str(year))])
                time.extend(datenum(tdata[:,0]))
                data.extend(tdata[:,1]) 
                if varname=='wind': data2.extend(tdata[:,2]) 
                station.extend(list(full([len(tdata),],stid)))
            
        S=zdata()
        
        # remove nan values
        station=array(station)
        data=array(data)
        data[data=='']='NaN'
        data=data.astype('float')
        fp=[data!=nan]
        data=data[fp]
        S.station=array(station)[fp]
        S.time=array(time)[fp]
        
        #remove spiking values within a given station
        if varname=='water_temperature':
            print(f'remove spiking vlaue for {varname}')
            for m,value in enumerate(data):
                if m==0: continue
                if abs(value-data[m-1])>5 and station[m]==station[m-1]:
                    data[m]='NaN'                    
        if varname=='wind':
            data2=array(data2)
            data2[data2=='']='NaN'
            S.wdir=data2.astype('float')
            S.speed=data.astype('float')
        elif varname=='water_temperature':
            S.temp=data.astype('float')
        elif varname=='hourly_height':
            S.wl=data.astype('float')
        else:
            S.data=data.astype('float')
        savez(sname,S)
        print(f'data saved into {sname}.npz')



#%% Data processing

def get_daily_mean(times,flow,lpf=False):
    #get daily mean
    if lpf:
        print('low-pass filter with cut-off frequency of 1/(7*24)')
        flow=lpfilt(flow,median(diff(times))*24,1/(7*24))
    # use low pass filter 
    print('in get daily mean')
    nrec=0
    sumval=0
    preday=floor(times[0])
    newtime,newflow=[],[]
    for itime,iflow in zip(times,flow):
        if itime<preday+1:
            nrec+=1
            sumval+=iflow
        else:
            newtime.append(preday)
            newflow.append(sumval/nrec)
            preday=floor(itime)
            nrec=1
            sumval=iflow
    newtime=array(newtime)+0.5
    newflow=array(newflow)
    return newtime,newflow

def remove_spike(times,data,ds=2,hw=5,rmnan=True,inter_spike=False):
    #ds: delta standard deviation
    #hw: half window
    #inter_spike: replace spiking values with the mean over [hw,hw]
    newdata=data.copy()
    newtime=times.copy()
    for m,itime in enumerate(times):
        fp=(abs(times-itime)<=hw)*(times!=itime)
        if sum(fp)<hw: continue
        tmp=data[fp]
        idata=data[m]
        if abs(idata-tmp.mean())>ds*tmp.std():
            newdata[m]=nan
            if inter_spike: #there will be no nan value
                newdata[m]=tmp.mean()
    if rmnan: #remove nan
        fp=isnan(newdata)
        newdata=newdata[~fp]
        newtime=newtime[~fp]
    return newtime,newdata
def linear_lag_analysis(time1,flow1,time2,flow2,daily=True,shiftings=arange(-5,6),degree=4):
    from sklearn.linear_model import LinearRegression
    from sklearn.preprocessing import PolynomialFeatures
    print('in linear_lag_analysis')
    if daily and median(diff(time1))!=1: time1,flow1=get_daily_mean(time1,flow1)
    if daily and median(diff(time2))!=1: time2,flow2=get_daily_mean(time2,flow2)
    best=zdata()
    best.R=0
    time0=time2.copy()
    best.pair2,best.pair1=[],[]
    best.r2=0
    for ishifting in shiftings:
        time2=time0+ishifting
        # find the pair
        pair1,pair2=[],[]
        iloc=0
        for m,itime1 in enumerate(time1):
            while(iloc<len(time2) and time2[iloc]<itime1):
                iloc+=1
            if iloc>=len(time2): break
            if time2[iloc]==itime1 and not isnan(flow1[m]) and not isnan(flow2[iloc]):
                pair1.append(flow1[m])
                pair2.append(flow2[iloc])
        pair1,pair2=array(pair1),array(pair2)
        if len(pair2)==0: print('no paired data')
        S=get_stat(pair1,pair2)
        print(S.R)
        if S.R>best.R:
            print(S.R,ishifting)
            best.R=S.R
            best.shifting=ishifting
            best.pair1=pair1
            best.pair2=pair2
    if len(best.pair2)==0: return best
    x=best.pair2.reshape(-1,1)
    y=best.pair1
    x_pred=arange(min(x),max(x)).reshape(-1,1)
    transformer = PolynomialFeatures(degree=degree, include_bias=False)
    transformer.fit(x)
    x_ = transformer.transform(x)
    model = LinearRegression().fit(x_, y)
    r_sq = model.score(x_, y)
    y_pred = model.predict(transformer.transform(x_pred))
    best.transformer=transformer
    best.model=model
    best.r2=r_sq
    best.x_pred=x_pred
    best.y_pred=y_pred
    return best
def seasonal_impose(time1,flow1,mtime,mmean,mstd,mmin,mmax,ds=None):
    #limit values to the lower and upper bound based on seasonal cycle
    #ds: if ds is not None, lower and upper bound will be mean+-ds*std; otherwise the lowe and upper will be min and max at each month
    y1=num2date(time1.min()).year
    y2=num2date(time2.max()).year
    ip=0
    flow2=flow1.copy()
    mtime2,mmean2,mstd2,mmin2,mmax2=[],[],[],[],[]
    for iyear in arange(y1,y2+1):
        for imon in arange(1,13):
            mtime2.append(datenum(iyear,imon,15))
            mmean2.append(mmean[imon-1])
            mstd2.append(mstd[imon-1])
            mmax2.append(mmax[imon-1])
            mmin2.append(mmin[imon-1])
    mmean3=interp(time1,mtime2,mmean2)
    mstd3=interp(time1,mtime2,mstd2)
    mmin3=interp(time1,mtime2,mmin2)
    mmax3=interp(time1,mtime2,mmax2)
    for m,flowi in enumerate(flow1):
        flow2[m]=min(mmax3[m],max(flowi,mmin3[m]))
        if not ds is None:
            flow2[m]=min(mmean3[m]+ds*mstd3[m],max(flowi,mmean3[m]-ds*mstd3[m]))
    if not ds is None:
        mlower=mmean3-ds*mstd3
        mupper=mmean3+ds*mstd3
    else:
        mlower=mmin3
        mupper=mmax3
    return time1,flow2,mlower,mupper
            
def find_missing_time(time1,time2,dt=0.1,block=0):
    #find those in time2 but not in time1
    #dt: not used noww
    #block: if block is not zero, only missing value within a continuous block will be regarded as missing time; 
    #       
    missing_time=[]
    for itime in time2:
        #if sum(abs(itime-time1)<dt)==0:
        if not itime in time1:
            missing_time.append(itime)
    #only return those missing one for a continuous block; remove thos scatter ones
    
    if block>0:
        difft=time2[1]-time2[0]
        if difft<0: ppppp
        ntime=missing_time.copy()
        for m,itime in enumerate(missing_time):
            if m==0: ibegin=0; iend=0; continue
            if abs(itime-missing_time[m-1]-difft)<0.1*difft:
                iend=m
            else:
                if iend-ibegin<block:
                    del ntime[ibegin:iend+1]
                ibegin=m
        print('before and after removing scattered one: {} {}'.format(len(missing_time),len(ntime)))
        missing_time=ntime
    
    return missing_time

def mean_profile(ttime,tdep,tdata):
    '''
    Find the mean profile given all availabe records 
    '''
    utime=unique(ttime)
    ttime,tdep,tdata=array(ttime),array(tdep),array(tdata)
    pdep=linspace(tdep.min(),tdep.max(),20)
    z=[]
    for timei in utime:
        fp=ttime==timei
        tmpdep=tdep[fp]
        tmpdata=tdata[fp]
        # avearge those at same depth
        udep=unique(tmpdep)
        udata=[nanmean(tmpdata[tmpdep==depi]) for depi in udep]
        z.append(interp(pdep,udep,udata))
    z=array(z)
    pdata=z.mean(axis=0)
    return pdep,pdata
#%% plotting related
def set_xtick(fmt=0):
    xlims=gca().get_xlim()
    begin_year=num2date(xlims[0]).year
    end_year=num2date(xlims[1]).year
    if begin_year==end_year or fmt==1:
        xts,xls=get_xtick(fmt=1) #get monthly ticks
    else:
        xts,xls=get_xtick() #get yearly tickes  
    setp(gca(),xticks=xts[0:-1:1],xticklabels=xls[0:-1:1],xlim=xlims)

def plot_vgrid(vgrid='vgrid.in',bname='transect.bp',hgrid='hgrid.gr3'):
    vd=read_schism_vgrid(vgrid)
    gd=read_schism_hgrid(hgrid)
    z=vd.compute_zcor(gd.dp)
    print('plot transect for '+bname)
    bp=read_schism_bpfile(str(bname))
    #compute dist
    dist=[0,]
    for i in arange(bp.nsta-1):
        disti=abs((bp.x[i+1]-bp.x[i])+1j*(bp.y[i+1]-bp.y[i]))+dist[i]
        dist.append(disti)
    dist=array(dist)

    sindp=near_pts(c_[bp.x,bp.y],c_[gd.x,gd.y]);
    zi=z[sindp]

    figure(figsize=[8,6])
    for k in arange(vd.nvrt):
        plot(dist,zi[:,k],'k-',lw=0.5)
    for i in arange(len(bp.x)):
        plot(ones(vd.nvrt)*dist[i],zi[i,:],'k-',lw=0.5)
    setp(gca(),ylim=[zi.min()-1,0.0],xlim=[0,dist.max()])
    gcf().tight_layout()
    # move_figure(gcf(),0,0)
    savefig('{}_vgrid'.format(bname.replace('.bp','')))

def plot_section(x,y,z,xn=40,yn=20,vmin=0,vmax=30,levels=None,cmap='jet',show_datapoint=False):
    if levels is None: levels=arange(vmin,vmax+1,1)
    xi = linspace(min(x), max(x), xn)
    yi = linspace(min(y), max(y), yn)
    
    #-- add more data points toward the bottom
    x,y,z=array(x),array(y),array(z)
    ux=unique(x)
    x2,y2,z2=[],[],[]
    for tx in ux:
        fp=x==tx
        x2.extend(tile(tx,len(yi)))
        y2.extend(yi)
        z2.extend(interp(yi,y[fp],z[fp]))
    x2,y2,z2=array(x2),array(y2),array(z2)
    #x,y,z=array(x),array(y),array(z)
    xi = linspace(min(x), max(x), xn)
    yi = linspace(min(y), max(y), yn)
    Z = griddata((x2, y2), z2, (xi[None,:], yi[:,None]), method='linear')
    X, Y = meshgrid(xi, yi)
    Z[Z>vmax]=vmax
    contourf(xi,yi,Z,levels=levels,vmin=vmin,vmax=vmax,cmap=cmap,linestyles='solid')
    if show_datapoint: plot(x,y,'x',ms=1)
    
    #fill the land at the bottom
    ux=unique(x)
    uy=[]
    for xi in ux:
        uy.append(y[x==xi].min())
    uy=array(uy)
    ux=[ux.max(),ux.min(),*ux]
    uy=[uy.min(),uy.min(),*uy]
    fill(ux,uy,color='w')

def plot_2d(run='run02j',stack=1,itime=0,ilay=51,grid='../run02j/hgrid.gr3',var='totalSuspendedLoad',cmap='jet',
            saveplot=False,showplot=True,odir='../run02j/outupts',clim=[0,30],use_bottom=True):
    if not os.path.exists(f'{run}/figs'): os.mkdir(f'{run}/figs')
    if grid.endswith('.npz'):
        gd=loadz(grid).hgrid
    else:
        gd=read_schism_hgrid(grid)

    # get the data first
    C=ReadNC(f'{odir}/out2d_{stack}.nc',1); vars_2d=[*C.variables]; C.close() #get 2d output variables
    if var in vars_2d:
        value= C[var][itime].ravel()
    else:
        if not os.path.exists(f'{odir}/{var}_{stack}.nc'): print(f'no output for {var}');
        try:
            C=ReadNC(f'{odir}/{var}_{stack}.nc',1);
            value= C[var][itime,:,ilay].ravel()
            if use_bottom and sum(value>1e5)>0: value0=C[var][itime,:,:] #entire column
            if var.find('sedConc')!=-1 or var=='totalSuspendedLoad': value=value*1000
        except:
            pass
    if var=='totalSuspendedLoad' and not os.path.exists(f'{odir}/{var}_{stack}.nc'):
        value=zeros(len(gd.x))
        if use_bottom: value0=zeros([len(gd.x),52])
        for var in ['sedConcentration_1','sedConcentration_2','sedConcentration_3','sedConcentration_4']:
            C=ReadNC(f'{odir}/{var}_{stack}.nc',1);
            value += C[var][itime,:,ilay].ravel()*1000
            if sum(value>1e5)>0 and use_bottom: value0 += C[var][itime,:,:]*1000

    value[value>1e5]=nan
    if use_bottom and sum(isnan(value))>0:
        print('use the bottom value')
        #for nan values, use the bottom one
        value0[value0>1e5]=nan
        fp=nonzero(isnan(value))[0]
        for ig in fp:
            tmp=value0[ig,:]
            ind=nonzero(~isnan(tmp))[0]
            value[ig]=tmp[ind[0]]

    if len(clim)==0: clim=[percentile(value[~isnan(value)],5),percentile(value[~isnan(value)],95)]
    gd.plot(fmt=1,value=value,clim=clim,cmap=cmap)
    gd.plot_bnd()
    show()
    if saveplot: savefig(f'{run}/figs/{var}_{stack}_{itime}_{ilay}.png')

