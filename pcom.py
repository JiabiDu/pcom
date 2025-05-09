#!/usr/bin/env python3

from pylib import *
import json
from scipy.stats import gaussian_kde
from scipy.interpolate import griddata
from math import atan
import imageio
#================================================================
# AUTHOR NOTES
# This package is used for Coastal Ocean Modelling
# It has several components: data doaloading, data processing, plotting 
# contact jiabi.du@gmail.com; jdu@tamug.edu
# our website: www.coastalmodeling.org
#===============================================================

#%% schism related
def in_domain(lons,lats,hgrid):
    lons,lats=array(lons),array(lats)
    if hgrid.endswith('npz'):
        gd=loadz(hgrid).hgrid
        gd.x,gd.y=gd.lon,gd.lat
    else:
        gd=read_schism_hgrid(hgrid)
    ie,ip,acor=gd.compute_acor(c_[lons,lats])
    return ie!=-1 #points inside model modmain

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
def pair_data(x1,y1,x2,y2,hw=0.2/24,remove_nan=True):
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
        if remove_nan and (isnan(iy) or isnan(nanmean(y2[fp]))): continue
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
def get_cbp_data(years=[1984,2022],stations=[1325,1328],pstr='63,78',sdir='data/cbp_nutrient/',sname='data/cbp_nutrient.npz',download_again=False):
    #------------------------------------------------------------------------------
    #download data
    #------------------------------------------------------------------------------
    nowtime=datetime.datetime.now().strftime("%Y-%m-%d")
    url0='http://data.chesapeakebay.net/api.CSV/WaterQuality/WaterQuality/1-1-{}/12-31-{}/0,1/2,4,6/12,13,15,35,36,2,3,7,33,34,23,24/Station'.format(*years)
    if not os.path.exists(sdir): os.mkdir(sdir)
    #download data for each station
    for m,station in enumerate(stations):
        url='{}/{}/{}'.format(url0,station,pstr)
        fname='{}/{}_{}_{}.csv'.format(sdir,station,years[0],years[-1])
        if os.path.exists(fname) and os.path.getmtime(fname)>=time.time()-3600*6 and not download_again: print(fname,'exists and updated'); continue
        print('download: {}, {}/{}'.format(station,m+1,len(stations)))
        try:
            urlsave(url,fname)
        except:
            print(f'no data found for {station}')
            pass
        if os.path.getsize(fname)<500: os.remove(fname)
    
    #------------------------------------------------------------------------------
    #read data
    #------------------------------------------------------------------------------
    fnames=array([sdir+i for i in os.listdir(sdir) if i.endswith('.csv')]) #get all data files
    if updated(sname,fnames): print(sname,'exists and updated'); return
    svars=array(['Station','SampleDate','SampleTime','TotalDepth','Depth','Layer','Parameter','MeasureValue','Unit','Latitude','Longitude'])
    mvars=array(['station','date','time','tdepth','depth','layer','var','data','unit','lat','lon'])
    S=zdata(); 
    for i in mvars: exec('S.{}=[]'.format(i))
    for m,fname in enumerate(fnames):
        print('reading: {}, {}/{}'.format(fname,m,len(fnames)))
        #read all fields
        try:
            lines=array([i.strip().replace('"','').split(',') for i in open(fname,'r').readlines() if len(i)>50])
        except:
            lines=array([i.strip().replace('"','').split(',') for i in open(fname,'r',encoding='cp850').readlines() if len(i)>50])
            
        #remove invalid line
        tl=lines[0]; fp=array([nonzero(tl==i)[0][0] for i in svars]); lines=lines[1:,fp]
        fpn=lines[:,7]!=''; lines=lines[fpn]
    
        #assign -999 to empty field
        fpn=lines[:,3]==''; lines[fpn,3]='-999'
        fpn=lines[:,4]==''; lines[fpn,4]='-999'
    
        #save each variables
        for n,mvar in enumerate(mvars): exec('S.{}.extend(lines[:,{}])'.format(mvar,n))
    
    #orgnaize data, change format
    for i in mvars: exec('S.{}=array(S.{})'.format(i,i))
    S.tdepth=S.tdepth.astype('float32'); S.depth=S.depth.astype('float32'); S.data=S.data.astype('float')
    S.lat=S.lat.astype('float'); S.lon=S.lon.astype('float')
    
    #get time
    S.time=datestr2num(array(['{} {}'.format(i,j) for i,j in zip(S.date,S.time)])); del S.date
    
    #try to shrink data
    uvars=unique(S.var); S.unit=dict(zip(uvars, array([unique(S.unit[S.var==i])[0] for i in uvars])))
    stations=unique(S.station); S.lon=dict(zip(stations, squeeze(array([unique(S.lon[S.station==i]) for i in stations]))))
    stations=unique(S.station); S.lat=dict(zip(stations, squeeze(array([unique(S.lat[S.station==i]) for i in stations]))))
    
    #save data
    savez(sname,S)

def get_usgs_flow(stations=None,StartT='1980-1-1',EndT='2022-1-1',sname=None,reRead=False,sdir=None, reDownload=False):
    y1=num2date(datenum(StartT)).year; y2=num2date(datenum(EndT)-1/24).year
    if sdir is None: sdir=f'usgs_{y1}_{y2}'
    if not os.path.exists(sdir): os.mkdir(sdir)
    if not os.path.exists('figs'): os.mkdir('figs')
    if sname is None: 
        sname=f'npz/usgs_flow_{y1}_{y2}'
        if not fexist('npz'): os.mkdir('npz')
    if os.path.exists(sname+'.npz') and not reDownload: return
    for m,station in enumerate(stations):
        #get links
        urls=['https://nwis.waterdata.usgs.gov/usa/nwis/uv/?cb_00060=on&format=rdb&site_no='+station+'&period=&begin_date='+StartT+'&end_date='+EndT,  #15-min flow
              'https://waterdata.usgs.gov/nwis/dv?cb_00060=on&format=rdb&site_no='+station+'&period=&begin_date='+StartT+'&end_date='+EndT] #daily flow
        tags=['daily','15min']
    
        #download usgs flow data; first, try 15min data; if fails, then, try daily data
        for url,tag in zip(urls,tags):
            fname='{}/{}_{}_{}_{}.txt'.format(sdir,station,y1,y2,tag) if y1!=y2 else '{}/{}_{}_{}.txt'.format(sdir,station,y1,tag)
            if fexist(fname) and not reDownload: print('exist '+fname); break
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
            if os.path.getsize(fname)<6000: continue
            if not reRead and os.path.exists(fname.replace('txt','npz')): #corresponding npz file, in case error occur for some specific stations
                print('loading existing '+fname.replace('txt','npz'))
                Z=loadz(fname.replace('txt','npz'))
                stationi,mti,flowi=Z.stationi,Z.mti,Z.flowi
                if type(flowi[0])==str or type(flowi[0])==str_:
                    fpn=array([set(i).issubset({'0','1','2','3','4','5','6','7','8','9','.'}) and len(i)!=0 for i in flowi])
                    stationi,mti,flowi=stationi[fpn],mti[fpn],flowi[fpn].astype('float')
            else:
                print('reading {}: {}/{}'.format(fname,m+1,len(fnames)))
                if 'daily' in fname:
                    stationi,mti,flowi=array([array(i.split('\t')[1:4]) for i in open(fname,'r').readlines() if i.startswith('USGS')]).T
                    fpn=array([set(i).issubset({'0','1','2','3','4','5','6','7','8','9','.'}) and len(i)!=0 for i in flowi])
                    if sum(fpn)==0: 
                       print('format is wrong: not in daily format; will try to read in 15min format')
                       stationi,mti,tz,flowi=array([array(i.split('\t')[1:5]) for i in open(fname,'r').readlines() if i.startswith('USGS')]).T
                       fpn=array([set(i).issubset({'0','1','2','3','4','5','6','7','8','9','.'}) and len(i)!=0 for i in flowi])
                       if sum(fpn)==0: sys.exit('There is still no records for this file',fname)
                    #compute time
                    mti=quickdatenum(mti)
                else:
                    stationi,mti,tz,flowi=array([array(i.split('\t')[1:5]) for i in open(fname,'r').readlines() if i.startswith('USGS')]).T
                    #compute time
                    mti=quickdatenum(mti)
                #remove records containing non-number values
                fpn=array([set(i).issubset({'0','1','2','3','4','5','6','7','8','9','.'}) and len(i)!=0 for i in flowi])
                stationi,mti,flowi=stationi[fpn],mti[fpn],flowi[fpn].astype('float')
                St=zdata()
                St.stationi=stationi
                St.mti=mti
                St.flowi=flowi
                print('save data into',fname.replace('txt','npz'))
                savez(fname.replace('txt','npz'),St)
            #save variables
            station.extend(stationi)
            mtime.extend(mti)
            flow.extend(flowi.astype('float'))
            
        #save data
        flow=array(flow); 
        #fpn=~((flow=='Ice')|(flow=='')|(flow=='Ssn')|(flow=='Eqp')|(flow=='***')|(flow=='Mnt')|(flow=='Dis'))
        #fpn=array([set(i).issubset({'0','1','2','3','4','5','6','7','8','9','.'}) for i in flow])
        S=zdata()
        S.station=array(station)
        S.time=array(mtime)
        S.flow=array(flow).astype('float')*0.028316847
          
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

def delta_mtime(fname,days):
    # return if the modification time of fname is within a defined days from now
    timestamp=os.path.getmtime(fname)
    return datetime.datetime.now()-datetime.datetime.fromtimestamp(timestamp)>datetime.timedelta(days=days)

def get_usgs_temp(stations=None,StartT='1980-1-1',EndT='2022-1-1',sname=None,reRead=False,sdir=None,reLoad=True):
    y1=num2date(datenum(StartT)).year; y2=num2date(datenum(EndT)-1/24).year
    if sdir is None: sdir=f'usgs_temp_{y1}_{y2}'
    os.makedirs(sdir,exist_ok=True)
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
            if fexist(fname): 
               if delta_mtime(fname,2/24): break #downloaded within 2 hr from now
               if not reLoad: break
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
def get_noaa_tide_current(stations=['8637689'],years=arange(2007,2022),varnames=['hourly_height'],sdir='data',load_again=False,datum='navd'):
    print('== get data from noaa tide & current ==' )
    if not os.path.exists(sdir): os.mkdir(sdir)
    #url0='https://api.tidesandcurrents.noaa.gov/api/prod/datagetter?units=metric&time_zone=gmt&application=NCCOOS&format=csv&interval=h&datum=msl'
    for station in stations:
        ### get lon lat from json txt output; optional
        # url='https://api.tidesandcurrents.noaa.gov/mdapi/prod/webapi/stations/{}.json?expand=details,products&units=metric'.format(station)
        # try: 
        #     urlsave(url,'tmp.json')
        # except:
        #     print('no data found for station ',station)
        #     continue
        # with open('tmp.json') as f: jdata = json.load(f)
        # lat,lon,name=jdata['stations'][0]['lat'],jdata['stations'][0]['lng'],jdata['stations'][0]['name']
        # os.remove('tmp.json')
        # print('{}:{}, {}N, {}W'.format(station,name,lon,-lat))
        
        # get data
        print('get noaa data')
        for product_name in varnames:
            for year in years:
                if product_name=='currents':
                    url0=f'https://api.tidesandcurrents.noaa.gov/api/prod/datagetter?units=metric&time_zone=gmt&application=NCCOOS&format=csv&interval=h&datum={datum}&vel_type=default'
                else:
                    url0=f'https://api.tidesandcurrents.noaa.gov/api/prod/datagetter?units=metric&time_zone=gmt&application=NCCOOS&format=csv&interval=h&datum={datum}'
                url='{}&begin_date={}0101&end_date={}1231&station={}&product={}'.format(url0,year,year,station,product_name)
                #print(url)
                fname='{}/{}_{}_{}.csv'.format(sdir,product_name,station,year)
                if os.path.isfile(fname) and not load_again: print('exists',fname); continue
                try: #when there is no data at the given station, error may occurs. So use "try"
                    urlsave(url,fname)
                    print('.. save into {}'.format(fname))
                except:
                    print('... bad request for {} at {} (likely data do not exist)'.format(product_name,station))
                    break
                    pass
            
    ''' Some Notes
    1. interval=h is necesary when download period between one month and a year
    2. for each api option see: https://api.tidesandcurrents.noaa.gov/api/prod/
    3. api helper in https://tidesandcurrents.noaa.gov/api-helper/url-generator.html
    '''
def get_noaa_predict_tide(sname='data/noaa_ptide.npz',year=2022,stations=['8637689',],sdir='./',load_again=False):
    if os.path.exists(sname) and not load_again: print(sname,'exists'); return
    time,ptide,station=[],[],[]
    for stationi in stations:
        url=f'https://api.tidesandcurrents.noaa.gov/api/prod/datagetter?product=predictions&application=NOS.COOPS.TAC.WL&begin_date={year}0101&end_date={year}1231&datum=MSL&station={stationi}&time_zone=GMT&units=metric&interval=h&format=csv'
        fname=f'{sdir}/tide_prediction_{stationi}.csv'
        urlsave(url,fname)
        print('complete downloading data into',fname)
        data=[i.strip().split(',') for i in open(fname,'r').readlines()[1:]]; data=array(data)
        time.extend(datenum(data[:,0]))
        ptide.extend(data[:,1].astype('float'))
        station.extend(tile(stationi,len(data)))
    S=zdata()   
    S.time,S.ptide,S.station=array(time),array(ptide),array(station)
    savez(sname,S)
    
def process_noaa_tide_current(stations=['8637689'],years=arange(2007,2022),varnames=['hourly_height'],sname_pre='npz/',sdir='data/',read_again=False,sname=None):
    print('== read noaa tide and current data ==')
    if not sdir.endswith('/'): sdir=sdir+'/'
    for varname in varnames: #['wind','hourly_height','water_temp','air_temp','conductivity']:
        if len(stations)==1 and sname==None: sname='{}{}_{}_{}_{}'.format(sname_pre,stations[0],varname,years[0],years[-1])
        if len(stations)>1 and sname==None: sname='{}{}_{}_{}'.format(sname_pre,varname,years[0],years[-1])
        if os.path.isfile(sname+'.npz') and not read_again: print('file exist '+sname+'.npz'); continue
        time=[] #initialize the list, list is extentable
        data=[]
        station=[]
        if varname=='wind': data2=[] #for wind direction
        for m,stid in enumerate(stations):
            for year in years:
                fname='{}{}_{}_{}.csv'.format(sdir,varname,stid,year)
                print(m+1,'reading '+fname)
                if not os.path.isfile(fname) or os.path.getsize(fname)<1e3: print('no data'); continue #if not exist or file size less than 1kb
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
        fp=~isnan(data)
        data=data[fp]
        if varname=='wind': data2=array(data2); data2=data2[fp]
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
        os.makedirs(sname[0:sname.rfind('/')],exist_ok=True) #create the folder 
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
        pair1,pair2=array(pair1).ravel(),array(pair2).ravel()
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
    y2=num2date(time1.max()).year
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

def combine_data(dataA='data/noaa_tide.npz',dataB='data/noaa_ptide.npz',sname='data/tide.npz',varA=None,varB=None,add_diff=False,redo=False):
    if updated(sname,[dataA,dataB]) and not redo: print(sname+' exists and updated');return
    A=loadz(dataA); B=loadz(dataB)
    if not (hasattr(A,'time') and hasattr(A,'time')): sys.exit('time not exist in the input file')
    fp=B.time>A.time.max() #must have time
    if varA==None: #find the data data name automatically, suitable when there is just one variable
        atts=[var for var in vars(A)]
        varA=[i for i in atts if not i in ['station','time','VINFO']][0]
        print(f'{varA} in {dataA} to be combiend')
        #remove nan values
    #exec(f'tfp=~isnan(A.{varA}); A.{varA}=A.{varA}[tfp]; A.time=A.time[tfp]')
    if varB==None:
        atts=[var for var in vars(B)]
        varB=[i for i in atts if not i in ['station','time','VINFO']][0]
        print(f'{varB} in {dataB} to be combiend')
    #exec(f'tfp=~isnan(B.{varB}); B.{varB}=B.{varB}[tfp]; B.time=B.time[tfp]')
    if add_diff: #add a systematic difference in historical data
        T=zdata()
        x1=A.time; exec(f'T.y1=A.{varA}')
        x2=B.time; exec(f'T.y2=B.{varB}')
        y1,y2=pair_data(x1,T.y1,x2,T.y2,hw=0.2/24)
        mdiff=nanmean(y1)-nanmean(y2) #y1 is the usually the true value (eg. observation)
        exec(f'B.{varB}=B.{varB}+mdiff')
        print(f'!! add {mdiff} to {varB} in {dataB}; len y1 and y2=',len(y1),len(y2))
    S=zdata()
    exec('S.data=concatenate([A.{},B.{}[fp]])'.format(varA,varB))
    S.time=concatenate([A.time,B.time[fp]])
    
    fp=~isnan(S.data)
    S.time,S.data=S.time[fp],S.data[fp]
    savez(sname,S)
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

def gplot(x,y,gfactor=10,**kwargs): 
    ''' 
    plot gappy data, process the data first, and then call plot
    '''
    x,y=gappy_data(x,y,gfactor=gfactor)
    plot(x,y,**kwargs)

def gappy_data(x,y,gfactor=10):
    x,y=array(x),array(y)
    if len(x)==0: return array([])
    ind=argsort(x)
    x,y=x[ind],y[ind]
    mdt=median(diff(x))
    tmpx=[]
    for m,xi in enumerate(x):
        if m==0: continue
        dt=xi-x[m-1]
        if dt>mdt*gfactor: #gap larger than multiple times of the median dt
           tmpx.append(x[m-1]+mdt)
    x=concatenate([x,array(tmpx)])
    y=concatenate([y,ones(len(tmpx))*nan])
    ind=argsort(x)
    x,y=x[ind],y[ind]
    return x,y

def download_ndbc(year=2012):
    url='https://www.ndbc.noaa.gov/data/historical/stdmet/'
    cmd=f'wget -r -nc -np -nH -nd -A *{year}.txt.gz {url}; mkdir gz; mv *.txt.gz gz; rm robots.txt.tmp'
    print(cmd); os.system(cmd)
    
    fnames=[i for i in os.listdir('gz') if i.endswith('.txt.gz')]
    for fname in fnames:
        cmd=f'gzip -d gz/{fname}'
        print(cmd); os.system(cmd)
    
def get_skill(x1,y1,x2,y2):
    py1,py2=pair_data(x1,y1,x2,y2,hw=0.2/24)
    return get_stat(py1,py2,fmt=0)

def get_cbibs_data(station='YS',sname='data/cbibs_YS/YS_salt.npz',start_time='2022-01-01T00:00:00',end_time='2032-01-01T00:00:00',
                   just_update=True,isplot=False,var='sea_water_salinity'):
    '''
    For API details, see https://buoybay.noaa.gov/data/api
    '''
    print(f'--- get_cbibs_data for {var} at station {station}',start_time[0:10],end_time[0:10])
    if os.path.exists(sname) and just_update:
        #change the start_time
        S=loadz(sname)
        start_time=num2date(S.time[-1]).strftime("%Y-%m-%dT%H:%M:%S")
        print('change start_time to ',start_time)
        
    apikey='f159959c117f473477edbdf3245cc2a4831ac61f'
    url=f'https://mw.buoybay.noaa.gov/api/v1/json/query/{station}?key={apikey}&sd={start_time}z&ed={end_time}z&var={var}'
    try: 
        urlsave(url,'tmp.json') #reading data
        print('finish get buoy data from remote')
    except:
        print('some error in getting data from',url)
        return
    #process the json data
    print('reading json')
    with open('tmp.json') as f: jdata = json.load(f)
    time,data=[],[]
    for i in jdata['stations'][0]['variable'][0]['measurements']:
        time.append(i['time'])
        data.append(i['value'])
    print('finish reading json')
    if len(time)==0: print('no data found'); return
    time=datenum(time)
    print('finish datenum')
    time,data=array(time),array(data)
    ind=argsort(time)
    time,data=time[ind],data[ind]
    data=data.astype('float')
    if os.path.exists(sname) and just_update:
        fp=time>S.time.max()
        if sum(fp)==0: print('no new records'); return
        S.time=concatenate([S.time,time[fp]]) #append the latest records
        S.data=concatenate([S.data,data[fp]])
        print('append records: ',sum(fp))
    else:
        print('save records')
        S=zdata()
        S.time,S.data=time,data.astype('float')
    if var=='sea_water_salinity':
        S.data[S.data>35]=nan
        S.data[S.data<=1]=nan
    savez(sname,S)
    print(sname)
    if isplot: figure(figsize=[10,3]);  plot(S.time,S.data); set_xtick(fmt=1)

def updated(sname,fnames):
    return os.path.exists(sname) and all([os.path.getmtime(sname)>os.path.getmtime(fname) for fname in fnames])

def get_direction(u,v):
    if v==0:
        if u>0: return 90
        if u<0: return 270
        if u==0: return 0
    else:
        direction=atan(u/v)*180/pi
        if v<0: direction=direction+180
    return direction%360

def find_peaks(x,y,hw=10):
    peaks,peaktimes=[],[]
    for m,yi in enumerate(y):
        if yi>=y[max(0,m-hw):min(m+hw,len(y))].max():
            peaks.append(yi)
            peaktimes.append(x[m])
    return array(peaktimes),array(peaks)

def find_troughs(x,y,hw=10):
    troughs,troughtimes=[],[]
    for m,yi in enumerate(y):
        if yi<=y[max(0,m-hw):min(m+hw,len(y))].min():
            troughs.append(yi)
            troughtimes.append(x[m])
    return array(troughtimes),array(troughs)

def nearest(xi,x):
    dx=abs(xi-x)
    return nonzero(dx==dx.min())[0]

def gap(rel=[1.05,0,0.05,1]): #get subset position based on current one
    #default rel is to add a colorbar
    a=gca().get_position()
    pos=[a.x0,a.y0,a.width,a.height]
    return [pos[0]+pos[2]*rel[0],pos[1]+pos[3]*rel[1],pos[2]*rel[2],pos[3]*rel[3]]

def rtext(x,y,note,xm=None,ym=None,ax=None,**args):
    '''
    add annoation to the current axes at relative location to x-axis and y-axis
       x: relative location in x-axis, tpycially 0-1, but can be negative or larger than 1 
       y: relative location in y-axis, typically 0-1, but can be negative or larger than 1
       note: annotation string
       xm: range of x-axis; gca().xlim() is used if not given
       ym: range of y-axis; gca().ylim() is used if not given
       ax: axes; gca() is used if not given
       **args: all other parameters used in matploblib.plot.text are also applicable

    E.q., rtext(0.1,0.9,"(a)",size=10,weight='bold') -> place label at top left corner

    note: suggest not to rearrange the x/ylim after using this function, or place this after
    '''
    if ax is None: ax=gca()
    sca(ax)
    if xm is None: xm=xlim()
    if ym is None: ym=ylim()
    text(xm[0]+x*diff(xm),ym[0]+y*diff(ym),note,**args)

def extract_sflux(olon= -94.752329,olat=29.351558,sdir='../run17c2/sflux/',sname='narr_galveston.npz'):
    '''
    Extract data from sflux nc files for a specific location; now for air only. 
      olon,olat: location of target
      sdir: sflux directory
      sname: name of npz file to store the data
    '''
    if not sdir.endswith('/'):sdir=sdir+'/'
    nc='sflux_air_1.0001.nc'
    c=ReadNC(sdir+nc)
    dist=sqrt((c.lon.val-olon)**2+(c.lat.val-olat)**2)
    fp=(dist==dist.min()); fp=nonzero(fp)
    ilon,ilat=fp[0][0],fp[1][0]
    #%
    fnames=[i for i in os.listdir(sdir) if 'air' in i]
    fnames.sort()
    u,v,t,h,p,time=[],[],[],[],[],[]
    for m,fname in enumerate(fnames):
        print(sdir+fname)
        c=ReadNC(sdir+fname,fmt=1)
        u.extend(c['uwind'][:,ilon,ilat])
        v.extend(c['vwind'][:,ilon,ilat])
        t.extend(c['stmp'][:,ilon,ilat])
        h.extend(c['spfh'][:,ilon,ilat])
        p.extend(c['prmsl'][:,ilon,ilat])
        time.extend(datenum(2018,1,1)+m+c['time'][:])
    u,v,t,time,p,h=array(u),array(v),array(t),array(time),array(p),array(h)
    C=zdata()
    C.u,C.v,C.t,C.time,C.p,C.h=u,v,t,time,p,h
    print(sname); savez(sname,C)

def make_gif(image_dir='20230730/figs/',begin_string='20230730_biomass',gif='20230730/biomass.gif',duration=100):
    #make gif
    from PIL import Image

    # Sort the PNG files in ascending order
    png_files = sorted([f for f in os.listdir(image_dir) if f.endswith('.png') and f.startswith(begin_string)])

    # Create a list to store the images
    images = []

    # Load each PNG image file and append it to the list
    for png_file in png_files:
        image_path = os.path.join(image_dir, png_file)
        image = Image.open(image_path)
        images.append(image)

    # Save the list of images as a GIF animation
    print('saving GIF to ',gif)
    images[0].save(gif, save_all=True, append_images=images[1:], optimize=True, duration=duration, loop=0)

def get_age(lon,lat,time,rn=30*6): 
    #rn: number of release time for the simulation
    print('Calculating age')
    difflon=diff(lon,axis=0)
    difflat=diff(lat,axis=0)
    fp=(difflon!=0)*(difflat!=0)
    npar=lon.shape[1]
    age=zeros(lon.shape)
    rtime=zeros(npar)
    for ipar in arange(npar):
        itime=nonzero(fp[:,ipar])[0][0] #the first move
        rtime[ipar]=time[itime] #release time
    #use the median value for every npar/rn
    for itime in arange(rn):
        rtime[int(itime*npar/rn):int((itime+1)*npar/rn)]=median(rtime[int(itime*npar/rn):int((itime+1)*npar/rn)])
    for ipar in arange(npar):
        age[:,ipar]=time-rtime[ipar]
    age[age<0]=0
    return age,rtime

def remove_after_hitting_boundary(lon,lat,south=26,east=-87):
    print('Removing track after the particle hit the boundary')
    hit=(lat<=south)|(lon>=east)
    hits=sum(hit,axis=0)
    for ipar in arange(lon.shape[1]):
        if hits[ipar]==0: continue
        itime=nonzero(hit[:,ipar])[0][0] #the first hit
        lon[itime:,ipar]=nan
        lat[itime:,ipar]=nan
    return lon,lat

def remove_sticking(lon,lat):
    print('Removing sticking points')
    difflon=diff(lon,axis=0)
    difflat=diff(lat,axis=0)
    fp=(difflon==0)*(difflat==0)
    lon,lat=lon[1:,:],lat[1:,:]
    lon[fp]=nan; lat[fp]=nan
    return lon,lat

def cal_flushing(nc,reg=None,run=None,iplot=True,rn=30*6,south=26,east=-87,debug=False,min_r2=0.9,begin_age=30,end_age=120):
    print('calculating flushing for',nc)
    if os.path.getsize(nc)<1024**2: print(nc,'file size < 1M');  return
    if reg==None: print('must provide a reg file'); return
    if run==None: print('must mprovide a run name to name npz file'); return
    import numpy as np
    if not fexist('npz/'+run+'_age.npz'):
        C=ReadNC(nc)
        lon=ma.getdata(C.lon.val)
        lat=ma.getdata(C.lat.val)
        mtime=ma.getdata(C.time.val)/86400
        # use nan for lon/lat values when they are masked
        print('removing track when time value is invalid') #occationally, at some writing steps, masked value is used. 
        for itime in nonzero([ma.is_masked(i) for i in C.time.val])[0]: lon[itime,:]=nan; lat[itime,:]=nan  

        #calculate the age, find the release date
        t0=time.time()
        age,rtime=get_age(lon,lat,mtime,rn=rn)
        if age is None:
            print('No particle moving')
            with open('error.out','a') as f:
                f.write(f'\nProblematic run: {run}')
            return

        #remove those hiting southern boundary
        lon,lat=remove_after_hitting_boundary(lon,lat,south=south,east=east)

        print('Cal inpolygon') #take about 20s
        bp=read_schism_reg(reg)
        ind_inpoly=reshape(inside_polygon(c_[lon.ravel(),lat.ravel()],bp.x,bp.y)==1,lon.shape)

        # get release index
        print('Get release index')
        ind_release=tile(np.round((rtime-rtime[0])/mean(diff(unique(rtime)))).astype('int'),[lon.shape[0],1])
        if rn==1: ind_release=tile(0,[lon.shape[0],1])

        # get age index
        print('Get age index')
        timeinterval=median(diff(age[:,0]))
        ind_age=np.round((age/timeinterval)).astype('int')

        # get the number inside the polygon
        print('Get particle number in polygon over time for each release');
        Z=zeros([rn,ind_age.max()+1])
        for i,j,k in zip(ind_release.ravel(),ind_age.ravel(),ind_inpoly.ravel()):
            if i<0 or j<0: continue
            if k: Z[i,j]+=1
        print('Complete calculating in {:.2f}s'.format(time.time()-t0))
        lage=arange(ind_age.max()+1)*timeinterval
        npar=lon.shape[1]/rn
        S=zdata()
        S.npar,S.lage,S.Z,S.rtime=npar,lage,Z,rtime
        savez('npz/'+run+'_age.npz',S)
    else:
        S=loadz('npz/'+run+'_age.npz')
        npar,lage,Z,rtime=S.npar,S.lage,S.Z,S.rtime

    from scipy.stats import linregress
    ftimes=[]
    etimes=[]
    ftimes2=[]
    print('Caculating flushing time')
    Fs,R2s=[],[]
    for ir in arange(rn):
        x=lage[1:]; y=Z[ir,1:]/npar
        x=x[y>0]; y=y[y>0]
        ly=y
        ly=lpfilt(y, 2/24, 0.48) #do low pass filter
        print('   at release No.',ir)
        etime=nan
        if sum(ly<exp(-1))>0: etime=x[nonzero(ly<exp(-1))[0][0]]
        # find efolding time
        ftime=nan
        R2=[]; F=[]
        for age in arange(begin_age,end_age):
            fp=(x<=age)
            m=linregress(x[fp],log(y[fp]))
            slope,rvalue,pvalue=m.slope,m.rvalue,m.pvalue
            R2.append(rvalue**2)
            F.append(-1/slope)
        R2=array(R2); F=array(F)
        Fs.append(F); R2s.append(R2) #save the R2 and F
        if R2[-1]>min_r2: ftime=F[R2==max(R2)][0] #if the linear trend also observable to end the simulation

        # find flushing time based on linear trend 
        ftime2=nan
        fp=(x<=end_age) #use all of them
        m=linregress(x[fp],log(y[fp]))
        slope,rvalue,pvalue=m.slope,m.rvalue,m.pvalue
        if rvalue**2>min_r2: ftime2=-1/m.slope

        etimes.append(etime)
        ftimes.append(ftime)
        ftimes2.append(ftime2)

    S=zdata()
    S.lage=lage; S.release_time=unique(rtime)
    S.Z=Z
    S.npar=npar
    S.ftime=array(ftimes)
    S.ftime2=array(ftimes2)
    S.etime=array(etimes)
    S.begin_age,S.end_age,S.F,S.R2=begin_age,end_age,Fs,R2s
    savez('npz/'+run+'_flushing.npz',S)
    print('Flushing calculation results saved into','npz/'+run+'_flushing.npz')
    print(f'{run} mean etime, ftime, and ftime2: {nanmean(S.etime):.1f},{nanmean(S.ftime):.1f},{nanmean(S.ftime2):.1f}') 
    if iplot: 
        close('all')
        print('Plotting')
        figure(figsize=[6,4])
        for ir in arange(rn):
            plot(lage[1:],Z[ir,1:]/npar*100,lw=0.5,color=[0,0,0,0.2])
        plot([0,90],[exp(-1)*100,exp(-1)*100],'k--')
        plot(lage[1:],Z[:,1:].mean(axis=0)/npar*100,lw=1.5,color=[1,0,0,1])
        xlim([0,90])
        ylabel('Percent of paticles retained (%)')
        xlabel('Days since release')
        rtext(0.2,0.83,'Mean flushing time (from LR)= {:.2f} \u00B1 {:.2f} d'.format(nanmean(ftimes),std(ftimes)))
        rtext(0.2,0.9,'Mean flushing time (from e-folding)= {:.2f} \u00B1 {:.2f} d'.format(nanmean(etimes),std(etimes)))
        title(run)
        tight_layout()
       
        grid('on')
        savefig('figs/'+run+'_flushing.png',dpi=300)

def cal_LET(nc,xlims=[-100,-80],ylims=[18,31],xinterval=0.01,yinterval=0.01,time_interval=2,fname=None,east=-81,south=17.5):
    import numpy as np
    print('Calculating LET'); t0=time.time()
    C=ReadNC(nc)
    lon=ma.getdata(C.lon.val)
    lat=ma.getdata(C.lat.val)
    #% remove sticking points
    lon,lat=remove_sticking(lon,lat)
    #% remove track after hitting boundary
    lon,lat=remove_after_hitting_boundary(lon,lat,east=east,south=south)
    # calculate the index
    indi=np.round((lon-xlims[0])/xinterval).astype('int')
    indj=np.round((lat-ylims[0])/yinterval).astype('int')
    X,Y=meshgrid(arange(xlims[0],xlims[1],xinterval),arange(ylims[0],ylims[1],yinterval))
    Z=zeros(X.shape)
    for i,j,loni,lati in zip(indi.ravel(),indj.ravel(),lon.ravel(),lat.ravel()):
        if i<0 or j<0: continue #they are due to nan value
        Z[j,i]+=1 #be mindful it is [j,i],not [i,j]
    Z=Z/lon.shape[1]*time_interval  #Z/npar*time_interval
    print('Complete calculating LET in {:.2f}s'.format(time.time()-t0))
    if fname!=None:
        S=zdata()
        S.X,S.Y,S.Z=X,Y,Z
        savez(fname,S); print(f'data saved into {fname}')
    return X,Y,Z

def plot_land(gd,color=[0.9,0.9,0.9,1]):
    xv,yv=gd.x[gd.ilbn[0]],gd.y[gd.ilbn[0]]
    xv=[*xv,min(xv)-1,min(xv)-1,max(xv)]
    yv=[*yv,min(yv),max(yv),max(yv)]
    fill(xv,yv,color=color)
    for ia,il in enumerate(gd.ilbn):
        if gd.island[ia]==1:fill(gd.x[il],gd.y[il],color=color)
    xlim([gd.x.min(),gd.x.max()])
    ylim([gd.y.min(),gd.y.max()])

def toDate(otime):
    return [num2date(i) for i in otime]

def make_mp4(mp4,fnames,fps=10):
    fnames.sort()
    images = [imageio.v2.imread(fname) for fname in fnames]
    imageio.mimsave(mp4, images, fps=fps)

def adj_figsize():
    #print('adjut figsize based on lon/lat limits')
    aw,ah=gca().get_position().width,gca().get_position().height
    fw,fh=gcf().get_size_inches()
    fh=(diff(ylim()))*fw*aw/(ah*diff(xlim())*cos(mean(ylim())*pi/180)); fh=fh[0]
    gcf().set_size_inches(fw,fh, forward=True)

def cal_reachtime(nc,sreg,dreg,rname):
    brun=nc.split('/')[-2]
    M=ReadNC(nc)
    dt=mean(diff(M.time.val))
    sregion=sreg.split('/')[-1].split('.')[-2]
    dregion=dreg.split('/')[-1].split('.')[-2]

    lon=M.lon.val.ravel()
    lat=M.lat.val.ravel()

    #taking about 1-2 minutes
    for ireg,reg in enumerate([sreg,dreg]):
        region=reg.split('/')[-1].split('.')[-2]
        sname=f'npz2/{brun}_inreg_{region}.npz'
        reg=read_schism_reg(reg)
        if  not fexist(sname):
            print(f'cal inpolygon {region} for {nc}')
            inreg=inside_polygon(c_[lon,lat],reg.x,reg.y)
            inreg=reshape(inreg,M.lon.val.shape)
            S=zdata(); S.inreg=inreg; savez(sname,S)
        if ireg==0: S=loadz(sname); in_sreg=S.inreg
        if ireg==1: S=loadz(sname); in_dreg=S.inreg

    ntime,npar=M.lon.val.shape
    reachtime=[]
    for ipar in range(npar):
        ins=in_sreg[:,ipar]
        ind=in_dreg[:,ipar]
        if sum(ind==1)==0: #never reach to the destination region
           reachtime.append(nan)
           continue
        t2=nonzero(ind==1)[0][0] #the first time in destination region
        t1=nonzero(ins[:t2]==1)[0][-1]  #last time in sregion, use t2 to avoid reentry case
        reachtime.append((t2-t1)*dt/86400)

    S=zdata()
    S.sregion=sregion; S.dregion=dregion
    S.reachtime=reachtime
    S.ntime,S.npar=ntime,npar
    savez(rname,S)

def running(jname='run19q3_2000'):
    out=sbp.check_output('squeue -u jdu',shell=True).decode().split('\n')
    for rec in out[1:-1]:
        if rec.split()[1]==jname: return True
    return False
