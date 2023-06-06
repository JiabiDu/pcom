#!/usr/bin/env python3
from pylib import *
from pcom import *
#%%
def plot_salinity(stations,M,B,O,row,col,margin=[0.05,0.02,0.06,0.05],dxy=[0.03,0.00],p=None):
    run,brun,gd,reftime=p.run,p.brun,p.gd,p.reftime
    istack=0
    ps=get_subplot_position2(margin=margin,dxy=dxy,ds=[row,col]);ps=reshape(ps,(row*col,4))
    figure(figsize=[6*col,2*row])
    xlims=[M.time.min()+reftime-2,M.time.max()*1.2+reftime]
    m=-1
    for station in stations:   
        fp=(O.station==station)*(O.var=='seawater_salinity')*(O.time>M.time.min()+reftime)*(O.time<M.time.max()+reftime)
        if sum(fp)==0 and station!=stations[-1]: continue
        m+=1
        #plot the model surface and bottom
        print('plot salt for station ',m+1,station)
        if m%(col*row)==0:  clf()
        axes(position=ps[m%(col*row)])
        maxs,mins=0,35 #for ylim purpose
        
        olon,olat=O.lon[station],O.lat[station]
        fp=(O.station==station)*(O.var=='seawater_salinity')*(O.time>=reftime)*(O.time<=M.time.max()+reftime)
        #plot(O.time[fp],O.data[fp],'.',color=[0.5,0.5,0.5,1],lw=1,ms=1)
        x,y=O.time[fp],O.data[fp]
        fp=~isnan(y)
        if sum(fp)>0:
            x,y=x[fp],y[fp]
            y=lpfilt(y, median(diff(x)), 0.48); y[-96:]=nan #get subtidal
            plot(x,y,'.',color=[0.5,0.5,0.5,1],lw=1,ms=1)
            maxs=max(maxs,max(y)); mins=min(mins,min(y))
            ox,oy=x,y
        
        #for bottom and surface in base run; only show the subtidal
        if B!=None and sum(B.bp.station==station)>0:
            fp=(B.bp.station==station)*(B.bp.z==0); #surface
            fp=nonzero(fp)[0][0] 
            y=lpfilt(B.salinity[fp].ravel(), median(diff(B.time)), 0.48); y[-96:]=nan #get subtidal
            plot(B.time+reftime,y,'-',color=[0.3,0.3,1,1],lw=0.5)
            maxs=max(maxs,max(y)); mins=min(mins,min(y))
            Sb=get_skill(ox,oy,B.time+reftime,y)
            
            fp=(B.bp.station==station)*(B.bp.z==max(B.bp.z[B.bp.station==station])); #bottom
            fp=nonzero(fp)[0][0] 
            y=lpfilt(B.salinity[fp].ravel(), median(diff(B.time)), 0.48); y[-96:]=nan #get subtidal
            plot(B.time+reftime,y,'-',color=[0,0,0.7,1],lw=0.5)
            maxs=max(maxs,max(y)); mins=min(mins,min(y))

        # for surface
        fp=(M.bp.station==station)*(M.bp.z==0); fp=nonzero(fp)[0][0] #surface
        #plot(M.time+reftime,M.salinity[fp].ravel(),'-',color=[0,0.0,0.8,0.5],lw=0.2)
        # get the subtidal
        y=lpfilt(M.salinity[fp].ravel(), median(diff(M.time)), 0.48); y[-96:]=nan
        plot(M.time+reftime,y,'-',color=[1,0.3,0.3,1],lw=1)
        maxs=max(maxs,max(y)); mins=min(mins,min(y))
        Sm=get_skill(ox,oy,M.time+reftime,y)
        
        # for bottom
        fp=(M.bp.station==station)*(M.bp.z==max(M.bp.z[M.bp.station==station])); 
        fp=nonzero(fp)[0][0] #bottom
        #plot(M.time+reftime,M.salinity[fp].ravel(),'-',color=[0.8,0.0,0.0,0.5],lw=0.5)
        # get the subtidal
        y=lpfilt(M.salinity[fp].ravel(), median(diff(M.time)), 0.48); y[-96:]=nan
        plot(M.time+reftime,y,'-',color=[0.7,0.0,0,1],lw=1)
        maxs=max(maxs,max(y)); mins=min(mins,min(y))
        
        xlim(xlims) #times 1.2 to allow some space showing the monitoring station in a subset
        ylim([floor(mins),ceil(maxs)])
        set_xtick(fmt=0)
        rtext(0.02,0.9,'MAE: ')
        rtext(0.13,0.9,'{:.2f}'.format(Sm.MAE),color='r')
        rtext(0.2,0.9,'{:.2f}'.format(Sb.MAE),color='b')
        
        if m%(row*col)<(row-1)*col and m!=len(stations)-1: setp(gca(),xticklabels=[]) #only the last row will have time stamp
        #if m%(row*col)==0: legend([run+'-S',run+'-B',brun+'-S',brun+'-B','OBS'],ncol=3)
        if m%(row*col)==1 and B!=None: title(f'{run}(red); {brun}(blue); observation(gray)')
        if m%(row*col)==1 and B==None: title(f'{run}(red); observation(gray)',fontsize=12,fontweight='bold')
        
        # show the monitoring station
        pa=gca().get_position()
        pa=[pa.x0,pa.y0,pa.width,pa.height]
        pa=[pa[0]+0.84*pa[2],pa[1]+pa[3]*0.0,pa[2]*0.16,pa[3]*1.0]
        axes(position=pa) 
        fill([olon-0.2,olon+0.2,olon+0.2,olon-0.2],[olat-0.4,olat-0.4,olat+0.4,olat+0.4],color=[0.6,0.6,0.6,1])
        gd.plot_bnd(lw=1)
        #gO.plot(fmt=1,cmap='Blues',clim=[0,15],cb=False)
        gd.plot_bnd(color='c')
        plot(olon,olat,'mx',ms=5,markerfacecolor='m')
        xlim([olon-0.2,olon+0.2]); ylim([olat-0.4,olat+0.4])
        rtext(0.1,0.9,station,fontweight='bold')
        
        setp(gca(),xticks=[],yticks=[])
        axis('off')
        gca().set_facecolor((1.0, 0.47, 0.42))
        if m%(row*col)==(row*col)-1 or station==stations[-1]:
            istack+=1
            #tight_layout()
            figname=run+'/figs/salt_{:02d}_{}_vs_{}.png'.format(istack,run,brun)
            print(figname); savefig(figname,dpi=300)

def get_skill(x1,y1,x2,y2):
    py1,py2=pair_data(x1,y1,x2,y2,hw=0.2/24)
    return get_stat(py1,py2,fmt=0)

def plot_salinity_one_layer(M,B,O,fmt=1,stations=['BOLI','MIDG','FISH','TRIN'],deps=[3,3,1.5,1.5],bay='Galveston',p=None):
    lw=0.5; 
    if fmt==1: lw=1
    run,brun,gd,reftime=p.run,p.brun,p.gd,p.reftime
    #fmt=0 and 1 for full signal and subtidal signal
    figure(figsize=[6,2*len(stations)]); offset=0;
    for m,station in enumerate(stations):
        fp=(O.station==station)*(O.var=='seawater_salinity')*(O.time>M.time.min()+reftime)*(O.time<M.time.max()+reftime)
        if sum(fp)==0: continue
        print('plot salt for station ',station) 
        subplot(len(stations),1,m+1)
        ox,oy=O.time[fp],O.data[fp]
        if fmt==1: oy=lpfilt(O.data[fp],median(diff(O.time[fp])),0.48)
        plot(ox,oy,'x',color=[0,0,0,0.5],lw=lw,ms=1)
    
        if  brun!=None: 
            fp=B.bp.station==station
            fp0=nonzero(fp)[0]
            tmpz=B.bp.z[fp0]
            tmp=abs(tmpz-deps[m])
            tfp=fp0[tmp==tmp.min()]
            y=B.salinity[tfp].ravel()
            if fmt==1: y=lpfilt(B.salinity[tfp].ravel(),median(diff(B.time)),0.48) #0.48==50hour
            plot(B.time+reftime,y,'-',color=[0,0.8,0.8,0.8],lw=lw)
            mxB=B.time+reftime
            myB=lpfilt(B.salinity[tfp].ravel(),median(diff(M.time)),0.48) #0.48==50hour
            
        fp=M.bp.station==station
        fp0=nonzero(fp)[0]
        tmpz=M.bp.z[fp0]
        tmp=abs(tmpz-deps[m]) #find the depth closest to the needed depth. 
        tfp=fp0[tmp==tmp.min()]
        mx=M.time+reftime
        my=M.salinity[tfp].ravel()
        if fmt==1: my=lpfilt(M.salinity[tfp].ravel(),median(diff(M.time)),0.48) #0.48==50hour
        plot(mx,my,'-',color=[1,0,0,0.7],lw=lw)    
    
        xlim([M.time.min()+reftime-2,M.time.max()+reftime])
        ylim([0,36])
        set_xtick(fmt=2)
        #bs_diff=M.salinity[bfp].ravel().mean()-M.salinity[sfp].ravel().mean()
        S=get_skill(ox,oy,mx,my)
        if brun!=None: SB=get_skill(ox,oy,mxB,myB)
        r"$\bf{" + str(number) + "}$"
        if brun!=None: 
            rtext(0.01,0.90,'('+chr(97+m)+') '+ r"$\bf{"+station+"}$"+'; MAE={:.2f} vs {:.2f} in {}'.format(S.MAE,SB.MAE,brun),backgroundcolor=[1,1,1,0.5])
        else:
            rtext(0.01,0.90,'('+chr(97+m)+') '+ r"$\bf{"+station+"}$"+'; MAE={:.2f}'.format(S.MAE),backgroundcolor=[1,1,1,0.5])
        grid('on',ls=':')
        if m==0: 
            title(run)
            if brun!= None: legend(['obs',brun,run])
            if brun==None: legend(['obs',run])
        if m==0: ylabel('Salinity (psu)')
    tight_layout()
    tmp=''
    if fmt==1: tmp='_sub_tidal'
    figname=f'{run}/figs/salt_{bay}{tmp}_{run}_vs_{brun}.png'
    print(figname); savefig(figname,dpi=300)
    
#%%
def plot_elevation(stations,M,B,O,row,col,margin=[0.05,0.02,0.06,0.1],dxy=[0.05,0.00],fmt=0,p=None,xlims=None,ylims=None):
    #fmt=0: full signal
    #fmt=1: subtidal
    run,brun,gd,reftime,station_name=p.run,p.brun,p.gd,p.reftime,p.station_name
    istack=0
    ps=get_subplot_position2(margin=margin,dxy=dxy,ds=[row,col]);ps=reshape(ps,(row*col,4))
    figure(figsize=[6*col,2*row])
    if xlims==None: xlims=[M.time.min()+reftime-2,M.time.max()*1.2+reftime]
    m=-1
    for station in stations:
        m+=1
        #plot the model surface and bottom
        print('plot elevation for station ',m+1,station)
        if m%(col*row)==0:  clf()
        axes(position=ps[m%(col*row)])
        if B!=None: #old model run
            fp=B.bp.station==station
            if fmt==0: plot(B.time+reftime,B.elevation[fp].ravel(),'-',color=[0,0.8,0.8,0.6],lw=0.5)
            foyi=lpfilt(B.elevation[fp].ravel(),diff(B.time).mean(),0.2)
            if fmt==1: plot(B.time+reftime,foyi,'-',color=[0.0,0.8,0.8,1],lw=1);btime,bfoyi=B.time+reftime,foyi
            if fmt==2: plot(B.time+reftime,B.elevation[fp].ravel()-foyi,'-',color=[0.0,0.8,0.8,1],lw=0.5) #tidal signal
       
        #new model run
        fp=M.bp.station==station
        if fmt==0: plot(M.time+reftime,M.elevation[fp].ravel(),'-',color=[1,0,0,0.6],lw=0.5) #full signal
        foyi=lpfilt(M.elevation[fp].ravel(),diff(M.time).mean(),0.2)
        meanModel=mean(foyi)
        if fmt==1: plot(M.time+reftime,foyi,'-',color=[0.8,0,0,1],lw=1);mtime,mfoyi=M.time+reftime,foyi #subtidal
        if fmt==2: plot(M.time+reftime,M.elevation[fp].ravel()-foyi,'-',color=[0.8,0,0,1],lw=0.5) #tidal signal
        
        #observation
        fp=(O.station==station)*(O.time>M.time.min()+reftime)*(O.time<min(M.time.max()+reftime,xlims[1]))*(~isnan(O.wl))
        if sum(fp)>0:
            foyi=lpfilt(O.wl[fp],1/24,0.2)
            offset=meanModel-nanmean(foyi) #to correct the datum; TODO: use the median diff in subtidal signal
            if isnan(offset): offset=0
            if fmt==0: plot(O.time[fp],O.wl[fp]+offset,'-',color=[0,0,0,0.6],lw=0.5)
            if fmt==1: 
                plot(O.time[fp],foyi+offset,'-',color=[0,0,0,1],lw=1)
                Sm=get_skill(O.time[fp],foyi+offset,mtime,mfoyi)
                Sb=get_skill(O.time[fp],foyi+offset,btime,bfoyi)
            if fmt==2: 
                plot(O.time[fp],O.wl[fp]-foyi,'-',color=[0,0,0,0.6],lw=0.5); 
                if ylims==None:ylims=[min(O.wl[fp]-foyi)-0.2,max(O.wl[fp]-foyi)+0.2]
                ylim(ylims)
            print(offset)
        set_xtick(fmt=2)
        xlim(xlims)
        rtext(0.01,0.9,station+':'+station_name[station]+' offset={:.2f}'.format(offset))
        if fmt==1 and sum(fp)>0: rtext(0.01,0.7,'MAE: {:.2f}m vs {:.2f}m'.format(Sm.MAE,Sb.MAE))
        
        if m%(row*col)<col*(row-1) and station != station[-1]: setp(gca(),xticklabels='')
        
        # show the monitoring station
        pa=gca().get_position()
        pa=[pa.x0,pa.y0,pa.width,pa.height]
        pa=[pa[0]+0.84*pa[2],pa[1]+pa[3]*0.0,pa[2]*0.16,pa[3]*1.0]
        axes(position=pa) 
        gd.plot_bnd(lw=0.5)
        #gO.plot(fmt=1,cmap='Blues',clim=[0,15],cb=False)
        gd.plot_bnd(color=[0.4,0.4,0.4,0.8],lw=1)
        olon,olat=p.lon[station],p.lat[station]
        plot(olon,olat,'mo',ms=5,markerfacecolor=[1,0,1,1])
        xlim([olon-1,olon+1]); ylim([olat-1,olat+1])
        setp(gca(),xticks=[],yticks=[])
        axis('off')
        
        if m%(row*col)==(row*col)-1 or station==stations[-1]:
            istack+=1
            #tight_layout()
            if fmt==0: figname=run+'/figs/elev_{:02d}_{}_vs_{}.png'.format(istack,run,brun)
            if fmt==1: figname=run+'/figs/elev_subtidal_{:02d}_{}_vs_{}.png'.format(istack,run,brun)
            if fmt==2: figname=run+'/figs/elev_tidal_{:02d}_{}_vs_{}.png'.format(istack,run,brun)
            print(figname); savefig(figname,dpi=300)

#%% for current 
def get_direction(u,v):
    if v==0:
        if u>0: return 90
        if u<0: return 270
        if u==0: return 0
    else:
        direction=atan(u/v)*180/pi
        if v<0: direction=direction+180
    return direction%360

#% show the current station and name
def plot_obs_station(M,gd=None,var='current',p=None):
    run=p.run
    xlims=[M.bp.x.min()-0.2,M.bp.x.max()+0.2]
    ylims=[M.bp.y.min()-0.2,M.bp.y.max()+0.2]
    figure(figsize=[16,*8*diff(ylims)*2/diff(xlims)])
    if gd!=None: gd.plot(fmt=1,clim=[0,20],cb=False)
    plot(M.bp.x,M.bp.y,'ro')
    for lon,lat,station in zip(M.bp.x,M.bp.y,M.bp.station):
        text(lon,lat,station,fontsize=6)
    xlim(xlims)
    ylim(ylims)
    tight_layout()
    savefig(f'{run}/figs/{var}_obs_locations.png',dpi=300)

def plot_current_direction_hist(M,O,p=None):
    run=p.run
    figure(figsize=[20,9])
    for m,station in enumerate(M.bp.station): #M.bp.station:  #g06010 entrance, bin1 at the surface
        # plot the observational data
        ofp=O.station==station
        ubin=unique(O.bin[ofp])
        ofp=(O.station==station)*(O.bin==1)
        ox,oy=gappy_data(O.time[ofp],O.speed[ofp])
        subplot(5,7,m+1)
        hist(O.direction[ofp],bins=arange(5,360,10),color='k')
        mfp=M.bp.station==station
        direction=array([get_direction(u,v) for u,v in zip(M.horizontalVelX[mfp].ravel(),M.horizontalVelY[mfp].ravel())])
        hist(direction,bins=arange(5,360,10),color=[1,0.3,0.3,0.8])
        rtext(0.1,0.9,station+': {:.0f} vs {:.0f}'.format(median(O.direction[ofp]),median(direction)))
        #axis('off')
    tight_layout()
    savefig(f'{run}/figs/current_direction_hist.png',dpi=300)

def plot_current_speed_hist(M,O,p=None):
    run=p.run
    figure(figsize=[20,9])
    for m,station in enumerate(M.bp.station): #M.bp.station:  #g06010 entrance, bin1 at the surface
        # plot the observational data
        ofp=O.station==station
        ubin=unique(O.bin[ofp])
        ofp=(O.station==station)*(O.bin==1)
        subplot(5,7,m+1)
        hist(O.speed[ofp].ravel(),bins=arange(0,100,5),color='k')
        mfp=M.bp.station==station
        hist(M.speed[mfp].ravel()*100,bins=arange(0,100,5),color=[1,0.3,0.3,0.8])
        rtext(0.1,0.9,station+': {:.0f} vs {:.0f}'.format(median(O.speed[ofp]),median(M.speed[mfp]*100)))
        #axis('off')
    tight_layout()
    savefig(f'{run}/figs/current_speed_hist.png',dpi=300)

def plot_current(M,O,B,p=None,zoom=0):
    run,brun,reftime,gd=p.run,p.brun,p.reftime,p.gd
    row,col=5,1
    margin=[0.05,0.02,0.06,0.1]; dxy=[0.03,0.05]
    ps=get_subplot_position2(margin=margin,dxy=dxy,ds=[row,col]);ps=reshape(ps,(row*col,4))
    figure(figsize=[12*col,2*row])
    istack=0;m=-1
    stations=M.bp.station
    for station in stations: #M.bp.station:  #g06010 entrance, bin1 at the surface
        m+=1
        if m%(col*row)==0:  clf()
        axes(position=ps[m%(col*row)])
        
        # plot the observational data
        ofp=O.station==station
        ubin=unique(O.bin[ofp])
        ofp=(O.station==station)*(O.bin==1)
        ox,oy=gappy_data(O.time[ofp],O.speed[ofp])
        if zoom!=0:  #zoom in the first few days (defined by zoom)
            fp=ox<=ox.min()+zoom
            ox,oy=ox[fp],oy[fp]
            
        #find the main direction, use histogram, or the medi
        plot(ox-reftime,oy,'-',color=[0,0,0,1],lw=0.5)
        if B!=None:
            mfp=B.bp.station==station
            plot(B.time,B.speed[mfp].ravel()*100,'-',color=[0.3,0.3,1,0.5],lw=0.5)
        mfp=M.bp.station==station
        tfp=M.time+reftime<=ox.max()
        plot(M.time[tfp],M.speed[mfp,tfp].ravel()*100,'-',color=[1,0.3,0.3,1],lw=0.5)
        xlim([ox.min()-reftime,ox.min()+(ox.max()-ox.min())*1.2-reftime])
        #set_xtick(fmt=1)
        rtext(0.05,0.9,station,backgroundcolor='w',fontweight='bold')
        
        #if m%(row*col)<(row-1)*col and m!=len(stations)-1: setp(gca(),xticklabels=[]) #only the last row will have time stamp
        if m%(row*col)==0: title(f'{run}(red); {brun}(blue)')
        if m%(row*col)==(row-1)*col or station==M.bp.station[-1]: xlabel('Days since {}'.format(str(num2date(p.reftime))[0:10]))
        
        #-- add station location
        # show the monitoring station
        olon,olat=M.bp.x[mfp],M.bp.y[mfp]
        pa=gca().get_position()
        pa=[pa.x0,pa.y0,pa.width,pa.height]
        pa=[pa[0]+0.84*pa[2],pa[1]+pa[3]*0.0,pa[2]*0.16,pa[3]*1.0]
        axes(position=pa) 
        gd.plot_bnd(lw=0.5)
        #gO.plot(fmt=1,cmap='Blues',clim=[0,15],cb=False)
        gd.plot_bnd(color=[0,0.8,0.8,0.8])
        plot(olon,olat,'mo',ms=5,markerfacecolor=[1,0,1,0.5])
        xlim([olon-0.5,olon+0.5]); ylim([olat-0.5,olat+0.5])
        setp(gca(),xticks=[],yticks=[])
        axis('off')
        
        
        if m%(row*col)==(row*col)-1 or station==stations[-1]:
            istack+=1
            #tight_layout()
            xlabel('Days since {}'.format(str(num2date(p.reftime))[0:10]))
            figname=run+f'/figs/current_staack{istack}_zoom{zoom}_{run}_vs_{brun}.png'
            print(figname); savefig(figname)
    
#%%
def plot_temp(stations,M,B,O,row,col,margin=[0.05,0.02,0.06,0.05],dxy=[0.03,0.00],p=None):
    run,brun,gd,reftime=p.run,p.brun,p.gd,p.reftime
    istack=0
    ps=get_subplot_position2(margin=margin,dxy=dxy,ds=[row,col]);ps=reshape(ps,(row*col,4))
    figure(figsize=[6*col,2*row])
    xlims=[M.time.min()+reftime-2,M.time.max()*1.2+reftime]
    m=-1
    for station in stations:   
        fp=(O.station==station)*(O.time>=reftime)*(O.time<=M.time.max()+reftime)
        if sum(fp)==0 and station!=stations[-1]: continue #no observation
        m+=1
        print('plot temperature for station ',m+1,station)
        if m%(col*row)==0:  clf()
        axes(position=ps[m%(col*row)])
        maxs,mins=0,35 #for ylim purpose
        olon,olat=M.bp.x[M.bp.station==station][0],M.bp.y[M.bp.station==station][0]
        
        x,y=O.time[fp],O.temp[fp]
        fp=~isnan(y)
        if sum(fp)>0:
            x,y=x[fp],y[fp]
            y=lpfilt(y, median(diff(x)), 0.48); y[-96:]=nan #get subtidal
            plot(x,y,'.',color=[0.5,0.5,0.5,1],lw=1,ms=1)
            maxs=max(maxs,max(y)); mins=min(mins,min(y))
        
        if B!=None and sum(B.bp.station==station)>0:
            fp=(B.bp.station==station)*(B.bp.z==0); #surface
            fp=nonzero(fp)[0][0] 
            y=lpfilt(B.temperature[fp].ravel(), median(diff(B.time)), 0.48); y[-96:]=nan #get subtidal
            plot(B.time+reftime,y,'-',color=[0.3,0.3,1,1],lw=0.5)
            maxs=max(maxs,max(y)); mins=min(mins,min(y))

        # for surface
        fp=(M.bp.station==station)*(M.bp.z==0); fp=nonzero(fp)[0][0] #surface
        #plot(M.time+reftime,M.salinity[fp].ravel(),'-',color=[0,0.0,0.8,0.5],lw=0.2)
        # get the subtidal
        y=lpfilt(M.temperature[fp].ravel(), median(diff(M.time)), 0.48); y[-96:]=nan
        plot(M.time+reftime,y,'-',color=[1,0.3,0.3,1],lw=1)
        maxs=max(maxs,max(y)); mins=min(mins,min(y))

        xlim(xlims) #times 1.2 to allow some space showing the monitoring station in a subset
        ylim([floor(mins),ceil(maxs)])
        set_xtick(fmt=0)
        if maxs<20:
            rtext(0.05,0.9,station+'({:.2f},{:.2f})'.format(olon,olat),backgroundcolor='w')
        else:
            rtext(0.05,0.15,station+'({:.2f},{:.2f})'.format(olon,olat),backgroundcolor='w')
        if m%(row*col)<(row-1)*col and m!=len(stations)-1: setp(gca(),xticklabels=[]) #only the last row will have time stamp
        if m%(row*col)==0: ylabel('Water temp (C)')
        if m%(row*col)==1 and B!=None: title(f'{run}(red); {brun}(blue); observation(gray)')
        if m%(row*col)==1 and B==None: title(f'{run}(red); observation(gray)',fontsize=12,fontweight='bold')
        
        # show the monitoring station
        pa=gca().get_position()
        pa=[pa.x0,pa.y0,pa.width,pa.height]
        pa=[pa[0]+0.84*pa[2],pa[1]+pa[3]*0.0,pa[2]*0.16,pa[3]*1.0]
        axes(position=pa) 
        gd.plot_bnd(lw=0.5)
        #gO.plot(fmt=1,cmap='Blues',clim=[0,15],cb=False)
        gd.plot_bnd(color='c')
        plot(olon,olat,'mx',ms=3,markerfacecolor='m')
        xlim([olon-0.2,olon+0.2]); ylim([olat-0.4,olat+0.4])
        setp(gca(),xticks=[],yticks=[])
        axis('off')
        if m%(row*col)==(row*col)-1 or station==stations[-1]:
            istack+=1
            #tight_layout()
            figname=run+'/figs/temp_{:02d}_{}_vs_{}.png'.format(istack,run,brun)
            print(figname); savefig(figname,dpi=300)