import matplotlib.pyplot as plt
import numpy as np
import os
import datetime
import jdutil
import pdb
import signalsmooth
from taskinit import *
def get_dspec(mspath=None,msfile=None,specfile=None,bl='19&22',timeran=None,chanran=None,savetxt=None):
    # Note: antennas specified in "bl" is antennas INDEX but not antenna NAME.
    ''' test notes '''
    # Open the ms and plot dynamic spectrum
    print 'retrieving spectral data...'
    ms.open(mspath+'/'+msfile)
    ms.selectinit(datadescid=0)
    if bl and (type(bl)==str):
        ms.msselect({'baseline':bl})
    if timeran and (type(timeran)==str):
        ms.msselect({'time':timeran})
    if chanran and (type(chanran)==str):
        [bchan,echan]=chanran.split('~')
        nchan=int(echan)-int(bchan)+1
        ms.selectchannel(nchan,int(bchan),1,1)
    specdata=ms.getdata(['data','time','axis_info'],ifraxis=True)
    npol = specdata['data'].shape[0]
    nchan = specdata['data'].shape[1]
    ntime = specdata['data'].shape[3]
    print 'shape of the spectral data: ', specdata['data'].shape
    print 'number of time pixels: ', specdata['time'].shape[0]
    print 'number of frequency channels selected: ', specdata['axis_info']['freq_axis']['chan_freq'].shape[0]
    spec = specdata['data'].reshape(npol,nchan,ntime)
    tim = specdata['time']
    freq = specdata['axis_info']['freq_axis']['chan_freq'].reshape(nchan)
    #tb.open(mspath+msfile)
    #col='DATA'
    #pol=0 #0 is rcp and 1 is lcp"
    #nspw=1 #number of spectral windows to display
    #sort='casa'
    #subtb=tb.query("ANTENNA1=="+ant1+" && ANTENNA2=="+ant2+" && DATA_DESC_ID>=0 && DATA_DESC_ID<="+str(nspw-1),columns=col)
    #spec=subtb.getcol(col)
    #dim=spec.shape
    #newdim=(dim[0],dim[1]*nspw,dim[2]/nspw)
    #if sort == 'aips':
    #    print 'ms is in aips output order (spw varies faster)'
    #    spec_res=spec.reshape(newdim,order='F')
    #elif sort == 'casa':
    #    print 'ms is in casa sort order (spw varies slower)'
    #    spec_tmp=np.split(spec,nspw,axis=2)
    #    spec_res=np.concatenate(spec_tmp,axis=1)

    # Retrieve time info
    #print 'retrieving time info...'
    #subtb_t=tb.query("ANTENNA1=="+ant1+" && ANTENNA2=="+ant2+" && DATA_DESC_ID==0")
    #tim=subtb_t.getcol("TIME") #time in seconds
    #tb.close()

    # Get the spectral channel information
    #print 'retrieving frequency info...'
    #tb.open(mspath+msfile+'/SPECTRAL_WINDOW')
    #freq=tb.getcol('CHAN_FREQ')
    #tb.close()

    # Save variables
    if not specfile:
        specfile=msfile+'.bl'+bl.replace('&','-')+'.spec.npz'
    np.savez(mspath+'/'+specfile,spec=spec,tim=tim,freq=freq,bl=bl)
    if savetxt:
        spec_r=np.absolute(spec[0,:,:])
        specrtxt=msfile+'.bl'+bl.replace('&','-')+'.spec_r.txt'
        np.savetxt(mspath+'/'+specrtxt,spec_r)
        spec_l=np.absolute(spec[1,:,:])
        specltxt=msfile+'.bl'+bl.replace('&','-')+'.spec_l.txt'
        np.savetxt(mspath+'/'+specltxt,spec_l)
        timtxt=msfile+'.bl'+bl.replace('&','-')+'.tim.txt'
        timstrs=[]
        for t in tim:
            timstrs.append(qa.time(qa.quantity(t,'s'),prec=9)[0])
        np.savetxt(mspath+'/'+timtxt,timstrs,fmt='%s')
        freqtxt=msfile+'.bl'+bl.replace('&','-')+'.freq.txt'
        np.savetxt(mspath+'/'+freqtxt,freq)

def get_med_dspec(mspath=None,msfile=None,specfile=None,\
                  smoothmode='base',tsmoothwidth=40,tsmoothtype='hanning',\
                  fsmoothwidth=4,fsmoothtype='hannning',\
                  fsmooth=1,timeran=None,chanran=None,dmin=0.5,dmax=2.0):
    # Open the ms and plot dynamic spectrum
    print 'retrieving spectral data...'
    ms.open(mspath+'/'+msfile)
    ms.selectinit(datadescid=0)
    if timeran and (type(timeran)==str):
        ms.msselect({'time':timeran})
    if chanran and (type(chanran)==str):
        [bchan,echan]=chanran.split('~')
        nchan=int(echan)-int(bchan)+1
        ms.selectchannel(nchan,int(bchan),1,1)
    specdata=ms.getdata(['data','time','axis_info'],ifraxis=True)
    npol = specdata['data'].shape[0]
    nchan = specdata['data'].shape[1]
    nbl = specdata['data'].shape[2]
    ntime = specdata['data'].shape[3]
    print 'shape of the spectral data: ', specdata['data'].shape
    print 'number of time pixels: ', specdata['time'].shape[0]
    print 'number of frequency channels selected: ', specdata['axis_info']['freq_axis']['chan_freq'].shape[0]
    print 'number of baselines: ', specdata['axis_info']['ifr_axis']['baseline'].shape[0]
    # ratio over time-average
    if smoothmode=='hipass':
        specsmooth=np.copy(specdata['data'])
        # time smooth
        for i in range(npol):
            for j in range(nchan):
                for k in range(nbl):
                    specsmooth[i,j,k,:] = signalsmooth.smooth(specdata['data'][i,j,k,:],tsmoothwidth,tsmoothtype)
        # frequency smooth
        for i in range(npol):
            for j in range(nbl):
                for k in range(ntime):
                    specsmooth[i,:,j,k] = signalsmooth.smooth(specsmooth[i,:,j,k],fsmoothwidth,fsmoothtype)
        specratio=specdata['data']/specsmooth
    if smoothmode=='base':
        specsmooth=np.mean(specdata['data'],axis=3)
        specratio=specdata['data']/specsmooth[:,:,:,None]
        
    spec_med=np.median(specratio,axis=2)
    if fsmooth:
        for i in range(npol):
            for j in range(ntime):
                    spec_med[i,:,j] = signalsmooth.smooth(spec_med[i,:,j],fsmoothwidth,fsmoothtype)

    tim = specdata['time']
    freq = specdata['axis_info']['freq_axis']['chan_freq'].reshape(nchan)
    tim0 = tim[0]
    tim0str=qa.time(qa.quantity(tim0,'s'),prec=8)[0]
    tim_ = tim-tim[0]
    freqghz = freq/1e9
    f=plt.figure(figsize=(8,8),dpi=100)
    ax1=f.add_subplot(211)
    f.subplots_adjust(hspace=0.4)
    ax1.pcolormesh(tim_,freqghz,np.abs(spec_med[0,:,:]),cmap='jet',vmin=dmin,vmax=dmax)
    ax1.set_xlim([tim_[0],tim_[-1]])
    ax1.set_ylim([freqghz[0],freqghz[-1]])
    ax1.set_title('Median-Filtered Dynamic Spectrum (RCP)')
    ax1.set_xlabel('Time (seconds) since '+tim0str)
    ax1.set_ylabel('Frequency (GHz)')
    ax2=f.add_subplot(212)
    ax2.pcolormesh(tim_,freqghz,np.abs(spec_med[1,:,:]),cmap='jet',vmin=dmin,vmax=dmax)
    ax2.set_xlim([tim_[0],tim_[-1]])
    ax2.set_ylim([freqghz[0],freqghz[-1]])
    ax2.set_title('Median-Filtered Dynamic Spectrum (LCP)')
    ax2.set_xlabel('Time (seconds) since '+tim0str)
    ax2.set_ylabel('Frequency (GHz)')
    if specfile:
        np.savez(mspath+'/'+specfile,spec_med=spec_med,tim=tim,freq=freq)
    


def plt_dspec(mspath=None,specfile=None,pol='I', dmin=None, dmax=None, fig=None):
    # Set up variables 
    if pol != 'RR' and pol != 'LL' and pol !='RRLL' and pol != 'I' and pol != 'V' and pol != 'IV':
        print "Please enter 'RR', 'LL', 'RRLL', 'I', 'V', 'IV' for pol"
        return 0

    # import the spec data, time, and frequency information
    #for files in os.listdir(mspath):
    #    if files.endswith('.spec.npz'):
    #        infile=files

    specdata=np.load(mspath+'/'+specfile)
    spec=specdata['spec']
    tim=specdata['tim']
    freq=specdata['freq']

    # setup plot parameters
    print 'ploting dynamic spectrum...'
    # mask the channels from 512 to 519 (no observation)
    spec=np.ma.array(spec)
    spec[:,512:519,:]=np.ma.masked
    spec_med=np.median(np.absolute(spec))
    # set the time axis
    ntim=tim.shape
    if ntim[0] < 20:
        xticks=np.arange((ntim[0]-1)/2+1)*2
    elif ntim[0] < 100:
        xticks=np.arange((ntim[0]-1)/20+1)*20
    elif ntim[0] < 600:
        xticks=np.arange((ntim[0]-1)/100+1)*100
    elif ntim[0] < 2400:
        xticks=np.arange((ntim[0]-1)/400+1)*400
    elif ntim[0] < 12000:
        # 1 min per step
        tstart=np.fix(tim[0]/60.)+1
        xstart=np.abs(tim-tstart*60.).argmin()
        xticks=np.arange(ntim[0]/1200)*1200+xstart
    elif ntim[0] > 12000:
        xticks=np.arange(ntim[0]/6000+1)*6000
    nticks=xticks.shape
    xticktims=[]
    for i in range(nticks[0]):
        #xticktim0=qa.time(qa.quantity(tim[xticks[i]],'s'))
        tims=tim[xticks[i]] #in seconds
        tims_jd=jdutil.mjd_to_jd(tims/3600./24.) #to julian date
        tims_dt=jdutil.jd_to_datetime(tims_jd)
        tims_dt2=tims_dt+datetime.timedelta(seconds=round(tims_dt.microsecond/1e6))
        tims_char=tims_dt2.strftime('%H:%M:%S')
        xticktims.append(tims_char)
    xticks=list(xticks)
    # do the plot
    f=plt.figure(figsize=(10,6),dpi=100)
    if not dmin:
        dmin=spec_med/20.
    if not dmax:
        dmax=spec_med*5.
    if pol != 'RRLL' and pol != 'IV':
        ax=f.add_subplot(111)
        if pol=='RR':
            spec_plt=np.absolute(spec[0,:,:])
        elif pol=='LL':
            spec_plt=np.absolute(spec[1,:,:])
        elif pol=='I':
            spec_plt=(np.absolute(spec[0,:,:])+np.absolute(spec[1,:,:]))/2.
        elif pol=='V':
            spec_plt=(np.absolute(spec[0,:,:])-np.absolute(spec[1,:,:]))/2.
        ax.imshow(spec_plt,aspect='auto',origin='lower',extent=[0,ntim[0]-1,np.min(freq)/1e6,np.max(freq)/1e6],\
                  vmin=dmin,vmax=dmax,interpolation='none')
        ax.set_xticks(xticks)
        ax.set_xticklabels(xticktims)
        ax.set_xlabel('Universal Time (50 ms/pixel)')
        ax.set_ylabel('Frequency (MHz)')
        ax.set_title('VLA dynm spec for pol '+pol)
        ax.set_autoscale_on(False)
    else:
        R_plot=np.absolute(spec[0,:,:])
        L_plot=np.absolute(spec[1,:,:])
        I_plot=(Rplot+Lplot)/2.
        V_plot=(Rplot-Lplot)/2.
        if pol == 'RRLL':
            spec_plt_1=R_plot
            spec_plt_2=L_plot
        if pol == 'IV':
            spec_plt_1=I_plot
            spec_plt_2=V_plot

        ax1=f.add_subplot(211)
        ax1.imshow(spec_plt_1,aspect='auto',origin='lower',extent=[0,ntim[0]-1,np.min(freq)/1e6,np.max(freq)/1e6],\
                   vmin=dmin,vmax=dmax,interpolation='none')
        ax1.set_xticks(xticks)
        ax1.set_xticklabels(xticktims)
        ax1.set_xlabel('Universal Time (50 ms/pixel)')
        ax1.set_ylabel('Frequency (MHz)')
        ax2.set_title('VLA dynm spec for pol '+pol[0])
        ax1.set_autoscale_on(False)
        ax2=f.add_subplot(212)
        ax2.imshow(spec_plt_2,aspect='auto',origin='lower',extent=[0,ntim[0]-1,np.min(freq)/1e6,np.max(freq)/1e6],\
                   vmin=dmin,vmax=dmax,interpolation='none')
        ax2.set_xticks(xticks)
        ax2.set_xticklabels(xticktims)
        ax2.set_xlabel('Universal Time (50 ms/pixel)')
        ax2.set_ylabel('Frequency (MHz)')
        ax2.set_title('VLA dynm spec for pol '+pol[-1])
        ax2.set_autoscale_on(False)
    if not fig:
        figfile=specfile[:(specfile.find('spec'))]+pol+'.pdf'
    else:
        figfile=fig
    f.savefig(mspath+'/'+figfile)
    #plt.close()
    return (f,ax)
