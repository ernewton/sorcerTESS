from astropy.io import fits
import NoiseModel
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d
import random


# interpolate flare data onto subexposures
# default is 2 sec
def interp_lc(time, flux, binsize=2, kind='linear', rand_offset = True):

    fin = np.isfinite(time) * np.isfinite(flux)
    f = interp1d(time[fin], flux[fin], kind=kind)

    tmin = time[0]
    tmax = time[-1]

    tot = tmax - tmin ## total time spanned in seconds
    nbins = tot/binsize

    # randomly offset new bins
    if rand_offset:
    	rand = random.uniform(0,binsize)
    else:
    	rand = 0.

    # interpolate data onto evenly spaced subexposures
    xnew = np.linspace(tmin+rand, tmax, num=int(nbins))
    fnew = f(xnew) ## interpolated flux

    return xnew, fnew    

# draw subsamples of 10 data points and apply cosmic ray rejection algorithm
def bin_and_reject_cosmics(time, flux, algorithm='central_n2', n=10):
    start = random.randint(0,n)
    if algorithm=='central_n2':
        i = start
        tout = np.zeros(len(flux[start:])/n)
        fout = np.zeros(len(flux[start:])/n)
        step = 0
        while i <= (len(time)-n):
            x = np.sort(flux[i:i+n]) ## order flux in this bin
            tout[step] = np.mean(time[i:i+n]) ## record time stamp as mean
            fout[step] = np.mean(x[1:-1])*n ## remove lowest and highest and take mean, multiply by n data points
            step += 1
            i = i+n
            
    return tout, fout

# coadd subsamples to get from 2 sec to 2 min or 30 min cadence data
def coadd_to_texp(time, flux, cadence=120):

    # how many data points do I need to average to get to requested cadence?
    # everything here should be an integer...
    ntoaverage, check = divmod(cadence, round(time[2]-time[1]))
    print "avfsdf", ntoaverage, check
    assert check == 0

    # didn't necessarily pass in an array that divides perfectly into requested cadence
    start = random.randint(0,ntoaverage)
    trash = (len(time)-start) % ntoaverage
    print trash, -1*trash

    if trash == 0:
        tout = np.mean(time[start:].reshape(-1, int(ntoaverage)), 1)
        fout = np.sum(flux[start:].reshape(-1, int(ntoaverage)), 1)
    else:
        tout = np.mean(time[start:int(-1*trash)].reshape(-1, int(ntoaverage)), 1)
        fout = np.sum(flux[start:int(-1*trash)].reshape(-1, int(ntoaverage)), 1)
    
    return tout, fout



class LightCurve(NoiseModel.NoiseModel):
    
    def __init__(self, Tmag, Texp, Aper, fracflux=0.8, coord=None):
        NoiseModel.NoiseModel.__init__(self, Tmag, Texp, Aper, fracflux=fracflux, coord=coord, unit='counts')
        self.debug = False
            
    def add_data(self, x, y, xunit='days', yunit='mag', 
                        rebin=True, 
                        reject_cosmics=False,
                        add_noise=False):
        
        if reject_cosmics:
            rebin = True
            saveTexp = self.Texp
            self.Texp = 2 # 2 second needed to do CR rejection
        
        # need to rebin in time to 2 min or 30 min, but worrying 
        # about interpolating isn't import for this purpose because  
        # the flux and noise levels are set from Tmag: and
        # I will assume the same Tmag throughout a given exposure
        
        if xunit=='days':
            xsec = np.array(x-x[0])*24.*60.*60.
        elif xunit=='hours':
            xsec = np.array(x)*60.*60.
        else:
            print "Assuming time is in seconds."
            xsec = np.array(x)
            
        # assume the base Tmag previously provided corresponds to the median flux
        # this program will overwrite the float Tmag with an array
        if yunit=='flux':
            y = self.Tmag - 2.5*np.log10(y/np.nanmedian(y)) ## <-- NEED TO CHECK THIS
        
        if rebin:
            # If reject_cosmics is True, this will interpolate to 2 sec cadence
            print "Rebinning to ", self.Texp
            self.Time, self.Tmag = interp_lc(xsec, y, self.Texp)
            if not reject_cosmics and self.debug:
                plt.plot(xsec, y, label='inpt')
                plt.plot(self.Time, self.Tmag, c='indianred', label='Coadded ({0} sec)'.format(self.Texp))
        else:
            self.Time = np.array(x)
            self.Tmag = np.array(y)
        
        if add_noise:
            noise = self.tot_noise() # noise is in counts
            simdata = np.random.normal(self.star_counts(), noise)
            self.CountsObserved = simdata
          
        else:
            self.CountsObserved = self.star_counts()

        self.TimeObserved = np.array(self.Time)
        self.SignalToNoise = self.star_counts()/self.tot_noise()
        self.DeltaSN = (self.star_counts()-np.min(self.star_counts()))/self.tot_noise()
        if reject_cosmics:
            self.Texp = saveTexp
            if self.debug:
                plt.plot(self.TimeObserved, self.CountsObserved*self.Texp/2, c='k', label='Raw stream (2 sec)')
            
            tout, fout = bin_and_reject_cosmics(self.TimeObserved, self.CountsObserved)
            if self.debug:
                plt.plot(tout, fout*self.Texp/20, '--', c='steelblue', label='Cosmic ray rejected (20 sec)')
                
            tout, fout = coadd_to_texp(tout, fout, self.Texp)
            self.TimeObserved = np.array(tout)
            self.CountsObserved = np.array(fout)
            
            if self.debug:
                plt.plot(self.TimeObserved, self.CountsObserved, c='indianred', label='Coadded ({0} sec)'.format(self.Texp))
                plt.legend()
        
        if rebin:
			if xunit=='days':
				self.Time = self.Time/(24.*60.*60.)
				self.TimeObserved = self.TimeObserved/(24.*60.*60.)
			elif xunit=='hours':
				self.Time = self.Time/(60.*60.)
				self.TimeObserved = self.TimeObserved/(24.*60.*60.)
		
        	

            
    
     