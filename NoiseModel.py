import numpy as np
import math
import matplotlib.pyplot as plt

class NoiseModel():
    
    def __init__(self, Tmag, Texp, Aper, fracflux=0.8, coord=None, unit='ppm'):
        self.Tmag = Tmag # TESS magnitude of target
        self.Texp = Texp # exposure time
        self.Aper = Aper # number of pixels in aperture
        self.FracFlux = fracflux # fraction of total flux enclosed in aperture
        assert (fracflux < 1)
        
        # use coordinates to calculate ecliptic latitute
        if coord is None: 
            self.EcLat = 30.
            
        self.noise_unit = unit
     
    def cts2unit(self):
        
        if self.noise_unit == 'ppm':
            factor = 1./self.star_counts()*1e6
        elif self.noise_unit == 'percent':
            factor = 1./self.star_counts()*1e2
        elif self.noise_unit == 'counts':
            factor = 1.
        elif self.noise_unit == 'fraction':
        	factor = 1./self.star_counts()
        else:
            raise           
        return factor

    def ctsperpix2unit(self):
        
        factor = self.cts2unit() * self.Aper
        return factor

    # Star flux (returns counts if all light included in aperture)        
    def star_counts(self):
        
        # counts for a Tmag = 0 source
        SB0 = 1.514e6 # photons / s / cm2
        F0 = SB0 * 69 # TESS aperture is 69 cm2

        Fstar = F0 * np.power(10, -0.4 * self.Tmag ) # photons / s for star
        return Fstar * self.Texp * self.FracFlux

    # Poisson noise from star (returns counts)
    def star_noise(self):
        # star noise
        Fstar = self.star_counts()
        #return np.sqrt(Fstar) # photons
        return np.sqrt(Fstar) * self.cts2unit() # ppm

    # Verify star noise calculation against plot in Sullivan paper
    def check_star_noise(self):
        
        Tmags = np.linspace(0,16,100)
        plt.plot(Tmags, star_noise(Tmags))
        plt.yscale('log')
        plt.show()
    
    # Read noise (returns e- /pix)
    def read_noise(self):
        
        nexp = self.Texp/2. 
        #print "N exposures", nexp
        #print 10.*np.sqrt(nexp) # e- / pix
        return 10.*np.sqrt(self.Aper*nexp) * self.cts2unit() # ppm
    
    # Systematic noise floor
    def sys_noise(self):
        
        noise_floor = 60 # ppm per hour
        tmp = 60 * np.power(self.Texp/60./60., -0.5) # ppm
        if self.noise_unit == 'ppm':
            sys = tmp
        elif self.noise_unit == 'percent':
            sys = tmp*1e-4
        elif self.noise_unit == 'counts':
        	sys = tmp*1e-6*self.star_counts()
        elif self.noise_unit == 'fraction':
        	sys = tmp*1e-6
        return sys
    
    # Zodiacal Light contribution (returns e- / pix)
    def zodi_noise(self):

        # V = V-band surface brightness of zodiacal light in mag/arcsec2
        Vmax = 23.345
        deltaV = 1.148
        
        V = Vmax - deltaV * (self.EcLat/90 - 1)**2 
        ZL_flux = 2.56e-3 * np.power(10, -0.4*(V-22.8)) # photons / s / cm2 / arcsec2
        ZL_rec = ZL_flux * 69 * 21.1**2 # TESS aperture is 69 cm2, pixels are 21 arcsec/cm

        # should get 95-270 e- / pix in a 2 sec exposure - checks out
        #return ZL_rec*self.Texp # e- / pix    
        return np.sqrt(ZL_rec*self.Texp * self.Aper) * self.cts2unit() # ppm

    
    # Total noise (relative to flux from star)
    def tot_noise(self):
    
        star = self.star_noise()
        read = self.read_noise()
        zodi = self.zodi_noise()
        sys = self.sys_noise()
        
        return star + read + zodi + sys
	
	
	
    # Transit properties
    def add_transit(self, 
			prad=1., # planet radius in Jupiter radii
			srad=1., # stellar radius in solar radii
			smass=1., # stellar mass in solar masses
			period=5., # planet period in days
			b=0 # impact parameter
			):
	
		prad = prad*.10 # convert Jupiter radii to solar radii
		k = prad/srad # transit depth (srad is in solar radii)
		period = period/365. # Convert days to years
		a = np.power(period**2*smass, 1./3)*215. # calculated AU converted to solar radii
		
		T14 = period/math.pi*math.asin( srad/a*np.sqrt( (1.+k)**2. -b**2) )
		
		self.transit_duration = T14*365*24 # now in hrs
		self.transit_depth = k**2
		self.planet_period = period*365
    

    def calc_snr(self, contratio=None, obswindow=27, noisefromcontam=False, pointprecision=False):
	
        # contratio is the TIC contamination ratio (counts neighbors/ counts stars)
        if contratio is None:
			self.dilution = 1.
        else:
			self.dilution = 1./(contratio + 1)
			
        old_texp = self.Texp # save for later
        old_unit = self.noise_unit #
		
        if not pointprecision:
        	self.Texp = self.transit_duration * 60 * 60 # now in seconds
        # else, self.Texp is 2min or 30min.
        self.noise_unit = 'fraction'
        
        if not pointprecision:
        	ntransits = math.floor(obswindow/self.planet_period)
        else:
        	ntransits = 1.
        self.ntransits = ntransits
        
        intrnoise = 60.*1e-6 ## approxmately upper ~70% boundary of distribution for hottest stars in Sullivan Fig. 7
        
        if noisefromcontam:
        	contamnoise = contratio*intrnoise # assume same contribution from contaminants
        	intrnoise = np.sqrt(intrnoise**2 + contamnoise**2)
        	
        noise2 = self.tot_noise()**2 + intrnoise**2
        snr = np.sqrt(ntransits) * self.transit_depth * self.dilution / np.sqrt(noise2)
        
		
        self.Texp = old_texp
        self.noise_unit = old_unit
		
        return snr
