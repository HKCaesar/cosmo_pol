# -*- coding: utf-8 -*-
"""
Created on Tue Aug  4 12:45:21 2015

@author: wolfensb
"""

from cosmo_pol.utilities import cfg
import numpy as np

##############################################################################
# Parameters
##############################################################################

class Constant_class(object):
    def __init__(self):
        
        #,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
        # Numerical parameters
        self.EPS = np.finfo(np.float32).eps
        self.SIMULATED_VARIABLES = ['ZH','DSPECTRUM','RVEL','ZV','PHIDP','ZDR','RHOHV','KDP']
        
        #,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
        # COSMO parameters
        self.MAX_HEIGHT_COSMO=35000 # Maximum height of the COSMO domain
        
        #,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
        # Physical parameters
        self.C=299792458. # m/s
        self.RHO_W=1000./(1000**3) # kg/mm3 at 10°C
        self.RHO_I=916./(1000**3)  # kg/mm3 at 0°C
        self.M_AIR = 1 # dielectric constant of air
        self.VALID_D_RANGE=[0.2,10.] # range of valid hydrometeor diameters
        self.FREEZING_LEVEL=273.15 # Water freezing temperature
        self.A=1.6 # Average value of the turbulence range proposed by Doviak and Zrnic (p. 409)
        self.KW = 0.93 # Dielectric factor for water at weather radar frequencies
        self.RHO_0 = 1.225 # Air density at ground
        
        if cfg.CONFIG != {}:
            #,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
            # Radar parameters
            self.NVEL=self.C/(4*1E-6*cfg.CONFIG['radar']['PRI']*cfg.CONFIG['radar']['frequency']*1E09) # Nyquist velocity
            self.VRES=2*self.NVEL/cfg.CONFIG['radar']['FFT_length'] # size of velocity bin
            self.VARRAY=np.arange(-self.NVEL,self.NVEL+self.VRES,self.VRES) # array of velocity bins
            self.WAVELENGTH=self.C/(cfg.CONFIG['radar']['frequency']*1E09)*1000 # in mm
            self.PULSE_WIDTH = 2*cfg.CONFIG['radar']['radial_resolution']
            try:
                if len(cfg.CONFIG['radar']['sensitivity']) == 3:
                    self.RADAR_CONSTANT_DB = 180 -10*np.log10(np.pi**3 * cfg.CONFIG['radar']['3dB_beamwidth']**2 \
                                        * self.KW**2*self.PULSE_WIDTH/1000) - 2*cfg.CONFIG['radar']['sensitivity'][1]\
                                        + 10*np.log10(1024 * np.log(2) * (self.WAVELENGTH/1000.)**2) # Radar equation, Probert-Jones
            except:
                # If it fails we won't need it anyway
                pass
            if cfg.CONFIG['radar']['type']=='ground':
                self.RANGE_RADAR=np.arange(cfg.CONFIG['radar']['radial_resolution']/2,
                                      cfg.CONFIG['radar']['range'], 
                                      cfg.CONFIG['radar']['radial_resolution'])
                
            self.GPM_SN_RATIO=14.5 # Toyoshima et al, 2015
            self.GPM_RADIAL_RES=125
            
    def update(self):
        self.__init__()
        
constants = Constant_class()