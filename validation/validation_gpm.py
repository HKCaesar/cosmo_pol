# -*- coding: utf-8 -*-
"""
Created on Fri Jul 22 10:09:02 2016

@author: wolfensb
"""

import glob
import numpy as np
import os
import h5py
from cosmo_pol import pycosmo
from cosmo_pol.radar_operator import RadarOperator
from cosmo_pol.gpm.GPM_simulator import compare_operator_with_GPM
from scipy.interpolate import griddata

FOLDER_COSMO_1MOM = '/ltedata/COSMO/GPM_1MOM/'
FOLDER_COSMO_2MOM = '/ltedata/COSMO/GPM_2MOM/'
FOLDER_GPM = '/media/wolfensb/Storage/Dropbox/GPM_Analysis/DATA/FILES/GPM_DPR/'



files_GPM = glob.glob(FOLDER_GPM + '*.HDF5')

s = RadarOperator()

# Create output folders

if not os.path.exists('./GPM/Measured/Ka/'):
    os.makedirs('./GPM/Measured/Ka/')
if not os.path.exists('./GPM/Measured/Ku/'):
    os.makedirs('./GPM/Measured/Ku/')
if not os.path.exists('./GPM/Simulated_1mom/Ka/'):
    os.makedirs('./GPM/Simulated_1mom/Ka/')
if not os.path.exists('./GPM/Simulated_1mom/Ku/'):
    os.makedirs('./GPM/Simulated_1mom/Ku/')
if not os.path.exists('./GPM/Simulated_2mom/Ka/'):
    os.makedirs('./GPM/Simulated_2mom/Ka/')
if not os.path.exists('./GPM/Simulated_2mom/Ku/'):
    os.makedirs('./GPM/Simulated_2mom/Ku/') 
if not os.path.exists('./GPM/Simulated_PRECIP_1mom/Ku/'):
    os.makedirs('./GPM/Simulated_PRECIP_1mom/Ku/')   
if not os.path.exists('./GPM/Simulated_PRECIP_1mom/Ka/'):
    os.makedirs('./GPM/Simulated_PRECIP_1mom/Ka/')   
if not os.path.exists('./GPM/Simulated_PRECIP_2mom/Ku/'):
    os.makedirs('./GPM/Simulated_PRECIP_2mom/Ku/')   
if not os.path.exists('./GPM/Simulated_PRECIP_2mom/Ka/'):
    os.makedirs('./GPM/Simulated_PRECIP_2mom/Ka/')  
    
for f in files_GPM:
    bname = os.path.basename(f).strip('.HDF5')
    
    # 1 MOM
    f_cosmo = FOLDER_COSMO_1MOM + os.path.basename(f).replace('.HDF5','.grb')
#    s.load_model_file(f_cosmo)
    if os.path.exists(f_cosmo):
        
        # KU band
        #########
        # ZH
#        gpm_s = s.get_GPM_swath(f,'Ku')
#        ground,b,c = compare_operator_with_GPM(gpm_s,f)
#        np.save('./GPM/Simulated_1mom/Ku/'+bname+'.npy',ground['simulated'])
#        np.save('./GPM/Measured/Ku/'+bname+'.npy',ground['measured'])        
        
        # Precip
        # Get precip
        COSMO = pycosmo.open_file(f_cosmo)
        if pycosmo.check_if_variables_in_file(COSMO,'PRR_GSP_GDS10_SFC'):
            rain_COSMO=pycosmo.get_variable(COSMO,'PREC_RATE')*3600
            lat_COSMO = rain_COSMO.coordinates['lat_2D']
            lon_COSMO = rain_COSMO.coordinates['lon_2D']        
            gpm = h5py.File(f,'r')
            lats = gpm['NS']['Latitude']
            lons = gpm['NS']['Longitude']    
        
            rain_COSMO_interp=griddata((lat_COSMO.ravel(),lon_COSMO.ravel()),rain_COSMO.data.ravel(),(lats,lons))
            
            np.save('./GPM/Simulated_PRECIP_1mom/Ku/'+bname+'.npy',rain_COSMO_interp)
            
        # KA band
        #########
        # ZH
#        gpm_s = s.get_GPM_swath(f,'Ka')
#        ground,b,c = compare_operator_with_GPM(gpm_s,f)
#        np.save('./GPM/Simulated_1mom/Ka/'+bname+'.npy',ground['simulated'])
#        np.save('./GPM/Measured/Ka/'+bname+'.npy',ground['measured'])           

        # Precip
        # Get precip    
        if pycosmo.check_if_variables_in_file(COSMO,'PRR_GSP_GDS10_SFC'):
            gpm = h5py.File(f,'r')
            lats = gpm['HS']['Latitude']
            lons = gpm['HS']['Longitude']    
            
            rain_COSMO_interp=griddata((lat_COSMO.ravel(),lon_COSMO.ravel()),rain_COSMO.data.ravel(),(lats,lons))
            
            np.save('./GPM/Simulated_PRECIP_1mom/Ka/'+bname+'.npy',rain_COSMO_interp)
        
        
        
    # TODO!
#        
#    # 2 MOM
#    f_cosmo = FOLDER_COSMO_2MOM + os.path.basename(f).replace('.HDF5','.grb')
#    s.load_model_file(f_cosmo)
#    if os.path.exists(f_cosmo):
#
#        # KU band
#        #########
#        # ZH
##        gpm_s = s.get_GPM_swath(f,'Ku')
##        ground,b,c = compare_operator_with_GPM(gpm_s,f)
##        np.save('./GPM/Simulated_2mom/Ku/'+bname+'.npy',ground['simulated'])   
#        
#        # Precip
#        # Get precip
#        COSMO = pycosmo.open_file(f_cosmo)
#        rain_COSMO=pycosmo.get_variable(COSMO,'PREC_RATE')*3600
#        lat_COSMO = rain_COSMO.coordinates['lat_2D']
#        lon_COSMO = rain_COSMO.coordinates['lon_2D']        
#        gpm = h5py.File(f,'r')
#        lats = gpm['NS']['Latitude']
#        lons = gpm['NS']['Longitude']    
#        
#        rain_COSMO_interp=griddata((lat_COSMO.ravel(),lon_COSMO.ravel()),rain_COSMO.data.ravel(),(lats,lons))
#        
#        np.save('./GPM/Simulated_PRECIP_1mom/Ku/'+bname+'.npy',rain_COSMO_interp)
#            
#        # KA band
#        #########
#        # ZH
##        gpm_s = s.get_GPM_swath(f,'Ka')
##        ground,b,c = compare_operator_with_GPM(gpm_s,f)
##        np.save('./GPM/Simulated_2mom/Ka/'+bname+'.npy',ground['simulated'])     
#
#        # Precip
#        # Get precip    
#        gpm = h5py.File(f,'r')
#        lats = gpm['HS']['Latitude']
#        lons = gpm['HS']['Longitude']    
#        
#        rain_COSMO_interp=griddata((lat_COSMO.ravel(),lon_COSMO.ravel()),rain_COSMO.data.ravel(),(lats,lons))
#        
#        np.save('./GPM/Simulated_PRECIP_2mom/Ka/'+bname+'.npy',rain_COSMO_interp)
#        
#        
#        