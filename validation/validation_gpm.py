# -*- coding: utf-8 -*-
"""
Created on Fri Jul 22 10:09:02 2016

@author: wolfensb
"""

import glob
import datetime
import numpy as np
import os
import matplotlib.pyplot as plt
import h5py
from cosmo_pol import pycosmo
from cosmo_pol.radar_operator import RadarOperator
from cosmo_pol.gpm.GPM_simulator import compare_operator_with_GPM
from cosmo_pol.radar.read_Swiss_Radar_QPE import read_8bit_meteoswiss
from scipy.interpolate import griddata

FOLDER_COSMO_1MOM = '/ltedata/COSMO/GPM_1MOM/'
FOLDER_COSMO_2MOM = '/ltedata/COSMO/GPM_2MOM/'
FOLDER_GPM = '/media/wolfensb/Storage/Dropbox/GPM_Analysis/DATA/FILES/GPM_DPR/'
FOLDER_RADAR = '/ltedata/Radar_MeteoSwiss/Raw_data/Rain_product/8-bits/'

A_KU = 340.56
B_KU = 1.52
A_KA = 321.74
B_KA = 1.47

AVAIL_MIN = np.asarray([0,2,5,7])

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
if not os.path.exists('./GPM/Measured_Swiss/Ku/'):
    os.makedirs('./GPM/Measured_Swiss/Ku/')   
if not os.path.exists('./GPM/Measured_Swiss/Ka/'):
    os.makedirs('./GPM/Measured_Swiss/Ka/')  
    
for f in files_GPM:
    print(f)
    bname = os.path.basename(f).replace('.HDF5','')
    
#    # 1 MOM
    f_cosmo = FOLDER_COSMO_1MOM + os.path.basename(f).replace('.HDF5','.grb')
    
    if os.path.exists(f_cosmo):
        
        
        # KU band
        ########
        # ZH
        if not os.path.exists('./GPM/Simulated_1mom/Ku/'+bname+'.npy'):
            print(bname)
            s.load_model_file(f_cosmo)
            gpm_s = s.get_GPM_swath(f,'Ku')
            ground,b,c = compare_operator_with_GPM(gpm_s,f)
            np.save('./GPM/Simulated_1mom/Ku/'+bname+'.npy',ground['simulated'])
            np.save('./GPM/Measured/Ku/'+bname+'.npy',ground['measured'])        

#        # Ku
        try:
            if not os.path.exists('./GPM/Measured_Swiss/Ka/'+bname+'.npy'):
                gpm_f=h5py.File(f,'r')    
                lon_gpm=gpm_f['NS']['Longitude'][:]
                lat_gpm=gpm_f['NS']['Latitude'][:]  
            
                date=datetime.datetime.strptime(bname,'%Y-%m-%d-%H-%M')
                closest_match_min=np.argmin(np.abs(date.minute%10-AVAIL_MIN))
            
                subfolder=str(date.year)+'-'+str(date.month).zfill(2)+'/'
                filename='meteoswiss.radar.precip.sv.'+str(date.year)+str(date.month).zfill(2)+str(date.day).zfill(2)+str(date.hour).zfill(2)+str(int(np.floor(date.minute/10)))+str(AVAIL_MIN[closest_match_min])+'.gif'        # Precip
                lat,lon,QPE = read_8bit_meteoswiss(FOLDER_RADAR+subfolder+filename,200000)
                QPE_i = griddata((lat.ravel(),lon.ravel()),QPE.ravel(),(lat_gpm,lon_gpm))
                ZH_rad = 10*np.log10(A_KU * QPE_i**B_KU)
                
                np.save('./GPM/Measured_Swiss/Ku/'+bname+'.npy',ZH_rad)
                
                # Ka
                lon_gpm=gpm_f['HS']['Longitude'][:]
                lat_gpm=gpm_f['HS']['Latitude'][:]    
                QPE_i = griddata((lat.ravel(),lon.ravel()),QPE.ravel(),(lat_gpm,lon_gpm))
                print('ok')
                plt.hist(QPE_i.ravel()[np.isfinite(QPE_i.ravel())],bins=100)
                ZH_rad = 10*np.log10(A_KA * QPE_i**B_KA)        
                np.save('./GPM/Measured_Swiss/Ka/'+bname+'.npy',ZH_rad)
        except:
            pass
    
            
         # KA band
        ########
        # ZH
        if not os.path.exists('./GPM/Simulated_1mom/Ka/'+bname+'.npy'):
            if not len(s.dic_vars):
                s.load_model_file(f_cosmo)
            gpm_s = s.get_GPM_swath(f,'Ka')
            ground,b,c = compare_operator_with_GPM(gpm_s,f)
            np.save('./GPM/Simulated_1mom/Ka/'+bname+'.npy',ground['simulated'])
            np.save('./GPM/Measured/Ka/'+bname+'.npy',ground['measured'])           


    #     
    # 2 MOM
#    f_cosmo = FOLDER_COSMO_2MOM + os.path.basename(f).replace('.HDF5','.grb')
#    print(os.path.exists('./GPM/Simulated_2mom_1mom_MD/Ku/'+bname+'.npy'))
#    if os.path.exists(f_cosmo) and not os.path.exists('/storage/cosmo_pol/validation/GPM/Simulated_2mom_1mom_MD//Ku/'+bname+'.npy'):
#        s.load_model_file(f_cosmo)
#        # KU band
#        #########
#        # ZH
#        gpm_s = s.get_GPM_swath(f,'Ku')
#        ground,b,c = compare_operator_with_GPM(gpm_s,f)
#        np.save('./GPM/Simulated_2mom_1mom_MD/Ku/'+bname+'.npy',ground['simulated'])   
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
#        np.save('./GPM/Simulated_PRECIP_2mom/Ku/'+bname+'.npy',rain_COSMO_interp)
#            
#        # KA band
#        #########
#        # ZH
#        gpm_s = s.get_GPM_swath(f,'Ka')
#        ground,b,c = compare_operator_with_GPM(gpm_s,f)
#        np.save('./GPM/Simulated_2mom_1mom_MD/Ka/'+bname+'.npy',ground['simulated'])     
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