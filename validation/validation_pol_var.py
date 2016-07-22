# -*- coding: utf-8 -*-
"""
Created on Mon Jul  4 11:06:56 2016

@author: wolfensb
"""

from cosmo_pol.radar import small_radar_db, pyart_wrapper
from cosmo_pol.radar.pyart_wrapper import PyradRadop,PyradRadopVProf, RadarDisplay
from cosmo_pol.radar_operator import RadarOperator
from cosmo_pol import pycosmo

import numpy as np
import numpy.ma as ma
import glob, os
import matplotlib.pyplot as plt

EVENT = 'case2015081312_ONEMOM'

FOLDER_RADAR = '/storage/cosmo_pol/validation/data/' + EVENT[4:12]+'/radar/'
FOLDER_1MOM_MODEL = '/storage/cosmo_pol/validation/data/' + EVENT[4:12]+'/model_1mom/'
FOLDER_2MOM_MODEL = '/storage/cosmo_pol/validation/data/' + EVENT[4:12]+'/model_2mom/'

if not os.path.exists(FOLDER_RADAR):
    os.makedirs(FOLDER_RADAR)
if not os.path.exists(FOLDER_1MOM_MODEL):
    os.makedirs(FOLDER_1MOM_MODEL)
if not os.path.exists(FOLDER_2MOM_MODEL):
    os.makedirs(FOLDER_2MOM_MODEL)
    
files_c = pycosmo.get_model_filenames('/ltedata/COSMO/Validation_operator/'+EVENT+'/')
radop = RadarOperator(options_file='../option_files/MXPOL_PPI_dole.yml', diagnostic_mode=True) 
rad_db = aaaa = small_radar_db.CH_RADAR_db()

for i,f in enumerate(files_c['h'][0:]):
    print(i)
    if not os.path.exists(FOLDER_RADAR+'ZH_'+str(i)+'.npy'):
        time = pycosmo.get_time_from_COSMO_filename(f)
            
        # Radar
        files = rad_db.query(date=[str(time)],radar='D',res='L',angle=1)
        
        rad = pyart_wrapper.PyradCH(files[0][0].rstrip(),True)
        rad.snr_threshold(5)
#        rad.average(n_gates=6)
        rad.correct_velocity()
        rad.estimate_KDP()
        
        KDP_rad = ma.filled(rad.get_field(0,'KDP'),np.nan)
        ZH_rad = ma.filled(rad.get_field(0,'Z'),np.nan)
        ZDR_rad = ma.filled(rad.get_field(0,'ZDR'),np.nan)
        RHOHV_rad = ma.filled(rad.get_field(0,'ZDR'),np.nan)
                        
        np.save(FOLDER_RADAR+'ZH_'+str(i)+'.npy',ZH_rad)
        np.save(FOLDER_RADAR+'ZDR_'+str(i)+'.npy',ZDR_rad)
        np.save(FOLDER_RADAR+'KDP_'+str(i)+'.npy',KDP_rad)
        np.save(FOLDER_RADAR+'RHOHV_'+str(i)+'.npy',RHOHV_rad)            

    if not os.path.exists(FOLDER_1MOM_MODEL+'ZH_'+str(i)+'.npy'):        
        try:
    #        display = RadarDisplay(rad, shift=(0.0, 0.0)) 
    #        plt.figure()
    #        display.plot('Z', 0, 150000,colorbar_flag=True,title="ZH (radar)",vmin=0,mask_outside = True)
        #    plt.xlim([-150,150])
        #    plt.ylim([-150,150])
        #    
            # Simulation one mom
            radop.load_model_file(f,cfilename = files_c['c'][0])
            simul = radop.get_PPI(1)
            
            KDP_simul = ma.filled(simul.get_field(0,'KDP'),np.nan)
            ZH_simul = ma.filled(simul.get_field(0,'ZH'),np.nan)
            ZDR_simul = ma.filled(simul.get_field(0,'ZDR'),np.nan)    
            RHOHV_simul = ma.filled(simul.get_field(0,'RHOHV'),np.nan) 
    #          
    #        display = RadarDisplay(simul, shift=(0.0, 0.0)) 
    #        plt.figure()
    #        display.plot('ZH', 0, 150000,colorbar_flag=True,title="ZH (radar)",vmin=0,mask_outside = True)   
          
            np.save(FOLDER_1MOM_MODEL+'ZH_'+str(i)+'.npy',ZH_simul)
            np.save(FOLDER_1MOM_MODEL+'ZDR_'+str(i)+'.npy',ZDR_simul)
            np.save(FOLDER_1MOM_MODEL+'KDP_'+str(i)+'.npy',KDP_simul)
            np.save(FOLDER_1MOM_MODEL+'RHOHV_'+str(i)+'.npy',RHOHV_simul)            
        except:
            pass
    
#    if not os.path.exists(FOLDER_2MOM_MODEL+'ZH_'+str(i)+'.npy'):        
#        try:
#            f = f.replace('ONEMOM','TWOMOM')
#    #        display = RadarDisplay(rad, shift=(0.0, 0.0)) 
#    #        plt.figure()
#    #        display.plot('Z', 0, 150000,colorbar_flag=True,title="ZH (radar)",vmin=0,mask_outside = True)
#        #    plt.xlim([-150,150])
#        #    plt.ylim([-150,150])
#        #    
#            # Simulation one mom
#            radop.load_model_file(f,cfilename = files_c['c'][0])
#            simul = radop.get_PPI(1)
#            
#            KDP_simul = ma.filled(simul.get_field(0,'KDP'),np.nan)
#            ZH_simul = ma.filled(simul.get_field(0,'ZH'),np.nan)
#            ZDR_simul = ma.filled(simul.get_field(0,'ZDR'),np.nan)    
#            RHOHV_simul = ma.filled(simul.get_field(0,'RHOHV'),np.nan) 
#    #          
#    #        display = RadarDisplay(simul, shift=(0.0, 0.0)) 
#    #        plt.figure()
#    #        display.plot('ZH', 0, 150000,colorbar_flag=True,title="ZH (radar)",vmin=0,mask_outside = True)   
#          
#            np.save(FOLDER_2MOM_MODEL+'ZH_'+str(i)+'.npy',ZH_simul)
#            np.save(FOLDER_2MOM_MODEL+'ZDR_'+str(i)+'.npy',ZDR_simul)
#            np.save(FOLDER_2MOM_MODEL+'KDP_'+str(i)+'.npy',KDP_simul)
#            np.save(FOLDER_2MOM_MODEL+'RHOHV_'+str(i)+'.npy',RHOHV_simul)            
#        except:
#            pass    
#        