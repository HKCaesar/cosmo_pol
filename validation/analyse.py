# -*- coding: utf-8 -*-
"""
Created on Wed Jul  6 15:28:46 2016

@author: daniel
"""

import glob
import numpy as np
import copy
import matplotlib.pyplot as plt
plt.close('all')
def select(data):
    
    data_copy = copy.deepcopy(data)
    mask_zh = np.logical_and(data['ZH'] < 70,data['ZH']>-10).astype(int)
    mask_zdr = np.logical_and(data['ZDR'] < 10,data['ZDR']>-2).astype(int)
    mask_kdp = np.logical_and(data['KDP'] < 10,data['KDP']>-2).astype(int)
    

    data_copy['ZH'] = data['ZH'][(mask_zh + mask_zdr + mask_kdp )== 3]
    data_copy['ZDR'] = data['ZDR'][(mask_zh + mask_zdr + mask_kdp) == 3]
    data_copy['KDP'] = data['KDP'][(mask_zh + mask_zdr + mask_kdp ) == 3]
    
    return data_copy
    
EVENT = '20150813'
SCHEME = '1mom'

all_files = glob.glob('/storage/cosmo_pol/validation/DATA/CH/'+EVENT+'/radar/ZDR*.npy')

radar = {'ZH':[],'KDP':[],'ZDR':[]}
simul = {'ZH':[],'KDP':[],'ZDR':[]}

for f in all_files:
    try:
        radar['ZH'].extend(np.load(f.replace('ZDR','ZH')).ravel())
        radar['ZDR'].extend(np.load(f).ravel())
        radar['KDP'].extend(np.load(f.replace('ZDR','KDP')).ravel())

        simul['ZH'].extend(np.load(f.replace('radar','model_'+SCHEME).replace('ZDR','ZH')).ravel())
        simul['ZDR'].extend(np.load(f.replace('radar','model_'+SCHEME)).ravel())
        simul['KDP'].extend(np.load(f.replace('radar','model_'+SCHEME).replace('ZDR','KDP')).ravel())
    except:
        
        pass          
        
radar['ZH'] = np.array(radar['ZH'])
radar['ZDR'] = np.array(radar['ZDR'])+0.5
radar['KDP'] = np.array(radar['KDP'])
simul['ZH'] = np.array(simul['ZH'])
simul['ZDR'] = np.array(simul['ZDR'])
simul['KDP'] = np.array(simul['KDP'])

radar2 = select(radar)
simul2 = select(simul)

ZH_bins = np.arange(10,46,2)
plt.figure()
pts_mean_rad = np.zeros((len(ZH_bins)-1,))
pts_stdev_rad = np.zeros((len(ZH_bins)-1,))
pts_mean_simul = np.zeros((len(ZH_bins)-1,))
pts_stdev_simul = np.zeros((len(ZH_bins)-1,))

for j in range(len(ZH_bins)-1):
    mask1 = np.logical_and(radar2['ZH']<ZH_bins[j+1],radar2['ZH']>ZH_bins[j])
    pts_mean_rad[j]=np.nanmedian(radar2['KDP'][mask1])
    pts_stdev_rad[j]=np.nanstd(radar2['KDP'][mask1])
    mask2 = np.logical_and(simul2['ZH']<ZH_bins[j+1],simul2['ZH']>ZH_bins[j])
    pts_mean_simul[j]=np.nanmedian(simul2['KDP'][mask2])
    pts_stdev_simul[j]=np.nanstd(simul2['KDP'][mask2])    
    
plt.plot(0.5*(ZH_bins[0:-1]+ZH_bins[1:]),pts_mean_rad, 'r', 0.5*(ZH_bins[0:-1]+ZH_bins[1:]),pts_mean_simul,'g',linewidth=2)
plt.hold(True) 
plt.plot(0.5*(ZH_bins[0:-1]+ZH_bins[1:]),pts_mean_rad+pts_stdev_rad,'--r')   
plt.plot(0.5*(ZH_bins[0:-1]+ZH_bins[1:]),pts_mean_rad-pts_stdev_rad,'--r')   
plt.plot(0.5*(ZH_bins[0:-1]+ZH_bins[1:]),pts_mean_simul+pts_stdev_simul,'--g')   
plt.plot(0.5*(ZH_bins[0:-1]+ZH_bins[1:]),pts_mean_simul-pts_stdev_simul,'--g')   
plt.grid()
plt.xlabel('ZH [dBZ]')
plt.title('Heavy rainfall event: 13 August 2015')
plt.ylabel(r'Average KDP within ZH bin [$^{\circ}$/km]')
plt.legend(['Mt. Lema radar','Simulation'])
plt.savefig('./COMPARISON/CH/'+SCHEME+'/KDPvsZH_hist_'+EVENT+'.pdf',dpi=200)

plt.figure()
pts_mean_rad = np.zeros((len(ZH_bins)-1,))
pts_stdev_rad = np.zeros((len(ZH_bins)-1,))
pts_mean_simul = np.zeros((len(ZH_bins)-1,))
pts_stdev_simul = np.zeros((len(ZH_bins)-1,))

for j in range(len(ZH_bins)-1):
    mask1 = np.logical_and(radar2['ZH']<ZH_bins[j+1],radar2['ZH']>ZH_bins[j])
    pts_mean_rad[j]=np.nanmedian(radar2['ZDR'][mask1])
    pts_stdev_rad[j]=np.nanstd(radar2['ZDR'][mask1])
    mask2 = np.logical_and(simul2['ZH']<ZH_bins[j+1],simul2['ZH']>ZH_bins[j])
    pts_mean_simul[j]=np.nanmedian(simul2['ZDR'][mask2])
    pts_stdev_simul[j]=np.nanstd(simul2['ZDR'][mask2])    
    
plt.plot(0.5*(ZH_bins[0:-1]+ZH_bins[1:]),pts_mean_rad, 'r', 0.5*(ZH_bins[0:-1]+ZH_bins[1:]),pts_mean_simul,'g',linewidth=2)
plt.hold(True) 
plt.plot(0.5*(ZH_bins[0:-1]+ZH_bins[1:]),pts_mean_rad+pts_stdev_rad,'--r')   
plt.plot(0.5*(ZH_bins[0:-1]+ZH_bins[1:]),pts_mean_rad-pts_stdev_rad,'--r')   
plt.plot(0.5*(ZH_bins[0:-1]+ZH_bins[1:]),pts_mean_simul+pts_stdev_simul,'--g')   
plt.plot(0.5*(ZH_bins[0:-1]+ZH_bins[1:]),pts_mean_simul-pts_stdev_simul,'--g')   
plt.grid()
plt.xlabel('ZH [dBZ]')
plt.ylabel(r'Average ZDR within ZH bin [dB]')
plt.title('Heavy rainfall event: 13 August 2015')
plt.legend(['Mt. Lema radar','Simulation'])
plt.savefig('./COMPARISON/CH/'+SCHEME+'/ZDRvsZH_hist_'+EVENT+'.pdf',dpi=200)


plt.figure()
plt.hist(radar2['KDP'],bins=50,alpha=0.5,normed=True)
plt.hist(2*simul2['KDP'],bins=50,alpha=0.5,normed=True)
plt.xlim([-1,5])
plt.legend(['Radar','Simulation'])
plt.grid()
plt.title('Heavy rainfall event: 13 August 2015')
plt.ylabel('Density')
plt.xlabel(r'KDP [$^{\circ}$/km]')
plt.savefig('./COMPARISON/CH/'+SCHEME+'/KDP_hist_'+EVENT+'.pdf',dpi=200,bbox_inches='tight')


plt.figure()
plt.hist(radar2['ZDR'],bins=50,alpha=0.5,normed=True)
plt.hist(simul2['ZDR'],bins=50,alpha=0.5,normed=True)
plt.xlim([-1,5])
plt.legend(['Mt. Lema radar','Simulation'])
plt.grid()
plt.title('Heavy rainfall event: 13 Aug. 2015')
plt.ylabel('Density')
plt.xlabel(r'ZDR [dB]')
plt.savefig('./COMPARISON/CH/'+SCHEME+'/ZDR_hist_'+EVENT+'.pdf',dpi=200,bbox_inches='tight')

plt.figure()
plt.hist(radar2['ZH'],bins=50,alpha=0.5,normed=True)
plt.hist(simul2['ZH'],bins=50,alpha=0.5,normed=True)
plt.legend(['Mt. Lema radar','Simulation'])
plt.grid()
plt.title('Heavy rainfall event: 13 Aug. 2015')
plt.ylabel('Density')
plt.xlabel(r'ZH [dBZ]')
plt.savefig('./COMPARISON/CH/'+SCHEME+'/ZH_hist_'+EVENT+'.pdf',dpi=200,bbox_inches='tight')

