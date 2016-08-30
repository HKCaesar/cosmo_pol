# -*- coding: utf-8 -*-
"""
Created on Wed Aug 17 12:23:24 2016

@author: wolfensb
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats    

f = open('power_laws_ar','w')


plt.close('all')

agg = np.genfromtxt('stats_aggregatess_Davos_winter_2015-16.txt')
agg[:,2] = 1.0/agg[:,2]
grau = np.genfromtxt('stats_graupels_Davos_winter_2015-16.txt')
grau[:,2] = 1.0/grau[:,2]

agg_d_bins = np.arange(2.4,15.2,0.1)

ar = np.zeros((len(agg_d_bins)-1,))
ar_q10 = np.zeros((len(agg_d_bins)-1,))
ar_q90 = np.zeros((len(agg_d_bins)-1,))
ar_q25 = np.zeros((len(agg_d_bins)-1,))
ar_q75 = np.zeros((len(agg_d_bins)-1,))

ar_alpha = np.zeros((len(agg_d_bins)-1,))
ar_loc = np.zeros((len(agg_d_bins)-1,))
ar_scale = np.zeros((len(agg_d_bins)-1,))

for i in range(len(ar)):
    ar[i] = np.nanmedian(agg[:,2][np.logical_and(agg[:,1]<agg_d_bins[i+1],agg[:,1]>agg_d_bins[i])])
    ar_q10[i] = np.percentile((agg[:,2][np.logical_and(agg[:,1]<agg_d_bins[i+1],agg[:,1]>agg_d_bins[i])]),10)
    ar_q90[i] = np.percentile((agg[:,2][np.logical_and(agg[:,1]<agg_d_bins[i+1],agg[:,1]>agg_d_bins[i])]),90)
    ar_q25[i] = np.percentile((agg[:,2][np.logical_and(agg[:,1]<agg_d_bins[i+1],agg[:,1]>agg_d_bins[i])]),25)
    ar_q75[i] = np.percentile((agg[:,2][np.logical_and(agg[:,1]<agg_d_bins[i+1],agg[:,1]>agg_d_bins[i])]),75)    
    
    ar_alpha[i],ar_loc[i], ar_scale[i] = stats.gamma.fit(agg[:,2][np.logical_and(agg[:,1]<agg_d_bins[i+1],agg[:,1]>agg_d_bins[i])]-1)

ar_loc[ar_loc<1.]=1.

bins_centers = (agg_d_bins[0:-1]+ agg_d_bins[1:])/2.

fit_alpha = np.polyfit(np.log(bins_centers),np.log(ar_alpha),1)
fit_loc = np.polyfit(np.log(bins_centers),np.log(ar_loc),1)
fit_scale = np.polyfit(np.log(bins_centers),np.log(ar_scale),1)

f.write('Aggregates: alpha: a = '+str(np.exp(fit_alpha[1]))+', b = '+str(fit_alpha[0])+' \n') # python will convert \n to os.linesep
f.write('Aggregates: loc: a = '+str(np.exp(fit_loc[1]))+', b = '+str(fit_loc[0])+' \n') # python will convert \n to os.linesep
f.write('Aggregates: scale: a = '+str(np.exp(fit_scale[1]))+', b = '+str(fit_scale[0])+' \n') # python will convert \n to os.linesep
#
plt.plot(bins_centers,ar,'--b',bins_centers,ar_q25,'--g',bins_centers,ar_q75,'--g',bins_centers,ar_q10,'--r',bins_centers,ar_q90,'--r',linewidth=1.5)

plt.hold(True)
ar_alpha_fit = lambda D: np.exp(fit_alpha[1])*D**fit_alpha[0]
ar_loc_fit = lambda D: np.exp(fit_loc[1])*D**fit_loc[0]
ar_scale_fit = lambda D: np.exp(fit_scale[1])*D**fit_scale[0]
#
#
ar_bins = np.arange(0.6,3.5,0.01)

ar_pdf  = np.zeros((len(ar_bins),len(bins_centers)))
ar_med =  np.zeros((len(agg_d_bins)-1,))
ar_q10 =  np.zeros((len(agg_d_bins)-1,))
ar_q90 =  np.zeros((len(agg_d_bins)-1,))
ar_q25 =  np.zeros((len(agg_d_bins)-1,))
ar_q75 =  np.zeros((len(agg_d_bins)-1,))

for i,d in enumerate(bins_centers):

    gamm = stats.gamma(ar_alpha_fit(d),ar_loc_fit(d),ar_scale_fit(d))
    ar_pdf[:,i] = gamm.pdf(ar_bins)
    ar_med[i] = gamm.median()
    ar_q10[i] = gamm.ppf(0.1)
    ar_q90[i] = gamm.ppf(0.9)
    ar_q25[i] = gamm.ppf(0.25)
    ar_q75[i] = gamm.ppf(0.75)
                
plt.plot(bins_centers,ar_med,bins_centers,ar_q25,'g',bins_centers,ar_q75,'g',bins_centers,ar_q10,'r',bins_centers,ar_q90,'r',linewidth=2)
x,y = np.meshgrid(bins_centers,ar_bins)
plt.contourf(x,y,ar_pdf,levels=np.linspace(np.nanmin(ar_pdf),np.nanmax(ar_pdf),50))
plt.colorbar(label='Density [-]')

plt.figure()
i=20
plt.hist(agg[:,2][np.logical_and(agg[:,1]<agg_d_bins[i+1],agg[:,1]>agg_d_bins[i])],normed=True,bins=50)
gamm = stats.gamma(ar_alpha_fit(bins_centers[i]),ar_loc_fit(bins_centers[i]),ar_scale_fit(bins_centers[i]))
plt.plot(ar_bins,gamm.pdf(ar_bins))


####################################

#plt.figure()
#
#grau_d_bins = np.arange(0.8,5.,0.05)
#
#ar = np.zeros((len(grau_d_bins)-1,))
#ar_q10 = np.zeros((len(grau_d_bins)-1,))
#ar_q90 = np.zeros((len(grau_d_bins)-1,))
#ar_q25 = np.zeros((len(grau_d_bins)-1,))
#ar_q75 = np.zeros((len(grau_d_bins)-1,))
#
#ar_alpha = np.zeros((len(grau_d_bins)-1,))
#ar_loc = np.zeros((len(grau_d_bins)-1,))
#ar_scale = np.zeros((len(grau_d_bins)-1,))
#
#for i in range(len(ar)):
#    ar[i] = np.nanmedian(grau[:,2][np.logical_and(grau[:,1]<grau_d_bins[i+1],grau[:,1]>grau_d_bins[i])])
#    ar_q10[i] = np.percentile((grau[:,2][np.logical_and(grau[:,1]<grau_d_bins[i+1],grau[:,1]>grau_d_bins[i])]),10)
#    ar_q90[i] = np.percentile((grau[:,2][np.logical_and(grau[:,1]<grau_d_bins[i+1],grau[:,1]>grau_d_bins[i])]),90)
#    ar_q25[i] = np.percentile((grau[:,2][np.logical_and(grau[:,1]<grau_d_bins[i+1],grau[:,1]>grau_d_bins[i])]),25)
#    ar_q75[i] = np.percentile((grau[:,2][np.logical_and(grau[:,1]<grau_d_bins[i+1],grau[:,1]>grau_d_bins[i])]),75)    
#    
#    ar_alpha[i],ar_loc[i], ar_scale[i] = stats.gamma.fit(grau[:,2][np.logical_and(grau[:,1]<grau_d_bins[i+1],grau[:,1]>grau_d_bins[i])]-1)
#
#ar_loc[ar_loc<1.]=1.
#
#bins_centers = (grau_d_bins[0:-1]+ grau_d_bins[1:])/2.
#
#fit_alpha = np.polyfit(np.log(bins_centers),np.log(ar_alpha),1)
#fit_loc = np.polyfit(np.log(bins_centers),np.log(ar_loc),1)
#fit_scale = np.polyfit(np.log(bins_centers),np.log(ar_scale),1)
#
#f.write('graupel: alpha: a = '+str(np.exp(fit_alpha[1]))+', b = '+str(fit_alpha[0])+' \n') # python will convert \n to os.linesep
#f.write('graupel: loc: a = '+str(np.exp(fit_loc[1]))+', b = '+str(fit_loc[0])+' \n') # python will convert \n to os.linesep
#f.write('graupel: scale: a = '+str(np.exp(fit_scale[1]))+', b = '+str(fit_scale[0])+' \n') # python will convert \n to os.linesep
##
#plt.plot(bins_centers,ar,'--b',bins_centers,ar_q25,'--g',bins_centers,ar_q75,'--g',bins_centers,ar_q10,'--r',bins_centers,ar_q90,'--r',linewidth=1.5)
#
#plt.hold(True)
#ar_alpha_fit = lambda D: np.exp(fit_alpha[1])*D**fit_alpha[0]
#ar_loc_fit = lambda D: np.exp(fit_loc[1])*D**fit_loc[0]
#ar_scale_fit = lambda D: np.exp(fit_scale[1])*D**fit_scale[0]
##
##
#ar_bins = np.arange(0.6,3.5,0.01)
#
#ar_pdf  = np.zeros((len(ar_bins),len(bins_centers)))
#ar_med =  np.zeros((len(grau_d_bins)-1,))
#ar_q10 =  np.zeros((len(grau_d_bins)-1,))
#ar_q90 =  np.zeros((len(grau_d_bins)-1,))
#ar_q25 =  np.zeros((len(grau_d_bins)-1,))
#ar_q75 =  np.zeros((len(grau_d_bins)-1,))
#
#for i,d in enumerate(bins_centers):
#
#    gamm = stats.gamma(ar_alpha_fit(d),ar_loc_fit(d),ar_scale_fit(d))
#    ar_pdf[:,i] = gamm.pdf(ar_bins)
#    ar_med[i] = gamm.median()
#    ar_q10[i] = gamm.ppf(0.1)
#    ar_q90[i] = gamm.ppf(0.9)
#    ar_q25[i] = gamm.ppf(0.25)
#    ar_q75[i] = gamm.ppf(0.75)
#                
#plt.plot(bins_centers,ar_med,bins_centers,ar_q25,'g',bins_centers,ar_q75,'g',bins_centers,ar_q10,'r',bins_centers,ar_q90,'r',linewidth=2)
#x,y = np.meshgrid(bins_centers,ar_bins)
#plt.contourf(x,y,ar_pdf,levels=np.linspace(np.nanmin(ar_pdf),np.nanmax(ar_pdf),50))
#plt.colorbar(label='Density [-]')
#
#plt.figure()
#i=20
#plt.hist(grau[:,2][np.logical_and(grau[:,1]<grau_d_bins[i+1],grau[:,1]>grau_d_bins[i])],normed=True,bins=50)
#gamm = stats.gamma(ar_alpha_fit(bins_centers[i]),ar_loc_fit(bins_centers[i]),ar_scale_fit(bins_centers[i]))
#plt.plot(ar_bins,gamm.pdf(ar_bins))
#

f.close()