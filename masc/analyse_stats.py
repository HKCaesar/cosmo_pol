# -*- coding: utf-8 -*-
"""
Created on Wed Aug 17 12:23:24 2016

@author: wolfensb
"""

import numpy as np

plt.close('all')
agg = np.genfromtxt('stats_aggregatess_Davos_winter_2015-16.txt')

grau = np.genfromtxt('stats_graupels_Davos_winter_2015-16.txt')


agg_d_bins = np.arange(2.4,15.2,0.1)

ar = np.zeros((len(agg_d_bins)-1,))
ar_std = np.zeros((len(agg_d_bins)-1,))

for i in range(len(ar)):
    ar[i] = np.nanmean(agg[:,2][np.logical_and(agg[:,1]<agg_d_bins[i+1],agg[:,1]>agg_d_bins[i])])
    ar_std[i] = np.nanstd(agg[:,2][np.logical_and(agg[:,1]<agg_d_bins[i+1],agg[:,1]>agg_d_bins[i])])


bins_centers = (agg_d_bins[0:-1]+ agg_d_bins[1:])/2.

fit_mean = np.polyfit(np.log(bins_centers),np.log(ar),1)
fit_std = np.polyfit(np.log(bins_centers),np.log(ar_std),1)

plt.plot(bins_centers,ar,'--b',bins_centers,ar+ar_std,'--g',bins_centers,ar-ar_std,'--r',linewidth=1.5)
plt.hold(True)
ar_fit = lambda D: np.exp(fit_mean[1])*D**fit_mean[0]
ar_std_fit = lambda D: np.exp(fit_std[1])*D**fit_std[0]
plt.plot(bins_centers,ar_fit(bins_centers),'b',bins_centers,ar_fit(bins_centers)+ar_std_fit(bins_centers),'g',bins_centers,ar_fit(bins_centers)-ar_std_fit(bins_centers),'r',linewidth=1.5)
plt.grid()
plt.xlabel('Agg. diameter [mm]')
plt.ylabel('Axis ratio [-]')
plt.title('Aggregates')



ar_bins = np.arange(0.3,1.3,0.01)

ar_pdf  = np.zeros((len(ar_bins),len(bins_centers)))

for i,d in enumerate(bins_centers):
    sigm = ar_std_fit(d)
    mu = ar_fit(d)
    func = lambda ar: np.exp(-(ar-mu)**2/(2*sigm**2))/np.sqrt(2*np.pi*sigm**2)
    ar_pdf[:,i] = func(ar_bins)
    
x,y = np.meshgrid(bins_centers,ar_bins)
plt.contourf(x,y,ar_pdf,levels=np.linspace(np.nanmin(ar_pdf),np.nanmax(ar_pdf),50))
plt.colorbar(label='Density [-]')

grau_d_bins = np.arange(0.8,5.,0.05)

ar = np.zeros((len(grau_d_bins)-1,))
ar_std = np.zeros((len(grau_d_bins)-1,))

for i in range(len(ar)):
    ar[i] = np.nanmean(grau[:,2][np.logical_and(grau[:,1]<grau_d_bins[i+1],grau[:,1]>grau_d_bins[i])])
    ar_std[i] = np.nanstd(grau[:,2][np.logical_and(grau[:,1]<grau_d_bins[i+1],grau[:,1]>grau_d_bins[i])])


bins_centers = (grau_d_bins[0:-1]+ grau_d_bins[1:])/2.
fit_mean = np.polyfit(np.log(bins_centers),np.log(ar),1)
fit_std = np.polyfit(np.log(bins_centers),np.log(ar_std),1)

plt.figure()
plt.plot(bins_centers,ar,'--b',bins_centers,ar+ar_std,'--g',bins_centers,ar-ar_std,'--r',linewidth=1.5)
plt.hold(True)
ar_fit = lambda D: np.exp(fit_mean[1])*D**fit_mean[0]
ar_std_fit = lambda D: np.exp(fit_std[1])*D**fit_std[0]
plt.plot(bins_centers,ar_fit(bins_centers),'b',bins_centers,ar_fit(bins_centers)+ar_std_fit(bins_centers),'g',bins_centers,ar_fit(bins_centers)-ar_std_fit(bins_centers),'r',linewidth=1.5)
plt.grid()
plt.xlabel('Agg. diameter [mm]')
plt.ylabel('Axis ratio [-]')
plt.title('Graupel')



ar_bins = np.arange(0.3,1.3,0.01)

ar_pdf  = np.zeros((len(ar_bins),len(bins_centers)))

for i,d in enumerate(bins_centers):
    sigm = ar_std_fit(d)
    mu = ar_fit(d)
    func = lambda ar: np.exp(-(ar-mu)**2/(2*sigm**2))/np.sqrt(2*np.pi*sigm**2)
    ar_pdf[:,i] = func(ar_bins)
    
x,y = np.meshgrid(bins_centers,ar_bins)
plt.contourf(x,y,ar_pdf,levels=np.linspace(np.nanmin(ar_pdf),np.nanmax(ar_pdf),50))
plt.colorbar(label='Density [-]')


