# -*- coding: utf-8 -*-
"""
Created on Wed Aug 17 12:23:24 2016

@author: wolfensb
"""

import numpy as np
import copy
import matplotlib.pyplot as plt
f = open('power_laws_cant','w')

plt.close('all')
agg = np.genfromtxt('stats_aggregatess_Davos_winter_2015-16.txt')
agg_sym = copy.deepcopy(agg[:,3])
agg_sym[agg_sym<0] = - agg_sym[agg_sym<0]
grau = np.genfromtxt('stats_graupels_Davos_winter_2015-16.txt')
grau_sym = copy.deepcopy(grau[:,3])
grau_sym[grau_sym<0] = - grau_sym[grau_sym<0]

agg_d_bins = np.arange(2.4,15.2,0.1)

canting = np.zeros((len(agg_d_bins)-1,))
canting_nosym = np.zeros((len(agg_d_bins)-1,))
canting_std = np.zeros((len(agg_d_bins)-1,))

for i in range(len(canting)):
    canting[i] = np.nanmean(agg_sym[np.logical_and(agg[:,1]<agg_d_bins[i+1],agg[:,1]>agg_d_bins[i])])
    canting_nosym[i] = np.nanmean(agg[:,3][np.logical_and(agg[:,1]<agg_d_bins[i+1],agg[:,1]>agg_d_bins[i])])
    canting_std[i] = np.nanstd(agg_sym[np.logical_and(agg[:,1]<agg_d_bins[i+1],agg[:,1]>agg_d_bins[i])])


bins_centers = (agg_d_bins[0:-1]+ agg_d_bins[1:])/2.

fit_mean = np.polyfit(np.log(bins_centers),np.log(canting),1)
fit_std = np.polyfit(np.log(bins_centers),np.log(canting_std),1)

f.write('Aggregates: std: a = '+str(np.exp(fit_std[1]))+', b = '+str(fit_std[0])+' \n') # python will convert \n to os.linesep


plt.plot(bins_centers,canting_nosym,'--b',linewidth=1.5)
plt.hold(True)
cant_fit = lambda D: np.exp(fit_mean[1])*D**fit_mean[0]
cant_std_fit = lambda D: np.exp(fit_std[1])*D**fit_std[0]

plt.xlabel('Agg. diameter [mm]')
plt.ylabel(r'Canting angle [$^{\circ}$]')
plt.title('Aggregates')

plt.figure()
plt.plot(bins_centers,canting_std,'--b',bins_centers,cant_std_fit(bins_centers), linewidth=1.5)




grau_d_bins = np.arange(0.8,5.,0.05)


canting = np.zeros((len(grau_d_bins)-1,))
canting_nosym = np.zeros((len(grau_d_bins)-1,))
canting_std = np.zeros((len(grau_d_bins)-1,))

for i in range(len(canting)):
    canting[i] = np.nanmean(grau_sym[np.logical_and(grau[:,1]<grau_d_bins[i+1],grau[:,1]>grau_d_bins[i])])
    canting_nosym[i] = np.nanmean(grau[:,3][np.logical_and(grau[:,1]<grau_d_bins[i+1],grau[:,1]>grau_d_bins[i])])
    canting_std[i] = np.nanstd(grau_sym[np.logical_and(grau[:,1]<grau_d_bins[i+1],grau[:,1]>grau_d_bins[i])])


bins_centers = (grau_d_bins[0:-1]+ grau_d_bins[1:])/2.

fit_mean = np.polyfit(np.log(bins_centers),np.log(canting),1)
fit_std = np.polyfit(np.log(bins_centers),np.log(canting_std),1)

f.write('-------\n')
f.write('Graupel: std: a = '+str(np.exp(fit_std[1]))+', b = '+str(fit_std[0])+' \n') # python will convert \n to os.linesep

plt.figure()
plt.plot(bins_centers,canting_nosym,'--b',linewidth=1.5)
plt.hold(True)
cant_fit = lambda D: np.exp(fit_mean[1])*D**fit_mean[0]
cant_std_fit = lambda D: np.exp(fit_std[1])*D**fit_std[0]

plt.xlabel('Graupel diameter [mm]')
plt.ylabel(r'Canting angle [$^{\circ}$]')
plt.title('Graupel')

plt.figure()
plt.plot(bins_centers,canting_std,'--b',bins_centers,cant_std_fit(bins_centers), linewidth=1.5)

f.close()




