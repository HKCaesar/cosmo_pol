
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 11 09:57:28 2016

@author: wolfensb
"""
# Setup DEFAULTS
global CONFIG


import numpy as np
import numbers
import yaml
from cosmo_pol.utilities.antenna_fit import optimize_gaussians

FLAG_NUMBER = 'ANY_NUMBER'


CONFIG = {}

DEFAULTS={
    'radar':
        {'coords':[46.42563,6.01,1672.7],\
        'type':'ground',\
        'frequency':5.6,\
        'range':150000,\
        'radial_resolution':500,\
        'PRI':700,\
        'FFT_length':256,\
        'sensitivity':[-5,10000],\
        '3dB_beamwidth':1.},\
    'attenuation':
        {'correction':1},\
    'refraction':
        {'scheme':1},\
    'integration':
        {'scheme':1,\
        'nv_GH':9,\
        'nh_GH':3,\
        'n_gaussians':7,\
        'antenna_diagram':'',\
        'nr_GH':7,\
        'na_GL':7},\
    'doppler':
        {'scheme':1,\
        'turbulence_correction':0},\
    'microphysics':
        {'scheme':1},\
    }

# Setup range of valid values
    
VALID_VALUES={
    'radar':
        {'coords':FLAG_NUMBER,\
        'type':['ground','GPM'],\
        'frequency':[2,'to',35.6],\
        'range':[5000,'to',500000],\
        'radial_resolution':[25,'to',5000],\
        'PRI':[10,'to',3000],\
        'FFT':[16,'to',2048],\
        'sensitivity':FLAG_NUMBER,\
        '3dB_beamwidth':[0.1,'to',10]},\
    'attenuation':
        {'correction':[0,1]},\
    'refraction':
        {'scheme':[1,2]},\
    'integration':
        {'scheme':[1,2],\
        'nv_GH':np.arange(1,31,2),\
        'nh_GH':np.arange(1,31,2),\
        'n_gaussians':np.arange(1,13,2),\
        'nr_GH':np.arange(1,31,2),\
        'na_GL':np.arange(1,31,2)},\
    'doppler':
        {'scheme':[1,2,3],\
        'turbulence_correction':[0,1]},\
    'microphysics':
        {'scheme':[1,2]},\
    }

def check_valid(section,key, value):

    flag_ok = False
    out = value
    message = ''
    
    if section in VALID_VALUES.keys() and key in VALID_VALUES[section].keys():
        valid_vals=VALID_VALUES[section][key]
        if isinstance(valid_vals, list):
            if 'to' in valid_vals:
                if value>=float(valid_vals[0]) and value<=float(valid_vals[2]):
                    flag_ok = True
                    message = ''
                else:
                    flag_ok = False
                    message = 'Invalid value for the '+section+': '+key+' parameter'+\
                    'Please choose one of the following values: '+\
                    'from '+str(valid_vals[0])+' to '+str(valid_vals[2])+\
                    'Using default option: '+section+': '+key+' = '+str(DEFAULTS[section][key]+'\n')
                    out = DEFAULTS[section][key]
            elif value in valid_vals:
                flag_ok = True
                message = ''
            else:
                flag_ok = False
                message = 'Invalid value for the "'+section+': '+key+'" parameter \n'+\
                    'Please choose one of the following values: '+\
                    '['+', '.join([str(s) for s in valid_vals])+'] \n'\
                    'Using default option: "'+section+': '+key+'" = '+str(DEFAULTS[section][key]+'\n')
        elif valid_vals == FLAG_NUMBER:
            if not hasattr(value, '__len__'):
                value = [value]
            for v in value:
                if not isinstance(v,numbers.Number):
                    flag_ok = False
                    message = 'Invalid value for the "'+section+': '+key+'" parameter \n'+\
                    'All values must be numbers \n' \
                    'Using default option: "'+section+': '+key+'" = '+str(DEFAULTS[section][key])+'\n'
                    out = DEFAULTS[section][key]
    else:
        flag_ok = True
        message='Parameter "'+section+': '+key+'" was not tested for valid range.\n'+\
        'It is your responsability to give relevant values \n'
    return flag_ok, out, message

def init(options_file):
    global CONFIG
    try:
        with open(options_file, 'r') as ymlfile:
            CONFIG = yaml.load(ymlfile)
    except:
        CONFIG = DEFAULTS
        print('Could not find '+options_file+' file, using default options')
        return
    # Parsing values
    for section in DEFAULTS:
        if section not in CONFIG.keys():
            CONFIG[section] = {}
        for key in DEFAULTS[section]:
            if key not in CONFIG[section].keys():
                CONFIG[section][key] = DEFAULTS[section][key]
            
            flag_ok, value, message = check_valid(section,key,CONFIG[section][key])
           
            if message != '':
                print(message)

            CONFIG[section][key] = value

        if section == 'integration': # Treat special case
            if CONFIG[section]['antenna_diagram'] != '' :

                try:
                    print('Trying to fit sum of gaussians on the provided antenna diagram...\n')
                    data = np.genfromtxt(CONFIG[section]['antenna_diagram'],delimiter=',')
                    x = data[:,0]
                    y = data[:,1]
                    gauss_params = optimize_gaussians(x,y,CONFIG[section]['n_gaussians'])
                    CONFIG[section]['antenna_params'] = gauss_params
                    print('Fit was successful !\n')
                except:
                    raise
                    print('Could not fit a sum of gaussians with the file provided in "\n'+\
                    '"'+section+':'+' antenna_diagram'+'"\n'+\
                    'Please proviide a comma-separated file with two columns'+\
                    ', one for the angles, one for the power in dB \n')       
    return

if __name__ == '__main__':
    init('/media/wolfensb/Storage/cosmo_pol/option_files/MXPOL_PPI.yml')  
    
    
    #print c
    #import yaml
    #
    #with open("test_yaml.yaml", 'r') as ymlfile:
    #    CONFIG = yaml.load(ymlfile)
    #
    #for section in CONFIG:
    #    print(section)
    #print(CONFIG['radar'])
    #print(CONFIG['doppler'])

