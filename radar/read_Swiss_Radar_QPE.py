# -*- coding: utf-8 -*-
"""
Created on Mon Jul 13 16:30:35 2015

@author: wolfensb
"""
import mahotas
import numpy as np
import matplotlib.pyplot as plt
import netCDF4


RADAR_POS=[(497100.4, 142465.1),(603688.7, 135469.1),(681202.8, 237606.0),(707957.0, 99764.5)] # Dole, Plaine Morte, Albis, Monte Lema

def CH1903toWGS84(x,y):
    # First, conversion from the military to the civillian system:
    y_p = (y - 600000)/1000000;
    x_p = (x - 200000)/1000000;

    # Step 2
    lambda_p = (2.6779094 + (4.728982 * y_p) + \
                + (0.791484 * y_p * x_p) + \
                + (0.1306 * y_p * pow(x_p, 2))) + \
                - (0.0436 * pow(y_p, 3))
    phi_p = (16.9023892 + (3.238272 * x_p)) + \
                - (0.270978 * pow(y_p, 2)) + \
                - (0.002528 * pow(x_p, 2)) + \
                - (0.0447 * pow(y_p, 2) * x_p) + \
                - (0.0140 * pow(x_p, 3))
    # Step 3 (conversion to degrees)
    lat = phi_p*100/36;
    lon = lambda_p*100/36;

    return lat,lon


def read_8bit_meteoswiss( filename, max_radar_radius ):
      
    data = mahotas.imread(filename)
    
    # Extract only the horizontal information:
    data = data[-640:, 0:710];
    
     # The co-ordinates in CH1903 form:
    y = np.linspace(255, 965, 710)*1000;
    x = np.linspace(480, -160,  640)*1000;
    
    #Construct the lookup table to map 8 bit integer value to rainfall rate in
    #mm/hr
    rr_lookup_table = np.zeros((256,))*float('nan');
    rr_lookup_table[1] = 0
    rr_lookup_table[2:21] = np.arange(0.1,1.05,0.05)
    rr_lookup_table[21] = 1.1
    rr_lookup_table[22:30] = np.arange(1.25,2.05,0.1)
    rr_lookup_table[30] = 2
    rr_lookup_table[31:168] = np.arange(3.05,140.05,1)
    rr_lookup_table[168:237] = np.arange(141.05,210.05,1)
    rr_lookup_table[237:251] = np.arange(210.05,350.05,10);
    
    rainfall_rates = rr_lookup_table[data];   #note that this means that 255 is incorrectly mapped, but still maps to a NaN, so not relevant

    # convert to WGS co-ordinates:
    Y,X=np.meshgrid(y,x)

    grid_scan=np.zeros((rainfall_rates.shape)) # Map which gives the number of radar by which every pixel is seen (for a given max radius from the radars)   
    for radar in RADAR_POS:
        dist=np.sqrt((Y-radar[0])**2+(X-radar[1])**2)
        grid_scan[dist<max_radar_radius]+=1
    rainfall_rates[grid_scan<1]=float('nan') # Remove pixels which are not scanned by any radar (for a given max radius from the radars)
    
    [ lat, lon ] = CH1903toWGS84( X, Y );    

    return lat, lon, rainfall_rates
    
def createNetCDF(lat,lon,rainfall_rates,filename):
  rootgrp=netCDF4.Dataset(filename, 'w', format='NETCDF3_CLASSIC') # Create handle for GPM netcdf

  [N,M]=lon.shape
  x_h = rootgrp.createDimension('x',N)
  y_h = rootgrp.createDimension('y',M)
  x_h = rootgrp.createVariable('lon2d','f4',('x','y',))
  y_h = rootgrp.createVariable('lat2d','f4',('x','y',))

  handle_var = rootgrp.createVariable('PRECIP','f4',('x','y',))
  #

  x_h[:] = lon
  y_h[:] = lat
  handle_var[:] = rainfall_rates
  rootgrp.close()

  return
  
