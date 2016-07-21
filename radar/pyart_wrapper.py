# -*- coding: utf-8 -*-
"""
Created on Wed Jan 27 10:41:15 2016

@author: wolfensb
"""

import cosmo_pol.constants.constants as constants

from kdp_estimation_kalman_v11 import estimate_kdp_kalman, estimate_kdp_vulpiani,filter_psidp

import read_radar_data
import pyart
import os, glob, datetime
import numpy as np
import re
import numpy.ma as ma

# Specify constants
# RADAR
RADAR_COORDS={}
RADAR_COORDS['Lema']=( 46.04079, 8.833219,1610.)
RADAR_COORDS['Dole']=(46.42562, 6.09997,1682.6)
RADAR_COORDS['Albis']=(47.28433, 8.51201,938.)
RADAR_COORDS['Plaine_Morte']=(46.37059, 7.48658,2930)

SENSITIVITY_FACTOR = -100
RADAR_CONSTANT = 75

# MODEL
UNITS_SIMUL={'ZH':'dBZ','KDP':'deg/km','PHIDP':'deg','RHOHV':'-','ZDR':'dB','RVEL':'m/s','DSPECTRUM':'dBZ','ZV':'dBZ',\
'U':'m/s','V':'m/s','W':'m/s','T':'K','RHO':'kg/m3','QR_v':'kg/m3','QS_v':'kg/m3',\
'QG_v':'kg/m3','QH_v':'kg/m3','std_ZH':r'mm^6 m^{-3}','std_KDP':'deg/km','std_ZDR':'-','std_RHOHV':'-','std_PHIDP':'deg',
'ZV':'dBZ'}

VAR_LABELS_SIMUL={'ZH':'Reflectivity','KDP':'Specific diff. phase','PHIDP':'Diff. phase',\
'RHOHV':'Copolar corr. coeff.','ZDR':'Diff. reflectivity','RVEL':'Mean doppler velocity','DSPECTRUM':'Doppler spectrum','ZV':'Vert. reflectivity',\
'U':'U-wind component','V':'V-wind component','W':'Vertical wind component','T':'Temperature','RHO':'Air density',\
'QR_v':'Mass density of rain','QS_v':'Mass density of snow','QG_v':'Mass density of graupel','QH_v':'Mass density of hail',
'std_ZH':'Stdev of hor. reflectivity','std_KDP':'Stdev of spec. diff. phase shift',
'std_ZDR':'Stdev of diff. reflectivity','std_RHOHV':'Stdev of copolor corr. coeff.',
'std_ZV':'Stdev of vert. reflectivity','PHIDP':'Diff. phase shift'}


VMIN_SIMUL={'ZH':0.,'KDP':0.,'PHIDP':0.,'RHOHV':0.6,'ZDR':0.,'RVEL':-25,'DSPECTRUM':-50,'ZV':0.,\
'U':-30,'V':-30,'W':-10,'T':200,'RHO':0.5,'QR_v':0.,'QS_v':0.,\
'QG_v':0.,'QH_v':0.}

VMAX_SIMUL={'ZH':55,'KDP':1,'PHIDP':20,'RHOHV':1,'ZDR':2,'RVEL':25,'DSPECTRUM':30
,'ZV':45,'U':30,'V':30,'W':10,'T':300,'RHO':1.4,'QR_v':1E-3,'QS_v':1E-3,\
'QG_v':1E-3,'QH_v':1E-2}


''' Inputs :
-------------------------------------------------------------------------------
  - filename : string:  Path of the radar PPI scan
  - high_res : boolean: True if high_res (83.3 m. radial res)
                        False if low_res (500 m. radial res)
  - vol_scan : boolean: True if all PPI scans for all elevations at that time-
                        step need to be loeaded
                        False if only that given scan (specified by filename)
                        needs to be loaded
  - max_range : float : maximum range from the radar to be considered (default
                        is infinity)
  - min_range : float : minimum range from the radar to be considered (default
                        is 10 km)
                        
'''

class RadarDisplay(pyart.graph.radardisplay.RadarDisplay):
    def __init__(self, radar, shift=(0.0, 0.0)):
        self.scan_type = None
        self._radar = None
        super(RadarDisplay,self).__init__(radar, shift=(0.0, 0.0))

    def plot(self, fields,sweep=0, max_range = 100000, **kwargs):
        filt = pyart.filters.GateFilter(self._radar)
        filt.exclude_above('range',max_range)
        
        if self.scan_type == 'vprof':
            self.plot_vprof(fields, **kwargs)
        else:
            super(RadarDisplay,self).plot(field = fields, sweep = sweep, gatefilter = filt, **kwargs)
    
    def plot_vprof(self, fields, sweep=0, mask_tuple=None,
            vmin=None, vmax=None, norm=None, cmap=None, mask_outside=False,
            title=None, title_flag=True,
            axislabels=(None, None), axislabels_flag=True,
            colorbar_flag=True, colorbar_label=None,
            colorbar_orient='vertical', edges=True, gatefilter=None,
            filter_transitions=True, ax=None, fig=None, **kwargs):
        
        if np.isscalar(fields):
            fields=[fields]
            
        ranges = self._radar.range['data']
            
        ax, fig = pyart.graph.common.parse_ax_fig(ax, fig)
        ax.hold(True)
        for field in fields:
            data = self._radar.get_field(sweep, field)
            if field == 'DSpectrum': # 1D data
                varray = constants.VARRAY
                pm = ax.pcolormesh(varray, ranges, data, vmin=vmin, vmax=vmax, cmap=cmap, norm=norm, **kwargs)
                self.set_limits(ylim=[0,10000],xlim=[np.min(varray),2])
                ax.set_xlabel('Velocity [m/s]')
                self.plot_colorbar(
                mappable=pm, label=colorbar_label, orient=colorbar_orient,
                field=field, ax=ax, fig=fig)
            else:
                if norm is None:
                    vmin, vmax = pyart.graph.common.parse_vmin_vmax(self._radar, field, vmin, vmax)
                else:
                    vmin = None
                    vmax = None
                if cmap is None:
                    cmap = pyart.config.get_field_colormap(field)
                    ax.plot(data,ranges,linewidth=2)
                self.set_limits(ylim=[0,10000])
                ax.set_xlabel(VAR_LABELS_SIMUL[field]+' ('+UNITS_SIMUL[field]+')')
        ax.grid()
        ax.set_ylabel('Height above radar [m]')
        
class PyradRadop(pyart.core.Radar):
    def __init__(self,scan_type, scan):
                    
        N_sweeps=len(scan['data'])
        
        fields={}
        fixed_angle={}
        fixed_angle['data']=np.zeros(N_sweeps,)
                
        sweep_start_ray_index={}
        sweep_start_ray_index['data']=[]
        sweep_stop_ray_index={}
        sweep_stop_ray_index['data']=[]
        
        varnames=scan['data'][0][0].values.keys()
    
        for i,k in enumerate(varnames):
            
            fields[k]={}
            fields[k]['data']=[]
            try: # No info found in header
                fields[k]['long_name']=VAR_LABELS_SIMUL[k]
            except:
                pass
            try:
                fields[k]['units']=UNITS_SIMUL[k]
            except:
                pass
            try:
                fields[k]['valid_min']=VMIN_SIMUL[k]
            except:
                pass
            try:
                fields[k]['valid_max']=VMAX_SIMUL[k]
            except: 
                pass
                
        # Initialize
        idx_start=0
        idx_stop=0
        elevations=[]
        azimuths=[]

        for i in range(N_sweeps):
            # Convert list of beams to array of data in polar coordinates
            polar_data_sweep={}
            for k in varnames:
                polar_data_sweep[k]=np.array([it.values[k] for it in scan['data'][i]])
                if k in ['ZDR','ZV','ZH','DSpectrum']: # Convert to dB
                    polar_data_sweep[k][polar_data_sweep[k]==0]=float('nan')
                    polar_data_sweep[k]=10*np.log10(polar_data_sweep[k])
                    
            [N_angles,N_ranges]=polar_data_sweep[varnames[0]].shape
            
            idx_stop=idx_start+N_angles-1
            sweep_start_ray_index['data'].append(idx_start)
            sweep_stop_ray_index['data'].append(idx_stop)
            idx_start=idx_stop+1
        
            if scan_type=='ppi':
                fixed_angle['data']=scan['elevations'][i]
                elevations.extend(list([scan['elevations'][i]]*N_angles))
                azimuths.extend(list(scan['azimuths']))
            elif scan_type=='rhi':
                fixed_angle['data']=scan['azimuths'][i]
                elevations.extend(list(scan['elevations']))
                azimuths.extend(list([scan['azimuths'][i]]*N_angles))
                
            for k in varnames:
                if not len(fields[k]['data']):
                    fields[k]['data']=polar_data_sweep[k]
                else:
                    fields[k]['data']=read_radar_data.row_stack(fields[k]['data'],polar_data_sweep[k])

        for k in varnames:
            fields[k]['data']=np.ma.array(fields[k]['data'],mask=np.isnan(fields[k]['data']))
        
        
        metadata={}
        
        
        # Position and time are obtained from the pos_time field of the scan dictionary
        latitude={'data' : np.array(scan['pos_time']['latitude'])}
        longitude={'data' :  np.array(scan['pos_time']['longitude'])}    
        altitude={'data' :  np.array(scan['pos_time']['altitude'])}   
        
        time_units='seconds since '+scan['pos_time']['time']
        time={'data' : np.zeros((N_angles,)),'units': time_units}
    
    
        sweep_number={'data' : np.arange(0,N_sweeps)}     
        sweep_mode={'data' : [scan_type]*N_sweeps}   

        metadata={}
        
        azimuth={'data' : np.array(azimuths)}
        rrange={'data': scan['ranges']}
        elevation={'data' :np.array(elevations)}

        # Finally add ranges as an additional variable, for convenience in order to
        # filter gates based on their range

        fields['range'] = {}
        fields['range']['data'] = np.tile(rrange['data'],(len(elevation['data']),1))
        
        # Create PyART instance
        super(PyradRadop,self).__init__(time,rrange,fields,metadata,
        scan_type,latitude,longitude,altitude,sweep_number,sweep_mode,
        fixed_angle,sweep_start_ray_index,sweep_stop_ray_index,azimuth, 
        elevation)
        
        
class PyradRadopVProf(PyradRadop):
    def __init__(self, scan):
        if not isinstance(scan['data'],list):
            scan['data']=[[scan['data']]]
        super(PyradRadopVProf,self).__init__('vprof',scan)
    
    def get_field(self,sweep_idx, variable):
        out = super(PyradRadopVProf,self).get_field(sweep_idx, variable)
        return np.squeeze(out)
        
class PyradMXPOL(pyart.core.Radar):
    def __init__(self,filename, max_range=np.Inf,min_range=10000):
        
        # Get directory name
        dirname=os.path.dirname(filename)+'/'
        if dirname=='/':
            dirname='./'
    
        all_files=[filename]
        
        
        fname_basename=os.path.basename(filename)

        # Specify loooots of inputs for the Pyart class
        if 'PPI' in fname_basename:
            scan_type='ppi'
        elif 'RHI' in fname_basename:
            scan_type='rhi'
            
        strdate=re.findall(r"([0-9]{8}-[0-9]{6})",fname_basename)[0] # Found out from the filename
        date=datetime.datetime.strptime(strdate,'%Y%m%d-%H%M%S')
        
        varnames=['Zh','Zdr','Kdp','Phidp','Rhohv','ZhCorr','ZdrCorr','RVel','Sw']
        labels=['Reflectivity','Diff. reflectivity','Spec. diff. phase','Diff. phase','Copolar corr. coeff','Att. corr reflectivity',\
        'Att corr. diff. reflectivity.','Mean doppler velocity','Spectral Width']
        
        units=['dBZ','dB','deg/km','deg','-','dBZ','dB','m/s','m2/s2']
        
    
        vmin=[0.,0.,0.,0.,0.6,0.,0.,-15.,0.]
        vmax=[55.,3.,4.,45.,1.,55.,3.,15.,3.]
                
        N_sweeps=len(all_files)
        fields={}
        fixed_angle={}
        fixed_angle['data']=np.zeros(N_sweeps,)
        
        sweep_start_ray_index={}
        sweep_start_ray_index['data']=[]
        sweep_stop_ray_index={}
        sweep_stop_ray_index['data']=[]
        
        for i,k in enumerate(varnames):
            fields[k]={}
            fields[k]['data']=[]
            fields[k]['long_name']=labels[i]
            fields[k]['units']=units[i]
            fields[k]['valid_min']=vmin[i]
            fields[k]['valid_max']=vmax[i]
            
        # Initialize
        idx_start=0
        idx_stop=0
        elevations=[]
        azimuths=[]
        nyquist=[]

        for i in range(N_sweeps):
      
            data=read_radar_data.readMXPOLRadData(all_files[i],varnames,max_range)
            if scan_type == 'rhi':
                fixed_angle['data']=data['elevation']
            elif scan_type == 'ppi':
                fixed_angle['data']=data['azimuth']
            

            [N_ranges,N_az]=data[varnames[0]].shape
            idx_stop=idx_start+N_az-1
            sweep_start_ray_index['data'].append(idx_start)
            sweep_stop_ray_index['data'].append(idx_stop)
            idx_start=idx_stop+1
            elevations.extend(list(data['elevation']))
            nyquist.extend([data['nyquist_vel']]*N_az)
            azimuths.extend(list(data['azimuth']))
            
            
            for j,v in enumerate(varnames):
                if v in data.keys():
                    if not len(fields[v]['data']):
                        fields[v]['data']=data[v]
                    else:
                        fields[v]['data']=read_radar_data.row_stack(fields[v]['data'],data[v])
                else:
                    print('Variable '+v+' was not found in file!')
        
        for v in varnames:
            fields[v]['data']=np.ma.array(fields[v]['data'],mask=fields[v]['data'] == -99900.0)

        metadata={}
        
        [a,N_ranges]=fields[varnames[0]]['data'].shape
        
        latitude={'data' : data['latitude']}
        longitude={'data' :data['longitude']}    
        altitude={'data' : data['altitude']}     
        sweep_number={'data' : np.arange(0,len(all_files))}     
        sweep_mode={'data' : [scan_type]*N_sweeps}   
        instrument_parameters={'nyquist_velocity': {'data':np.array(nyquist)}}
        
        metadata={}
        
        azimuth={'data' : np.array(azimuths)}
        rrange={'data':np.arange(N_ranges)*data['resolution']}
        elevation={'data' :np.array(elevations)}
    
        time_units='seconds since '+str(date)
        time={'data' : data['time'],'units': time_units}

        # Finally add ranges as an additional variable, for convenience in order to
        # filter gates based on their range

        fields['range'] = {}
        fields['range']['data'] = np.tile(rrange['data'],(len(elevation['data']),1))
        
        # Create PyART instance
        super(PyradMXPOL,self).__init__(time,rrange,fields,metadata,scan_type,latitude,longitude,altitude,sweep_number,sweep_mode,fixed_angle,\
        sweep_start_ray_index,sweep_stop_ray_index,azimuth, elevation,instrument_parameters=instrument_parameters)
    
    def estimate_KDP(self, method = 'vulpiani', minsize = 5, thresh_rhohv=0.65, max_discont = 90):
        band = 'C'
        dr = (self.range['data'][1] -  self.range['data'][0])/1000.

        kdp = {}
        kdp['data'] =  ma.masked_array(np.zeros(self.fields['Phidp']['data'].shape))*np.nan
        if method == 'kalman':
            stdev_kdp = {}
            stdev_kdp['data'] = ma.masked_array(np.zeros(self.fields['Phidp']['data'].shape))*np.nan
        phidp_filt = {}
        phidp_filt['data'] =ma.masked_array(np.zeros(self.fields['Phidp']['data'].shape))*np.nan
        
        # Filter Psidp
        idx_line = 0
        for n in self.sweep_number['data']:
            psidp_filt = filter_psidp(self.get_field(n,'Phidp'), self.get_field(n,'Rhohv'),
                         minsize, thresh_rhohv, max_discont)
            import matplotlib.pyplot as plt
            plt.imshow(psidp_filt)
            if method == 'vulpiani':
                a, b = estimate_kdp_vulpiani(psidp_filt, dr, windsize = 7, 
                                                band = band)
            elif method == 'kalman':
                a,c, b = estimate_kdp_kalman(psidp_filt, dr, band)
            dim = a.shape
            if method == 'kalman':
                stdev_kdp['data'][idx_line:idx_line+dim[0],0:dim[1]] = c
            kdp['data'][idx_line:idx_line+dim[0],0:dim[1]] = a
            phidp_filt['data'][idx_line:idx_line+dim[0],0:dim[1]] = b
            
            idx_line = idx_line + dim[0]
        
        kdp['units']='deg/km'
        kdp['valid_min']=np.nanmin(kdp['data'])
        kdp['valid_max']=np.nanmax(kdp['data'])    
        
        if method == 'kalman':
            stdev_kdp['units']='deg/km'
            stdev_kdp['valid_min']=np.nanmin(stdev_kdp['data'])
            stdev_kdp['valid_max']=np.nanmax(stdev_kdp['data'])    
            self.add_field('KDP_STD',stdev_kdp)
        phidp_filt['units']='deg'
        phidp_filt['valid_min']=np.nanmin(phidp_filt['data'])
        phidp_filt['valid_max']=np.nanmax(phidp_filt['data'])                
        
        self.add_field('Kdp',kdp)
        self.add_field('Phidp_filt',phidp_filt)
        
class PyradCH(pyart.core.Radar):
    def __init__(self,filename, high_res, vol_scan=False, max_range=np.Inf):
        
        # Get directory name
        dirname=os.path.dirname(filename)+'/'
        if dirname=='/':
            dirname='./'
    
        # If vol_scan : retrieve files with same timestamp
        fname_basename=os.path.basename(filename)
        if vol_scan:
            all_files=np.sort(glob.glob(dirname+fname_basename[0:15]+'*.h5'))
        else:
            all_files=[filename]
        
        # Get name of radar
        index_letter=fname_basename[2]
        if index_letter == 'A':
            radar_name='Albis'
        elif index_letter == 'L':
            radar_name='Lema'
        elif index_letter == 'D':
            radar_name='Dole'
        elif index_letter == 'P':
            radar_name='Plaine_Morte'
        
        # Get radar resolution
        if high_res:
            rres=83.3
        else:
            rres=500.
        
        # Specify loooots of inputs for the Pyart class
        
        scan_type='ppi'
        time=datetime.datetime.strptime(fname_basename[3:12],'%y%j%H%M')
        
        varnames=['Z','ZDR','ZV','V','W','RHO','CLUT','PHIDP']
        labels=['Reflectivity','Diff. reflectivity','Vert. reflectivity','Mean doppler velocity','Spectral Width','Copolar corr. coeff.','Clutter','Diff. phase']
        units=['dBZ','dB','dBZ','m/s','(m/s)^2','-','-','deg']
        vmin=[0,0,0,-15,0,0.6,0,0]
        vmax=[55.,3.,45.,15,3,1,100,150]
        
        N_sweeps=len(all_files)
        
        fields={}
        fixed_angle={}
        fixed_angle['data']=[]                
        
        sweep_start_ray_index={}
        sweep_start_ray_index['data']=[]
        sweep_stop_ray_index={}
        sweep_stop_ray_index['data']=[]
        
        for i,k in enumerate(varnames):
            fields[k]={}
            fields[k]['data']=[]
            fields[k]['long_name']=labels[i]
            fields[k]['units']=units[i]
            fields[k]['valid_min']=vmin[i]
            fields[k]['valid_max']=vmax[i]
            
        # Initialize
        idx_start=0
        idx_stop=0
        elevations=[]
        azimuths=[]
        nyquist=[]
        
        for i in range(N_sweeps):
            data=read_radar_data.readCHRadData(all_files[i],varnames,rres,max_range)
            fixed_angle['data'].append(data['elevation'])
            [N_az,N_ranges]=data[varnames[0]].shape
            idx_stop=idx_start+N_az-1
            sweep_start_ray_index['data'].append(idx_start)
            sweep_stop_ray_index['data'].append(idx_stop)
            idx_start=idx_stop+1
            elevations.extend([data['elevation']]*N_az)
            nyquist.extend([data['nyquist_vel']]*N_az)
            azimuths.extend(list(data['azimuth']))
            for j,v in enumerate(varnames):
                if not len(fields[v]['data']):
                    fields[v]['data']=data[v]
                else:
                    fields[v]['data']=read_radar_data.row_stack(fields[v]['data'],data[v])
                    
        for v in varnames:
            fields[v]['data']=np.ma.array(fields[v]['data'],mask=np.isnan(fields[v]['data']))
                
        metadata={}
        
        [a,N_ranges]=fields[varnames[0]]['data'].shape
        
        latitude={'data' : np.array([RADAR_COORDS[radar_name][0]])}
        longitude={'data' : np.array([RADAR_COORDS[radar_name][1]])}    
        altitude={'data' : np.array([RADAR_COORDS[radar_name][2]])}     
        sweep_number={'data' : np.arange(0,len(all_files))}     
        sweep_mode={'data' : ['ppi']*N_sweeps}   
        instrument_parameters={'nyquist_velocity': {'data':np.array(nyquist)}}
        
        metadata={}
        
        azimuth={'data' : [np.array(azimuths)]*N_sweeps}
        rrange={'data':np.arange(N_ranges)*data['resolution']}
        elevation={'data' :[np.array(elevations)]*N_sweeps}
        time_units='seconds since '+str(time)
        time={'data' : np.zeros((len(elevations),)),'units': time_units}
    
        # Finally add ranges as an additional variable, for convenience in order to
        # filter gates based on their range

        fields['range'] = {}
        fields['range']['data'] = np.tile(rrange['data'],(len(elevation['data']),1))
        fields['range']['data'] = ma.masked_array(fields['range']['data'],
                                    mask = np.isnan(fields['range']['data']))
        # Create PyART instance
        super(PyradCH,self).__init__(time,rrange,fields,metadata,scan_type,latitude,longitude,altitude,sweep_number,sweep_mode,fixed_angle,\
        sweep_start_ray_index,sweep_stop_ray_index,azimuth, elevation,instrument_parameters=instrument_parameters)
        
    def average(self,n_gates = 6):
        
        n_bins_initial = self.range['data'].shape[0]
        
        for k in self.fields.keys():
            if k == 'range':
                data = self.fields[k]['data']
                self.fields[k]['data'] = np.vstack([data \
                [:,n_gates*(i)] for i in  \
                range(n_bins_initial/n_gates)]).T
            elif k == 'Z': # Convert to linear
                data = 10*np.log10(self.fields[k]['data']) 
                self.fields[k]['data'] = 10**(0.1*(np.vstack(\
                [np.nanmean(data[:,n_gates*i:n_gates*(i+1)],axis=1) \
                for i in range(n_bins_initial/n_gates)]).T))
            elif k == 'ZDR':
                data = np.log10(self.fields[k]['data'])
                self.fields[k]['data'] = 10**(np.vstack( \
                [np.nanmean(data[:,n_gates*i:n_gates*(i+1)],axis=1) \
                for i in range(n_bins_initial/n_gates)]).T)
            elif k == 'PHIDP':
                data = self.fields[k]['data']
                self.fields[k]['data'] = np.vstack([data \
                [:,n_gates*(i)] for i in  \
                range(n_bins_initial/n_gates)]).T
            else:
                data = self.fields[k]['data']
                self.fields[k]['data'] = np.vstack([np.nanmean(data \
                [:,n_gates*i:n_gates*(i+1)],axis=1) for i in  \
                range(n_bins_initial/n_gates)]).T  

            
        self.range['data'] = np.array([self.range['data'] \
                [n_gates*(i)] for i in  \
                range(n_bins_initial/n_gates)])
        self.ngates = len(self.range['data'])

    def estimate_KDP(self, method = 'vulpiani', minsize = 5, thresh_rhohv=0.65, max_discont = 90):
        band = 'C'
        dr = (self.range['data'][1] -  self.range['data'][0])/1000.

        kdp = {}
        kdp['data'] =  ma.masked_array(np.zeros(self.fields['PHIDP']['data'].shape))*np.nan
        if method == 'kalman':
            stdev_kdp = {}
            stdev_kdp['data'] = ma.masked_array(np.zeros(self.fields['PHIDP']['data'].shape))*np.nan
        phidp_filt = {}
        phidp_filt['data'] =ma.masked_array(np.zeros(self.fields['PHIDP']['data'].shape))*np.nan
        
        # Filter Psidp
        idx_line = 0
        for n in self.sweep_number['data']:
            psidp_filt = filter_psidp(self.get_field(n,'PHIDP'), self.get_field(n,'RHO'),
                         minsize, thresh_rhohv, max_discont)
            if method == 'vulpiani':
                a, b = estimate_kdp_vulpiani(psidp_filt, dr, windsize = 7, 
                                                band = band)
            elif method == 'kalman':
                a,c, b = estimate_kdp_kalman(psidp_filt, dr, band)
            dim = a.shape
            if method == 'kalman':
                stdev_kdp['data'][idx_line:idx_line+dim[0],0:dim[1]] = c
            kdp['data'][idx_line:idx_line+dim[0],0:dim[1]] = a
            phidp_filt['data'][idx_line:idx_line+dim[0],0:dim[1]] = b
            
            idx_line = idx_line + dim[0]
        
        kdp['units']='deg/km'
        kdp['valid_min']=np.nanmin(kdp['data'])
        kdp['valid_max']=np.nanmax(kdp['data'])    
        
        if method == 'kalman':
            stdev_kdp['units']='deg/km'
            stdev_kdp['valid_min']=np.nanmin(stdev_kdp['data'])
            stdev_kdp['valid_max']=np.nanmax(stdev_kdp['data'])    
            self.add_field('KDP_STD',stdev_kdp)
        phidp_filt['units']='deg'
        phidp_filt['valid_min']=np.nanmin(phidp_filt['data'])
        phidp_filt['valid_max']=np.nanmax(phidp_filt['data'])                
        
        self.add_field('KDP',kdp)
        self.add_field('PHIDP_FILT',phidp_filt)
    
    def snr_threshold(self,threshold = 8):
        threshold_func = lambda r: SENSITIVITY_FACTOR + RADAR_CONSTANT + threshold + 20*np.log10(r)
        range_mat = np.tile(self.range['data'],(self.fields['Z']['data'].shape[0],1))/1000.
        mask = np.array(self.fields['Z']['data'])<threshold_func(range_mat)
        for k in self.fields.keys():
            self.fields[k]['data'].mask = mask

    def correct_velocity(self):
        corr_vel=pyart.correct.dealias_region_based(self,interval_splits=3,vel_field='V',rays_wrap_around=True)
        corr_vel['units']='m/s'
        corr_vel.pop('standard_name')
        corr_vel['valid_min']=np.nanmin(corr_vel['data'])
        corr_vel['valid_max']=np.nanmax(corr_vel['data'])    
        self.add_field('V_corr',corr_vel)
        
    @staticmethod
    def ang_to_sweep(angle):
        from read_radar_data import dic_elev, ELEVATION_ANGLES
        try:
            sweep=dic_elev[str(angle)]
        except:
            print('No PPI scan available for that elevation angle')
            print('Specify one of these')
            print(str(ELEVATION_ANGLES))
            return
            
        return sweep

if __name__=='__main__':
#    a=pyrad_RADOP('ppi',r)
#    f='/ltedata/HYMEX/SOP_2012/Radar/Proc_data/2012/09/29/MXPol-polar-20120929-123446-RHI-233_6.nc'
#    a=PyradMXPOL(f)
#    a.estimate_KDP()
    
#    
#    display = pyart.graph.RadarDisplay(a, shift=(0.0, 0.0))
#    display.plot('v_radial', 0, colorbar_flag=True,title="v_radial")
#    
    plt.close('all')
    filename='/ltedata/MeteoSwiss_Full_Radar_Data_LowRes/PLD14098/PLD1409805007U.005.h5'
    r=PyradCH('/ltedata/MeteoSwiss_Full_Radar_Data_LowRes/PLD14098/PLD1409805007U.005.h5',high_res=False,vol_scan=False)
    r.snr_threshold(8)
#    print(r.range['data'])
#    r.average(1)
    r.estimate_KDP(method='vulpiani')
#    print(r.range['data'])

    
#    rad_instance.correct_velocity()
#    plt.figure()
    display = pyart.graph.RadarDisplay(rad, shift=(0.0, 0.0))
    plt.figure()
    display.plot('KDP', 0, colorbar_flag=True,title="V_corr",vmin=-1,vmax=1.5) 

#    display.plot('PHIDP_FILT', 0, colorbar_flag=True,title="V_corr",vmin=-1) 
#    plt.figure()
#    display.plot('PHIDP', 0, colorbar_flag=True,title="V_corr",vmin=-1)
#        
#    plt.figure()
#    display.plot('Z', 0, colorbar_flag=True,title="V_corr",vmin=-1)        
#    import pickle
    
#    a = pickle.load(open('ex_rad.p','rb'))
    
#    a=pickle.load(open('../ex_vprof.txt'))
#    
#    b=PyradRadopVProf(a)
#    field=b.get_field(0,'DSpectrum')
#    
#    display = RadarDisplay(b, shift=(0.0, 0.0))
#    display.plot(['RVel','DSpectrum'], 0, colorbar_flag=True,title="ZH")