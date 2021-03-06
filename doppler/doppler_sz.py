# -*- coding: utf-8 -*-
"""
Created on Wed Aug 12 15:34:56 2015

@author: wolfensb

TODO : CORRECTION FOR AIR DENSITY RHO !!!
"""

import numpy as np
from scipy.ndimage.filters import gaussian_filter
import copy

import doppler_c
from cosmo_pol.hydrometeors import hydrometeors
from cosmo_pol.utilities.beam import Beam
from cosmo_pol.utilities import cfg
from cosmo_pol.utilities import utilities
from cosmo_pol.constants import constants

DEG2RAD=np.pi/180.0

def proj_vel(U,V,W,VH,theta,phi):
    return (U*np.sin(phi)+V*np.cos(phi))*np.cos(theta)+(W-VH)*np.sin(theta)

def get_doppler_velocity(list_beams, lut_sz = 0):
    ###########################################################################
    # Get setup
    global doppler_scheme
    global microphysics_scheme
    
    doppler_scheme=cfg.CONFIG['doppler']['scheme']
    microphysics_scheme=cfg.CONFIG['microphysics']['scheme']

    # Check if necessary variables for turbulence scheme are present
    if doppler_scheme == 4 and 'EDR' in list_beams.keys(): 
        add_turb=True
    else:
        add_turb=False
        if doppler_scheme == 4:
             print('No eddy dissipitation rate variable found in COSMO file, could not use doppler_scheme == 4, using doppler_scheme == 3 instead')
    
  # Get dimensions
    num_beams=len(list_beams) # Number of beams
    idx_0=int(num_beams/2) # Index of central beam
    len_beams=max([len(l.dist_profile) for l in list_beams]) # Beam length
    
    
    if microphysics_scheme == 1:
        hydrom_types=['R','S','G'] # Rain, snow and graupel
    elif microphysics_scheme == 2:
        hydrom_types=['R','S','G','H'] # Add hail 
            
    # Create dic of hydrometeors
    list_hydrom={}
    for h in hydrom_types:
        list_hydrom[h]=hydrometeors.create_hydrometeor(h,microphysics_scheme)
    
    
    ###########################################################################    
    # Get radial wind and doppler spectrum (if scheme == 1 or 2)
    if doppler_scheme == 1 or doppler_scheme ==2:
        
        rvel_avg=np.zeros(len_beams,)*float('nan') # average radial velocity 
        sum_weights=np.zeros(len_beams,) # mask of GH weights
        
        for beam in list_beams:
            if doppler_scheme == 1: # Weighting by PSD only
                v_hydro=get_v_hydro_unweighted(beam,list_hydrom)
                
            elif doppler_scheme == 2: # Weighting by RCS and PSD
                v_hydro=get_v_hydro_weighted(beam,list_hydrom,lut_sz)

            # Get radial velocity knowing hydrometeor fall speed and U,V,W from model
            theta = beam.elev_profile*DEG2RAD
            phi = beam.GH_pt[0]*DEG2RAD

            proj_wind=proj_vel(beam.values['U'],beam.values['V'],beam.values['W'],v_hydro,theta,phi)

            # Get mask of valid values
            sum_weights=utilities.sum_arr(sum_weights,~np.isnan(proj_wind)*beam.GH_weight)

            # Average radial velocity for all sub-beams
            rvel_avg=utilities.nansum_arr(rvel_avg,(proj_wind)*beam.GH_weight)
            
        # We need to divide by the total weights of valid beams at every bin
        rvel_avg/=sum_weights
        
    elif doppler_scheme == 3:
        rvel_avg=np.zeros(len_beams,)
        doppler_spectrum=np.zeros((len_beams,len(constants.VARRAY)))
        for beam in list_beams:
            beam_spectrum=get_doppler_spectrum(beam,list_hydrom,lut_sz)*beam.GH_weight # Multiply by GH weight
            if add_turb: # Spectrum spread caused by turbulence
                turb_std=get_turb_std(constants.RANGE_RADAR,beam.values['EDR'])
                beam_spectrum=turb_spectrum_spread(beam_spectrum,turb_std)
            doppler_spectrum+=beam_spectrum
        try:
            rvel_avg=np.sum(np.tile(constants.VARRAY,(len_beams,1))*doppler_spectrum,axis=1)/np.sum(doppler_spectrum,axis=1)
        except:
            rvel_avg*=float('nan')
            

    ###########################################################################    
    # Get mask
    # This mask serves to tell if the measured point is ok, or below topo or above COSMO domain
    mask=np.zeros(len_beams,)

    for i,beam in enumerate(list_beams):
        mask=utilities.sum_arr(mask,beam.mask) # Get mask of every Beam
    mask/=num_beams # Larger than 1 means that every Beam is below TOPO, smaller than 0 that at least one Beam is above COSMO domain
    mask[np.logical_and(mask>=0,mask<1)]=0
    
    # Finally get vectors of distances, height and lat/lon at the central beam
    idx_0=int(len(list_beams)/2)
    heights_radar=list_beams[idx_0].heights_profile
    distances_radar=list_beams[idx_0].dist_profile
    lats=list_beams[idx_0].lats_profile
    lons=list_beams[idx_0].lons_profile

    if doppler_scheme == 3:
        dic_vars={'RVEL':rvel_avg,'DSPECTRUM':doppler_spectrum}
    else:
        # No doppler spectrum is computed
        dic_vars={'RVEL':rvel_avg}
        
    beam_doppler=Beam(dic_vars,mask,lats, lons, distances_radar, heights_radar)    
    return beam_doppler
        
def get_v_hydro_unweighted(beam, list_hydrom):
    vh_avg = np.zeros(beam.values['T'].shape)
    vh_avg.fill(np.nan)
    n_avg = np.zeros(beam.values['T'].shape)
    n_avg.fill(np.nan)
    
    for i,h in enumerate(list_hydrom.keys()):
        if cfg.CONFIG['microphysics']['scheme'] == 1:
            if h=='S':
                list_hydrom['S'].set_psd(beam.values['T'],beam.values['Q'+h+'_v'])
            else:
                list_hydrom[h].set_psd(beam.values['Q'+h+'_v'])
        elif cfg.CONFIG['microphysics']['scheme'] == 2:
            list_hydrom[h].set_psd(beam.values['QN'+h+'_v'],beam.values['Q'+h+'_v'])
        
        # Get fall speed
        vh,n=list_hydrom[h].integrate_V()

        vh_avg = utilities.nansum_arr(vh_avg,vh)
        n_avg = utilities.nansum_arr(n_avg,n)
    
    v_hydro_unweighted = vh_avg/n_avg # Average over all hydrometeors

    return v_hydro_unweighted*(beam.values['RHO']/beam.values['RHO'][0])**(0.5)
    
def get_v_hydro_weighted(beam, list_hydrom, lut_sz):
    hydrom_scheme = cfg.CONFIG['microphysics']['scheme']
    
    vh_avg = np.zeros(beam.values['T'].shape)
    vh_avg.fill(np.nan)
    n_avg = np.zeros(beam.values['T'].shape)
    n_avg.fill(np.nan)
    
    for i,h in enumerate(list_hydrom.keys()):
        # Get list of diameters for this hydrometeor
        list_D = lut_sz[h].axes[lut_sz[h].axes_names['d']]
        
        valid_data = beam.values['Q'+h+'_v'] > 0
        
        # Get all elevations
        elev = beam.elev_profile
        # Since lookup tables are defined for angles >0, we have to check
        # if angles are larger than 90°, in that case we take 180-elevation
        # by symmetricity
        elev_lut = copy.deepcopy(elev)
        elev_lut[elev_lut>90]=180-elev_lut[elev_lut>90]
        # Also check if angles are smaller than 0, in that case, flip sign
        elev_lut[elev_lut<0]=-elev_lut[elev_lut<0]         

        T=beam.values['T']     

        # Get SZ matrix
        sz = lut_sz[h].lookup_line(e = elev[valid_data],
                                       t = T[valid_data])
              
        ''' 
        Part 1: Query of the SZ Lookup table  and RCS computation
        '''
        # Get SZ matrix
        sz = lut_sz[h].lookup_line(e = elev_lut[valid_data],
                                   t = T[valid_data])
        # get RCS
        rcs = 2*np.pi*(sz[:,:,0] - sz[:,:,1] - sz[:,:,2] + sz[:,:,3])
        rcs = rcs.T
        
        ''' 
        Part 2 : Get the PSD of the particles
        '''
        QM = beam.values['Q'+h+'_v'] # Get mass densities
        # 1 Moment case
        if hydrom_scheme == 1:
            if h != 'S' :
                list_hydrom[h].set_psd(QM[valid_data])
            else:
                # For snow N0 is Q and temperature dependent
                list_hydrom[h].set_psd(T[valid_data],QM[valid_data])
        # 2 Moment case     
        elif hydrom_scheme == 2:
            QN = beam.values['QN'+h+'_v']  # Get concentrations as well
            list_hydrom[h].set_psd(QN[valid_data],QM[valid_data])   
                      
        N=list_hydrom[h].get_N(list_D)
        if len(N.shape) == 1:
            N = np.reshape(N,[len(N),1]) # To be consistent with the einsum dimensions
        ''' 
        Part 3 : Integrate
        '''
        # Get fall speed
        v_f = list_hydrom[h].get_V(list_D)
        vh_w = np.trapz(N*v_f[:,np.newaxis]*rcs,axis=0)
        n_w = np.trapz(N*rcs, axis=0)
        

        vh_avg[valid_data] = utilities.nansum_arr(vh_avg[valid_data],
                                                  vh_w)
        n_avg[valid_data] = utilities.nansum_arr(n_avg[valid_data],
                                                  n_w)

        
    v_hydro_weighted = vh_avg/n_avg # Average over all hydrometeors

    return v_hydro_weighted*(beam.values['RHO']/beam.values['RHO'][0])**(0.5)

def get_doppler_spectrum(beam, list_hydrom, lut_sz):

    # Get dimensions
    len_beam = len(beam.dist_profile) # Length of the considered bea
    hydrom_types = list_hydrom.keys()
    
    # Initialize matrix of reflectivities
    refl=np.zeros((len_beam, len(constants.VARRAY)),dtype='float32')
     
    # Create matrix with every column being the set of diameters used for every
    # hydrometeor
    D=np.column_stack((lut_sz[h].axes[lut_sz[h].axes_names['d']] \
                                                    for h in hydrom_types))
    D_min = [list_hydrom[h].d_min for h in hydrom_types]
    
    step_D = D[1,:]-D[0,:] # The diameter step for every hydrometeor class
    
    # Get all elevations
    elev = beam.elev_profile
    # Since lookup tables are defined for angles >0, we have to check
    # if angles are larger than 90°, in that case we take 180-elevation
    # by symmetricity
    elev_lut = copy.deepcopy(elev)
    elev_lut[elev_lut>90]=180-elev_lut[elev_lut>90]
    # Also check if angles are smaller than 0, in that case, flip sign
    elev_lut[elev_lut<0]=-elev_lut[elev_lut<0]       
    
    # Get azimuth angle
    phi = beam.GH_pt[0] 
    
    # Correction of velocity for air density
    rho_corr = (beam.values['RHO']/beam.values['RHO'][0])**0.5 
    
    # Initialize matrix of radar cross sections for all diam. and hydrometeors
    rcs = np.zeros(D.shape,dtype='float32')
    rcs.fill(np.nan)
    
    for i in range(len_beam):  # Loop on all radar gates
        if beam.mask[i] == 0:
            # Get parameters of the PSD (lambda or lambda/N0)
            for j,h in enumerate(hydrom_types):
                if cfg.CONFIG['microphysics']['scheme'] == 1:
                    if h=='S':
                        list_hydrom['S'].set_psd(beam.values['T'][i],beam.values['Q'+h+'_v'][i])
                    else:
                        list_hydrom[h].set_psd(beam.values['Q'+h+'_v'][i])
                elif cfg.CONFIG['microphysics']['scheme'] == 2:
                    list_hydrom[h].set_psd(beam.values['QN'+h+'_v'][i],beam.values['Q'+h+'_v'][i])
    
            N0=np.array([list_hydrom[h].N0 for h in hydrom_types],dtype='float32')
            lambdas=np.array([list_hydrom[h].lambda_ for h in hydrom_types],dtype='float32')
            mu=np.array([list_hydrom[h].mu for h in hydrom_types],dtype='float32')
            nu=np.array([list_hydrom[h].nu for h in hydrom_types],dtype='float32')
            
            # Compute RCS for all hydrometeors
            for j,h in enumerate(hydrom_types):
                sz = lut_sz[h].lookup_line(e = elev_lut[i],t = beam.values['T'][i])                         
                # get RCS
                rcs[:,j]= (2*np.pi*(sz[:,0] - sz[:,1] - sz[:,2] + sz[:,3])).T
            # Import we use symetrical elevations only for lookup querying, not
            # for actual trigonometrical velocity estimation
            Da, Db, idx = get_diameter_from_rad_vel(list_hydrom,phi,
                          elev[i],beam.values['U'][i],
                          beam.values['V'][i],beam.values['W'][i],rho_corr[i]) 
            try:
                refl[i,idx] = doppler_c.get_refl(len(idx),Da, 
                                Db,D,rcs,N0,mu,lambdas,nu,step_D,D_min)[1]
                wavelength = constants.WAVELENGTH
                refl[i,idx] *= wavelength**4/(np.pi**5*constants.KW**2)
            except:
                print('An error occured in the Doppler spectrum calculation...')
                raise

    return refl
    
def get_turb_std(ranges,EDR):
    sigma_r=cfg.CONFIG['radar']['radial_resolution']
    sigma_theta=cfg.CONFIG['radar']['3dB_beamwidth']*DEG2RAD # Convert to rad
    
    turb_std=np.zeros((len(EDR),))
    
    # Method of calculation follow Doviak and Zrnic (p.409)
    # Case 1 : sigma_r << r*sigma_theta
    idx_r=sigma_r<0.1*ranges*sigma_theta
    turb_std[idx_r]=((ranges[idx_r]*EDR[idx_r]*sigma_theta*constants.A**(3/2))/0.72)**(1/3)
    # Case 2 : r*sigma_theta <= sigma_r
    idx_r=sigma_r>=0.1*ranges*sigma_theta
    turb_std[idx_r]=(((EDR[idx_r]*sigma_r*(1.35*constants.A)**(3/2))/(11/15+4/15*(ranges[idx_r]*sigma_theta/sigma_r)**2)**(-3/2))**(1/3))
    return turb_std
    
def turb_spectrum_spread(spectrum, turb_std):
    v=constants.VARRAY
    # Convolve spectrum and turbulence gaussian distributions
    # Get resolution in velocity
    v_res=v[2]-v[1]
    
    original_power=np.sum(spectrum,1) # Power of every spectrum (at all radar gates)
    spectrum=gaussian_filter(spectrum,[0, turb_std/v_res]) # Filter only columnwise (i.e. on velocity bins)
    convolved_power=np.sum(spectrum,1) # Power of every convolved spectrum (at all radar gates)
    
    spectrum=spectrum/convolved_power[:,None]*original_power[:,None]# Rescale to original power
    
    return spectrum

def get_diameter_from_rad_vel(list_hydrom, phi,theta,U,V,W,rho_corr):
    theta=theta*DEG2RAD
    phi=phi*DEG2RAD

    wh=1./rho_corr*(W+(U*np.sin(phi)+V*np.cos(phi))/np.tan(theta)\
        -constants.VARRAY/np.sin(theta))
        
    idx=np.where(wh>=0)[0] # We are only interested in positive fall speeds 

    wh=wh[idx]
    
    hydrom_types = list_hydrom.keys()
    
    D=np.zeros((len(idx),len(hydrom_types)), dtype='float32')
    
    # Get D bins from V bins
    for i,h in enumerate(list_hydrom.keys()): # Loop on hydrometeors
        D[:,i]=list_hydrom[h].get_D_from_V(wh)
        # Threshold to valid diameters
        D[D>=list_hydrom[h].d_max]=list_hydrom[h].d_max
        D[D<=list_hydrom[h].d_min]=list_hydrom[h].d_min
    
    # Array of left bin limits
    Da=np.minimum(D[0:-1,:],D[1:,:])
    # Array of right bin limits
    Db=np.maximum(D[0:-1,:],D[1:,:])

    # Get indice of lines where at least one bin width is larger than 0
    mask=np.where(np.sum((Db-Da)==0.0,axis=1)<len(hydrom_types))[0]
    
    return Da[mask,:], Db[mask,:],idx[mask]

 
if __name__=='__main__':
    
    
    import pickle
    import gzip
    from cosmo_pol.lookup.lut import Lookup_table
    beams = pickle.load(open('/media/wolfensb/Storage/cosmo_pol/ex_beams_rhi.txt','rb'))
#    lut = pickle.load(gzip.open('/media/wolfensb/Storage/cosmo_pol/lookup/final_lut/all_luts_sz_f_2_7_1mom.pz','rb'))
    
    from cosmo_pol.utilities import cfg
    cfg.init('./option_files/MXPOL_RHI.yml') # Initialize options with 'options_radop.txt'   
    
#    from cosmo_pol.utilities import tictoc
#    tictoc.tic()
    cfg.CONFIG['doppler']['scheme'] = 2
    tictoc.tic()    
    c = get_doppler_velocity(beams, lut)
    tictoc.toc()
    
#    from cosmo_pol.lookup.read_lut import get_lookup_tables
#    import pickle
#
#    l=pickle.load(open('../ex_beams_rhi.txt','rb'))
#    lut_pol,lut_rcs=get_lookup_tables('1mom',5.6,True)
#    
##
##    config['doppler_vel_method']=2
##    rvel=get_doppler_velocity(l,lut_rcs)
##    plt.figure()
##    plt.plot(rvel.values['v_radial'])   
##    
##    
#    cfg.CONFIG['doppler_scheme']=3
#    rvel=get_doppler_velocity(l,lut_rcs)
#
#    cfg.CONFIG['doppler_scheme']=2
#    rvel3=get_doppler_velocity(l,lut_pol)
#    
#    cfg.CONFIG['doppler_scheme']=1
#    rvel2=get_doppler_velocity(l)
#    import matplotlib.pyplot as plt
#    plt.plot(rvel.values['RVel'])
#    plt.hold(True)
#    plt.plot(rvel2.values['RVel'])
#    plt.plot(rvel3.values['RVel'])
#    plt.legend(['Dop','Unweighted','Weighted'])