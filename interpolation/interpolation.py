# -*- coding: utf-8 -*-
"""
Created on Fri Jul 24 17:19:22 2015

@author: wolfensb
"""
import numpy as np
import pyproj
import pickle
import interpolation_c
import cosmo_pol.pycosmo as pycosmo

from cosmo_pol.refraction import atm_refraction
from  cosmo_pol.utilities import cfg
from cosmo_pol.utilities.beam import Beam
from cosmo_pol.utilities import utilities

def integrate_GH_pts(list_GH_pts):
    num_beams=len(list_GH_pts)
    
    list_variables=list_GH_pts[0].values.keys()
    
    integrated_variables={}
    for k in list_variables:
        integrated_variables[k]=[float('nan')]
        sum_weights=0
        for i in list_GH_pts:
            sum_weights+=i.GH_weight
        for i in list_GH_pts:
            integrated_variables[k] = utilities.nansum_arr(integrated_variables[k],i.values[k]*i.GH_weight/sum_weights)
    
    # Get index of central beam
    idx_0=int(num_beams/2)
    
    # Sum the mask of all beams to get overall mask
    mask=np.zeros(num_beams,) # This mask serves to tell if the measured point is ok, or below topo or above COSMO domain
    for i,p in enumerate(list_GH_pts):
        mask=utilities.sum_arr(mask,p.mask) # Get mask of every Beam
    mask/=float(num_beams) # Larger than 1 means that every Beam is below TOPO, smaller than 0 that at least one Beam is above COSMO domain
    mask[np.logical_and(mask>=0,mask<1)]=0

    heights_radar=list_GH_pts[idx_0].heights_profile
    distances_radar=list_GH_pts[idx_0].dist_profile
    lats=list_GH_pts[idx_0].lats_profile
    lons=list_GH_pts[idx_0].lons_profile

    integrated_beam=Beam(integrated_variables,mask,lats, lons, distances_radar, heights_radar)
    return integrated_beam

def get_profiles_GH(dic_variables, azimuth, elevation, radar_range=0,N=0, list_refraction=0):
    
    list_variables=dic_variables.values()
    keys=dic_variables.keys()
    
    # Get options
    bandwidth_3dB=cfg.CONFIG['radar']['3dB_beamwidth']
    integration_scheme = cfg.CONFIG['integration']['scheme']
    refraction_method = cfg.CONFIG['refraction']['scheme']
    
    list_variables=dic_variables.values()
    keys=dic_variables.keys()
    
    if integration_scheme == 1: # Classical single gaussian scheme
        nh_GH = int(cfg.CONFIG['integration']['nh_GH'])
        nv_GH = int(cfg.CONFIG['integration']['nv_GH'])
        
        # Get GH points and weights
        sigma=bandwidth_3dB/(2*np.sqrt(2*np.log(2)))
    
        pts_hor, weights_hor=np.polynomial.hermite.hermgauss(nh_GH)
        pts_hor=pts_hor*np.sqrt(2)*sigma
       
        pts_ver, weights_ver=np.polynomial.hermite.hermgauss(nv_GH)
        pts_ver=pts_ver*np.sqrt(2)*sigma
    
        weights = np.outer(weights_hor,weights_ver)
        
        threshold=np.mean([(weights_hor[0]*weights_hor[int(nh_GH/2)])/(nv_GH*nh_GH), 
                           (weights_ver[0]*weights_ver[int(nv_GH/2)])/(nv_GH*nh_GH)])
        sum_weights=np.pi
        
        beam_broadening=nh_GH>1 or nv_GH>1 # Boolean for beam-broadening (if only one GH point : No beam-broadening)

    elif integration_scheme == 2: # Improved multi-gaussian scheme
        nr_GH = int(cfg.CONFIG['integration']['nr_GH']) 
        na_GL= int(cfg.CONFIG['integration']['na_GL']) 
        
        antenna_params = cfg.CONFIG['integration']['antenna_params']
        
        pts_ang, w_ang = np.polynomial.legendre.leggauss(na_GL)
        pts_rad, w_rad = np.polynomial.hermite.hermgauss(nr_GH)
        
        a_dB = antenna_params[:,0]
        mu = antenna_params[:,1]
        sigma = antenna_params[:,2]
        
        list_pts = []
        weights = []
        sum_weights = 0
        for i in range(nr_GH):
            for j in range(len(sigma)):
                for k in range(na_GL):
                    r = mu[j]+np.sqrt(2)*sigma[j]*pts_rad[i]
                    theta = np.pi * pts_ang[k] + np.pi
                    weight = np.pi*w_ang[k]*w_rad[i]*10**(0.1*a_dB[j])*np.sqrt(2)*sigma[j]*abs(r) # Laplacian
                    weights.append(weight)
                    
                    sum_weights+=weight

                    list_pts.append([r*np.cos(theta)+azimuth,r*np.sin(theta)+elevation])


        beam_broadening=nr_GH>1 or na_GL>1 # Boolean for beam-broadening (if only one GH point : No beam-broadening)
    
    elif integration_scheme == 3: # Real antenna, for testing only
        
        data = pickle.load(open('real_antenna_s.p','rb'))
        angles = data['angles']
        
        pts_hor = angles
        pts_ver = angles
        
        threshold = -np.Inf
        beam_broadening = True
        
        weights = data['data']
        sum_weights= np.sum(weights)
    
    
    list_beams = []     
    # create vector of bin positions
    rranges=np.arange(cfg.CONFIG['radar']['radial_resolution']/2,
                          cfg.CONFIG['radar']['range'],
                          cfg.CONFIG['radar']['radial_resolution'])       
                          
    if integration_scheme == 1 or integration_scheme == 3:
        
        if list_refraction == 0: # Calculate refraction for vertical GH points
            list_refraction=[]
                                  
            # Get coordinates of virtual radar
            radar_pos=cfg.CONFIG['radar']['coords']

            for pt in pts_ver:
                if cfg.CONFIG['radar']['type'] == 'GPM':
                    S,H, E = atm_refraction.get_GPM_refraction(pt+elevation)
                else:
                    S,H, E = atm_refraction.get_radar_refraction(rranges, pt+elevation, radar_pos, refraction_method, N)
                list_refraction.append((S,H,E))
        
        for i in range(len(pts_hor)): 
            for j in range(len(pts_ver)):
                if weights[i,j]>threshold or not beam_broadening:
                    
                    # GH coordinates
                    pt=[pts_hor[i]+azimuth,pts_ver[j]+elevation]
                    weight=weights[i,j]/sum_weights
                    # Interpolate beam
                    lats,lons,b=get_radar_beam_trilin(list_variables, pts_hor[i]+azimuth, list_refraction[j][0],list_refraction[j][1])
                    
                    # Create dictionary of beams
                    dic_beams={}
                    for k, bi in enumerate(b): # Loop on interpolated variables
                        if k == 0: # Do this only for the first variable (same mask for all variables)
                            mask_beam = np.zeros((len(bi)))
                            mask_beam[bi == -9999] =- 1 # Means that the interpolated point is above COSMO domain
                            mask_beam[np.isnan(bi)] = 1  # NaN means that the interpolated point is below COSMO terrain
                        bi[mask_beam!=0] = float('nan') # Assign NaN to all missing data
                        dic_beams[keys[k]] = bi # Create dictionary
                    list_beams.append(Beam(dic_beams, mask_beam, lats, lons, list_refraction[j][0],list_refraction[j][1],list_refraction[j][2],pt, weight))        
        
    elif integration_scheme == 2:
        # create vector of bin positions
        # Get coordinates of virtual radar
    
        radar_pos=cfg.CONFIG['radar']['coords']
        
        for i in range(len(list_pts)):
            if cfg.CONFIG['radar']['type'] == 'GPM':
                S,H, E = atm_refraction.get_GPM_refraction(list_pts[i][1])
            else:
                S,H, E = atm_refraction.get_radar_refraction(rranges, list_pts[i][1], radar_pos, refraction_method, N)
       
            lats,lons,b=get_radar_beam_trilin(list_variables, list_pts[i][0], S,H)
            # Create dictionary of beams
            dic_beams={}
            for k, bi in enumerate(b): # Loop on interpolated variables
                if k == 0: # Do this only for the first variable (same mask for all variables)
                    mask_beam=np.zeros((len(bi)))
                    mask_beam[bi==-9999]=-1 # Means that the interpolated point is above COSMO domain
                    mask_beam[np.isnan(bi)]=1  # NaN means that the interpolated point is below COSMO terrain
                bi[mask_beam!=0]=float('nan') # Assign NaN to all missing data
                dic_beams[keys[k]]=bi # Create dictionary
            list_beams.append(Beam(dic_beams, mask_beam, lats, lons, S,H,E,list_pts[i][1], weights[i]/sum_weights))       
    
    
    return list_beams
    
    
def get_radar_beam_trilin(list_vars, azimuth, distances_profile, heights_profile):
    # Get position of virtual radar from cfg.CONFIG
    radar_pos=cfg.CONFIG['radar']['coords']

    # Initialize WGS84 geoid
    g = pyproj.Geod(ellps='WGS84')

    # Get radar bins coordinates
    lons_rad=[]
    lats_rad=[]
    # Using the distance on ground of every radar gate, we get its latlon coordinates
    for d in distances_profile:
        lon,lat,ang=g.fwd(radar_pos[1],radar_pos[0],azimuth,d) # Note that pyproj uses lon, lat whereas I used lat, lon
        lons_rad.append(lon)
        lats_rad.append(lat)

    # Convert to numpy array
    lons_rad=np.array(lons_rad)
    lats_rad=np.array(lats_rad)
    
    # Initialize interpolated beams
    all_beams=[]
    
    # Get model heights and COSMO proj from first variable    
    ###########################################################################
    model_heights=list_vars[0].attributes['z-levels']
    rad_interp_values=np.zeros(len(distances_profile),)*float('nan')
    
    # Get COSMO local coordinates info
    proj_COSMO=list_vars[0].attributes['proj_info']
    # Get lower left corner of COSMO domain in local coordinates
    llc_COSMO=(float(proj_COSMO['Lo1']), float(proj_COSMO['La1']))
    res_COSMO=list_vars[0].attributes['resolution']

    # Get resolution 
    # Transform radar WGS coordinates into local COSMO coordinates

    coords_rad_loc=pycosmo.WGS_to_COSMO((lats_rad,lons_rad),[proj_COSMO['Latitude_of_southern_pole'],proj_COSMO['Longitude_of_southern_pole']])  
    llc_COSMO=np.asarray(llc_COSMO).astype('float32')
    

    # Now we interpolate all variables along beam using C-code file
    ###########################################################################
    for n,var in enumerate(list_vars):           

        model_data=var.data
        rad_interp_values=interpolation_c.get_all_radar_pts(len(distances_profile),coords_rad_loc,heights_profile,model_data,model_heights\
        , llc_COSMO,res_COSMO)
        all_beams.append(rad_interp_values[1][:])

    return lats_rad, lons_rad, all_beams


if __name__=='__main__':
#    import matplotlib.pyplot as plt
    import cfg.CONFIG
    cfg.CONFIG.init('./sample_option_files/options_ALBIS_radar_RHI.txt') # Initialize options with 'options_radop.txt'
#    from rad_wind import get_doppler_velocity
    file_h=pc.open_file('/ltedata/COSMO/case2014040802_PAYERNE_analysis_ONEMOM/lfsf00124000')

    dic_vars=pc.get_variables(file_h,['QR_v','QS_v','QG_v','U','V','W','T'],get_proj_info=True,shared_heights=True,assign_heights=True,c_file='/ltedata/COSMO/case2014040802_PAYERNE_analysis_ONEMOM/lfsf00000000c')
    list_GH_pts = get_profiles_GH(dic_vars,0, 3)
    import pickle
    pickle.dump(list_GH_pts,open('ex_bemas.txt','wb'))
#    results1=[]
#    results2=[]
#    list_GH_pts = get_profiles_GH(dic_vars,-90, 10)
#    for az in np.arange(0,1.5,1.5):
#        
#        list_GH_pts = get_profiles_GH(dic_vars,az, 10)
##        dop_vel, spectrum=get_doppler_velocity(list_GH_pts)
#        
#        results1.append(list_GH_pts[int(len(list_GH_pts)/2)].values['QR_v'])
#    a=np.asarray(results1)
#    plt.figure()
#
#    plt.imshow(a)
    
