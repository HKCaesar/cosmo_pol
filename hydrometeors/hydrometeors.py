
import numpy as np
from scipy.integrate import trapz

import dielectric

from cosmo_pol.constants import constants
from cosmo_pol.constants import constants_1mom
from cosmo_pol.constants import constants_2mom

###############################################################################

def create_hydrometeor(hydrom_type,scheme, sedi = False):
    if  hydrom_type == 'R':
       return Rain(scheme, sedi)
    elif hydrom_type == 'S':
        return Snow(scheme)
    elif hydrom_type == 'G':
        return Graupel(scheme)
    elif hydrom_type == 'H':
        return Hail(scheme)
        
###############################################################################

class Hydrometeor(object):
    def __init__(self, scheme):
        if scheme not in ['1mom','2mom']:
            scheme = '1mom'
            
        self.scheme = scheme
        self.d_max=None # Max diam in the integration over D
        self.d_min=None # Min diam in the integration over D
        
        # Power-law parameters
        self.a = None
        self.b = None
        self.alpha = None
        self.beta = None
        
        # PSD parameters
        self.lambda_ = None
        self.N0 = None
        self.mu = None
        self.nu = None
        
        # Scattering parameters
        self.canting_angle_std = None
        
        # Integration factors
        self.lambda_factor=None
        self.vel_factor=None
        self.ntot_factor=None
        
        if self.scheme == '2mom':
            self.x_max = None
            self.x_min = None

    def get_N(self,D): 
        try:
            siz = self.lambda_.shape
            D=np.repeat(D[:,None],siz[0],axis=1)
            if(len(siz)>1):
                D=(np.repeat(D[:,:,None],siz[1],axis=2))
            return np.squeeze(self.N0*D**self.mu*np.exp(-self.lambda_*D**self.nu))
        except:
            raise
            print('Wrong dimensions of input')
            
    def get_V(self,D):
        return self.alpha*D**self.beta
        
    def get_D_from_V(self,V):
        return (1.0/self.alpha*V)**(1.0/self.beta)
        
    def integrate_V(self):
        V_int = self.vel_factor*self.N0*self.alpha/self.nu*self.lambda_**(-(self.beta+self.mu+1)/self.nu)
        if self.scheme == '2mom':
            N_int = self.ntot
        else:
            N_int = self.ntot_factor*self.N0/self.nu*self.lambda_**(-(self.mu+1)/self.nu)
        return V_int, N_int
        
    def integrate_V_weighted(self,D,weights, method='sum'):
        try:
            if method not in ['sum','trapz']:
                print('Invalid integration method')
                print('Choosing sum instead')
                method = 'sum'
                
            V=self.get_V(D)
            N=self.get_N(D)
                
            dD = D[1]-D[0]
            
            if method == 'sum':
                # Fast version using einsum
                if len(N.shape) == 1:
                    V_weighted = np.sum(N*weights*V) * dD
                    N_weighted = np.sum(N*weights) * dD
                elif len(N.shape) == 2:
                    V_weighted = np.einsum('ij,i,i,->j',N,weights,V) * dD
                    N_weighted = np.einsum('ij,i->j',N,weights) * dD
                elif len(N.shape) == 3:
                    V_weighted = np.einsum('ijk,i,i->jk',N,weights,V) * dD
                    N_weighted = np.einsum('ijk,i->jk',N,weights) * dD 
                    
            elif method == 'trapz':
                # More precise but slow (using trapz)
                if len(N.shape) == 2:
                    V=np.repeat(V[:,None],N.shape[1],axis=1)
                    weights=np.repeat(weights[:,None],N.shape[1],axis=1)
                elif len(N.shape) == 3:
                    V=np.repeat(np.repeat(V[:,None,None],N.shape[1],axis=1),N.shape[2],axis=2)
                    weights=np.repeat(np.repeat(weights[:,None,None],N.shape[1],axis=1),N.shape[2],axis=2)  
                    
                V_weighted = trapz(N*weights*V,dx=dD,axis=0)# multiply by integration step
                N_weighted = trapz(N*weights,dx=dD,axis=0)# multiply by integration step
                    
            return V_weighted, N_weighted
        except:
            raise
            print('Wrong dimensions of input')
            
    def get_M(self,D):
        return self.a*D**self.b
            
    def set_psd(self,*args):
        if len(args) == 2 and self.scheme == '2mom':
            # First argument must be number density and second argument mass
            # density
            with np.errstate(divide='ignore'):
                qn = np.squeeze(np.array([args[0]]))
                q = np.squeeze(np.array([args[1]]))

                x_mean = np.minimum(np.maximum(q*1.0/(qn+constants.EPS),self.x_min),self.x_max)
                x_mean = np.squeeze(x_mean)

        
                _lambda = np.asarray((self.lambda_factor*x_mean)**(-self.nu/self.b))

                # This is done in order to make sure that lambda is always a 1D
                # array

                _lambda[q == 0] = float('nan')
                if _lambda.shape == ():
                    _lambda = np.array([_lambda])

                _N0 = np.asarray((self.nu/self.ntot_factor)*qn*_lambda**((self.mu+1)/self.nu))
                _N0 = _N0 * 1000**(-(1+self.mu))
                
                _lambda = _lambda * 1000**(-self.nu)
            
                self.N0 = _N0.T
                self.lambda_ = _lambda.T
                self.ntot = args[0]
        
            
###############################################################################
            
class Solid(Hydrometeor):
    def get_fractions(self,D):
        # Uses COSMO mass-diameter rule
        f_ice = 6*self.a/(np.pi*constants.RHO_I)*D**(self.b-3)
        f_air = 1-f_ice
        return [f_ice, f_air]
        
    def get_m_func(self,T,f):
        def func(D, T, f):
            frac = self.get_fractions(D)
            m_ice = dielectric.dielectric_ice(T,f)
            return dielectric.dielectric_mixture(frac, [m_ice, constants.M_AIR])
        return lambda D: func(D, T, f)

        
###############################################################################

class Rain(Hydrometeor):
    def __init__(self, scheme, sedi=True):

        if scheme not in ['1mom','2mom']:
            scheme = '1mom'
            
        self.scheme = scheme
        self.sedi = (sedi if self.scheme == '2mom' else False)
        self.d_max = constants_2mom.D_MAX_R
        self.d_min = constants_2mom.D_MIN_R        
        
        # Power-law parameters
        self.a = (constants_1mom.A_R if self.scheme == '1mom' else constants_2mom.A_R)
        self.b = (constants_1mom.B_R if self.scheme == '1mom' else constants_2mom.B_R)
        self.alpha = (constants_1mom.ALPHA_R if self.scheme == '1mom' else constants_2mom.ALPHA_R)
        self.beta = (constants_1mom.BETA_R if self.scheme == '1mom' else constants_2mom.BETA_R)
        
        # PSD parameters
        self.lambda_ = None
        self.N0 = (constants_1mom.N0_R if self.scheme == '1mom' else None)
        self.mu = (constants_1mom.MU_R if self.scheme == '1mom' else constants_2mom.MU_R)
        self.nu = (1 if self.scheme == '1mom' else constants_2mom.NU_R)
        
        # Scattering parameters
        self.canting_angle_std = 10.
        
        # Others
        self.lambda_factor = (constants_1mom.LAMBDA_FACTOR_R if self.scheme == '1mom' else constants_2mom.LAMBDA_FACTOR_R)
        self.vel_factor = (constants_1mom.VEL_FACTOR_R if self.scheme == '1mom' else constants_2mom.VEL_FACTOR_R)
        self.ntot_factor=(constants_1mom.NTOT_FACTOR_R if self.scheme == '1mom' else constants_2mom.NTOT_FACTOR_R)
                
        if self.scheme == '2mom':
            self.x_max = constants_2mom.X_MAX_R
            self.x_min = constants_2mom.X_MIN_R
        else:
            self.vel_factor = constants_1mom.VEL_FACTOR_R

    def get_V(self,D):
        if self.sedi:
            # Specific V-D relation for sedimentation
            return constants_2mom.C_1-constants_2mom.C_2*np.exp(-constants_2mom.C_3*D)
        else:
            return super(Rain,self).get_V(D) # Call general hydrometeor method for 

    def get_D_from_V(self,V):
        if self.sedi:
            # Specific V-D relation for sedimentation
            return -1/constants_2mom.C_3*np.log((constants_2mom.C_1-V)/constants_2mom.C_2)
        else:
            return super(Rain,self).get_D_from_V(V)
        
    def integrate_V(self,*args):
        if not self.sedi or not self.scheme == '2mom':
            return super(Rain,self).integrate_V()
        else:
            D=args[0]
            dD=D[1]-D[0]
            return np.einsum('ij,i->j',self.get_N(D),self.get_V(D))*dD,self.ntot
        
    def set_mu(self,QN,QM,QC):
        if self.sedi and self.scheme == '2mom':
            if QC<constants.EPS:
                D_mean=self.a*(QN/QN)**self.b*1000. # To get in mm (and not m)
                if D_mean <= constants_2mom.D_EQ:
                    mu = 9*np.tanh((constants_2mom.TAU_1*(D_mean-constants_2mom.D_EQ))**2)+1
                else:
                    mu = 9*np.tanh((constants_2mom.TAU_1*(D_mean-constants_2mom.D_EQ))**2)+1
                    
                self.mu = mu
                
    def set_psd(self,*args):
        if len(args) == 2 and self.scheme == '2mom':
            super(Rain,self).set_psd(*args)    
        elif self.scheme == '1mom':
            with np.errstate(divide='ignore'):
                _lambda = np.array((self.lambda_factor/args[0])**(1./(4.+self.mu)))
                _lambda[args[0]==0] = float('nan')
                self.lambda_ = _lambda
        else:
            print('Invalid call to function, if scheme == ''2mom'', input must be tuple of (QN,QM)')     
            print('if scheme == ''2mom'', input must be (QM)')
    
    def get_axis_ratio(self,D):
        ar = np.zeros((len(D),))
        
        ar[D<0.7] = 1.0
        mid_diam = np.logical_and(D<1.5,D>=0.7)
        ar[mid_diam] = 1.173 - 0.5165*D[mid_diam] + 0.4698*D[mid_diam]**2 - 0.1317*D[mid_diam]**3 - \
                8.5e-3*D[mid_diam]**4
        ar[D>=1.5] = 1.065 - 6.25e-2*D[D>=1.5] - 3.99e-3*D[D>=1.5]**2 + 7.66e-4*D[D>=1.5]**3 - \
                4.095e-5*D[D>=1.5]**4 
            
        return 1./ar

        
    def get_m_func(self,T,f):
        return lambda D: dielectric.dielectric_water(T, f)
    
###############################################################################
    
class Snow(Solid):
    def __init__(self, scheme):
        if scheme not in ['1mom','2mom']:
            scheme = '1mom'
            
        self.scheme = scheme
        self.d_max = constants_2mom.D_MAX_S
        self.d_min = constants_2mom.D_MIN_S  
        
        # Power-law parameters
        self.a = (constants_1mom.A_S if self.scheme == '1mom' else constants_2mom.A_S)
        self.b = (constants_1mom.B_S if self.scheme == '1mom' else constants_2mom.B_S)
        self.alpha = (constants_1mom.ALPHA_S if self.scheme == '1mom' else constants_2mom.ALPHA_S)
        self.beta = (constants_1mom.BETA_S if self.scheme == '1mom' else constants_2mom.BETA_S)
        
        # PSD parameters
        self.lambda_ = None
        self.N0 = None
        self.mu = (constants_1mom.MU_S if self.scheme == '1mom' else constants_2mom.MU_S)
        self.nu = (1 if self.scheme == '1mom' else constants_2mom.NU_S)
        
        # Scattering parameters
        self.canting_angle_std = 20.
        
        # Others
        self.lambda_factor = (constants_1mom.LAMBDA_FACTOR_S if self.scheme == '1mom' else constants_2mom.LAMBDA_FACTOR_S)
        self.vel_factor = (constants_1mom.VEL_FACTOR_S if self.scheme == '1mom' else constants_2mom.VEL_FACTOR_S)
        self.ntot_factor=(constants_1mom.NTOT_FACTOR_S if self.scheme == '1mom' else constants_2mom.NTOT_FACTOR_S)
        
        if self.scheme == '2mom':
            self.x_max = constants_2mom.X_MAX_S
            self.x_min = constants_2mom.X_MIN_S
    
            
    def set_psd(self,*args):
        if len(args) == 2 and self.scheme == '2mom':
            super(Snow,self).set_psd(*args)    
            
        elif len(args) == 2 and self.scheme == '1mom':
            # For N0 use relation by Field et al. 2005 (QJRMS)
            self.N0 = 13.5*(5.65*10**5*np.exp(-0.107*(args[0]-273.15)))/1000 # mm^-1 m^-3
            with np.errstate(divide='ignore'):
                _lambda=np.array((self.a*self.N0*self.lambda_factor/args[1])**(1./(self.b+1))) # in m-1
                _lambda[args[1]==0]=float('nan')
                self.lambda_ = _lambda
        else:
            print('Invalid call to function, if scheme == ''2mom'', input must be tuple of (QN,QM)')   
            print('if scheme == ''2mom'', input must be tuple of (T,QM)')
            
    def get_axis_ratio(self,D):
        ar=(0.01714*D+0.8467) # Brandes et al 2007 (Colorado snowstorms)
        return 1.0/ar
    
    def get_axis_ratio_pdf_masc(self,D):
        alpha = constants.A_AR_ALPHA_AGG*D**constants.B_AR_ALPHA_AGG
        loc = constants.A_AR_LOC_AGG*D**constants.B_AR_LOC_AGG
        scale = constants.A_AR_SCALE_AGG*D**constants.B_AR_SCALE_AGG
        return (alpha,loc,scale)
        

    def get_canting_angle_std_masc(self,D):
        cant_std = constants.A_CANT_STD_AGG*D**constants.B_CANT_STD_AGG
        return cant_std  
        
###############################################################################


class Graupel(Solid):
    def __init__(self, scheme):
        if scheme not in ['1mom','2mom']:
            scheme = '1mom'
            
        self.scheme = scheme
        self.d_max = constants_2mom.D_MAX_G
        self.d_min = constants_2mom.D_MIN_G  
        
        # Power-law parameters
        self.a = (constants_1mom.A_G if self.scheme == '1mom' else constants_2mom.A_G)
        self.b = (constants_1mom.B_G if self.scheme == '1mom' else constants_2mom.B_G)
        self.alpha = (constants_1mom.ALPHA_G if self.scheme == '1mom' else constants_2mom.ALPHA_G)
        self.beta = (constants_1mom.BETA_G if self.scheme == '1mom' else constants_2mom.BETA_G)
        
        # PSD parameters

        self.lambda_ = None
        self.N0 = constants_1mom.N0_G
        self.mu = (constants_1mom.MU_G if self.scheme == '1mom' else constants_2mom.MU_G)
        self.nu = (1 if self.scheme == '1mom' else constants_2mom.NU_G)
        
        # Scattering parameters
        self.canting_angle_std = 40.
        
        # Others
        self.lambda_factor = (constants_1mom.LAMBDA_FACTOR_G if self.scheme == '1mom' else constants_2mom.LAMBDA_FACTOR_G)
        self.vel_factor = (constants_1mom.VEL_FACTOR_G if self.scheme == '1mom' else constants_2mom.VEL_FACTOR_G)
        self.ntot_factor=(constants_1mom.NTOT_FACTOR_G if self.scheme == '1mom' else constants_2mom.NTOT_FACTOR_G)
        
        if self.scheme == '2mom':
            self.x_max = constants_2mom.X_MAX_G
            self.x_min = constants_2mom.X_MIN_G


    def set_psd(self,*args):
        if len(args) == 2 and self.scheme == '2mom':
            super(Graupel,self).set_psd(*args)    
        elif self.scheme == '1mom':
            with np.errstate(divide='ignore'):
                _lambda = np.array((self.lambda_factor/args[0])**(1./(4.+self.mu)))
                _lambda[args[0]==0] = float('nan')
                self.lambda_ = _lambda
        else:
            print('Invalid call to function, if scheme == ''2mom'', input must be tuple of (QN,QM)')     
            print('if scheme == ''2mom'', input must be (QM)')
            
    def get_axis_ratio(self,D): # Garrett, 2015 http://onlinelibrary.wiley.com/doi/10.1002/2015GL064040/full
        ar=0.9*np.ones(len(D),)
        return 1.0/ar
    
    def get_axis_ratio_pdf_masc(self,D):
        alpha = constants.A_AR_ALPHA_GRAU*D**constants.B_AR_ALPHA_GRAU
        loc = constants.A_AR_LOC_GRAU*D**constants.B_AR_LOC_GRAU
        scale = constants.A_AR_SCALE_GRAU*D**constants.B_AR_SCALE_GRAU
        return (alpha,loc,scale)
        
    def get_canting_angle_std_masc(self,D):
        cant_std = constants.A_CANT_STD_GRAU*D**constants.B_CANT_STD_GRAU
        return cant_std       
        
###############################################################################

        
class Hail(Solid):
    def __init__(self,scheme='2mom'):

        self.scheme = '2mom' # No 1-moment scheme for hail
        self.d_max = constants_2mom.D_MAX_H
        self.d_min = constants_2mom.D_MIN_H  
        
        # Power-law parameters
        self.a = (constants_1mom.A_H if self.scheme == '1mom' else constants_2mom.A_H)
        self.b = (constants_1mom.B_H if self.scheme == '1mom' else constants_2mom.B_H)
        self.alpha = (constants_1mom.ALPHA_H if self.scheme == '1mom' else constants_2mom.ALPHA_H)
        self.beta = (constants_1mom.BETA_H if self.scheme == '1mom' else constants_2mom.BETA_H)
        
        # PSD parameters
        self.lambda_ = None
        self.N0 = None
        self.mu = (constants_1mom.MU_H if self.scheme == '1mom' else constants_2mom.MU_H)
        self.nu = (1 if self.scheme == '1mom' else constants_2mom.NU_H)
        
        # Scattering parameters
        self.canting_angle_std = 40.
        
        # Others
        self.lambda_factor = constants_2mom.LAMBDA_FACTOR_H
        self.vel_factor = constants_2mom.VEL_FACTOR_H
        self.ntot_factor=(constants_1mom.NTOT_FACTOR_H if self.scheme == '1mom' else constants_2mom.NTOT_FACTOR_H)
        
        self.x_max = constants_2mom.X_MAX_H
        self.x_min = constants_2mom.X_MIN_H

    def get_axis_ratio(self,D):
        ar = 0.9 * np.ones((len(D),))
        return 1.0/ar
     
     

    
        
if __name__=='__main__':
    import pickle
    import gzip
    plt.close('all')
#    
#    lut1 = pickle.load(gzip.open('/media/wolfensb/Storage/cosmo_pol/lookup/stored_lut/lut_SZ_S_9_41_1mom.pz','rb'))
#    dd2 = lut1.value_table[10,10,:,0]
#    
#    lut2 = pickle.load(gzip.open('/media/wolfensb/Storage/cosmo_pol/lookup/stored_lut/lut_SZ_S_9_41_2mom.pz','rb'))
#    dd = lut2.value_table[10,10,:,0]
#        
#    plt.figure()
#    plt.plot(dd)
#    plt.hold(True)
#    plt.plot(dd2)
    
    
    
#    a=create_hydrometeor('R','2mom')
#    a.set_psd(np.ones(1,)*703,np.ones(1,)*0.0025270581)
    D = np.linspace(0.2,15,1024)    
#
##    v,n = a.integrate_V()
##    a.set_psd(np.array([10**-3]))
##    print(a.lambda_)
#    D = np.linspace(0.2,8,1024)
#    plt.plot(D,a.get_N(D))
#    plt.hold(True)
#    print(np.trapz(dd*a.get_N(D),D))
#    a=create_hydrometeor('R','1mom')
#    a.set_psd(np.ones(1,)*0.0025270581)
##    print(np.trapz(a.a*D**(a.b)*a.get_N(D),D))
#    plt.plot(D,a.get_N(D))    
##    print(np.trapz(a.a*D**(a.b)*a.get_N(D),D))
#    print(np.trapz(dd*a.get_N(D),D))

    
    #    frac = a.get_fractions(D)
#    plt.plot(D,frac[0],D,frac[1]    )
#    plt.hold(True)
#
    a1=create_hydrometeor('S','1mom')
    frac1 = a1.get_fractions(D)
    diel1 = a1.get_m_func(270,9.41)
    plt.plot(D,frac1[0],D,frac1[1]   )
    plt.hold(True)
    a2=create_hydrometeor('S','2mom')
    frac2 = a2.get_fractions(D)
    diel2 = a2.get_m_func(270,9.41)
    plt.plot(D,frac2[0],D,frac2[1]   ) 
    plt.legend(['ice1','air1','ice2','air2'])

    plt.figure()
    plt.plot(D,a1.a*D**a1.b)
    plt.plot(D,a2.a*D**a2.b)    
#    
    plt.figure()
    plt.plot(D,diel1(D),D,diel2(D))
    
    
    
    from pytmatrix import orientation
    from pytmatrix.tmatrix import Scatterer
    
    f=9.41
    wavelength=constants.C/(f*1E09)*1000 # in mm
    scatt = Scatterer(radius = 1.0, wavelength = wavelength)
    scatt.or_pdf = orientation.gaussian_pdf(std=20)
    scatt.orient = orientation.orient_averaged_fixed
    
    list_SZ_1=[]
    list_SZ_2=[]
    
    elevation = 5
    m_func1= a1.get_m_func(270,f)    
    m_func2= a2.get_m_func(270,f)     
    geom_back=(90-elevation, 180-(90-elevation), 0., 180, 0.0,0.0)
    geom_forw=(90-elevation, 90-elevation, 0., 0.0, 0.0,0.0)
    

    for i,d in enumerate(D):
        
        ar = a1.get_axis_ratio(d)
        
        scatt.radius = d/2.
        scatt.m = m_func1(d)
        scatt.axis_ratio = ar

        # Backward scattering (note that we do not need amplitude matrix for backward scattering)
        scatt.set_geometry(geom_back)
        Z_back = scatt.get_Z()
        # Forward scattering (note that we do not need phase matrix for forward scattering)
        scatt.set_geometry(geom_forw)
        S_forw=scatt.get_S()

        list_SZ_1.append(Z_back[0,0])

        ar = a2.get_axis_ratio(d)
        
        scatt.radius = d/2.
        scatt.m = m_func2(d)
        scatt.axis_ratio = ar

        # Backward scattering (note that we do not need amplitude matrix for backward scattering)
        scatt.set_geometry(geom_back)
        Z_back=scatt.get_Z()
        # Forward scattering (note that we do not need phase matrix for forward scattering)
        scatt.set_geometry(geom_forw)
        S_forw=scatt.get_S()

        list_SZ_2.append(Z_back[1,1])    
        
#    a.set_psd(np.array([10**-4]))
#    print(a.lambda_)
#    N2= a.get_N(np.linspace(0.1,6,100))    
#    
#    a.set_psd(np.array([0.5*(10**-4+10**-3)]))
#    N4= a.get_N(np.linspace(0.1,6,100))    
#    
#    a.lambda_ = np.array([(1.57084824+2.62033279)/2])
#    N5= a.get_N(np.linspace(0.1,6,100))    
#    plt.plot(D,0.5*(N+N2),D,N4,D,N5)
#    QN = 10**(np.linspace(0,8,150))
#    Q= 10**(np.linspace(-8,-2,150))
#    a.set_psd(QN,Q)
#    lambdas = a.lambda_
#    n0s = a.N0
#    N= a.get_N(np.linspace(0.1,6,100))
#    
#    uuu = a.integrate_V_weighted(np.linspace(0.1,6,100),np.ones((100,)))
#    
#    uuu=uuu[0]/uuu[1]
#    lambdas_2 = np.zeros((len(QN),len(Q)))
#    uuu_2 =  np.zeros((len(QN),len(Q)))
#    n0s_2 = np.zeros((len(QN),len(Q)))
#    for i,qi in enumerate(QN):
#        for j,qj in enumerate(Q):
#            a.set_psd(qi,qj)
#            lambdas_2[i,j]=a.lambda_
#            n0s_2[i,j]=a.N0
#            N= a.get_N(np.linspace(0.1,6,100))
#            fu= a.integrate_V_weighted(np.linspace(0.1,6,100),np.ones((100,)))
#            uuu_2[i,j]=fu[0]/fu[1]
#            
##    print a.get_D_from_V(np.array([6]*20))
#    a.set_psd(1000,0.01)
#    D=np.arange(0,10,0.01)
#    N=a.get_N(D)
#    W=np.random.rand(1000,)*0+1
##    
##    
#    a.set_psd(np.array([0.001]))
##    tic()
#    for i in range(30):
#        v=a.integrate_V_weighted(D,W)
#    toc()
##    ##print v[0]/v[1]
