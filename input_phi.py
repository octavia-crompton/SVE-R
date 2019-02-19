import numpy as np
import scipy as sp
import os, sys

def make_phi(params):
        """
        
        """
        import numpy as np
        vG = params['vG']
  
        zmax = params['zmax']
        dz = params['dz']
        zseal = params['zseal']
        
        z =  np.arange(0., zmax+dz, dz)                
        nz = z.shape[0]                                         
            
        n  = vG['lambda'] + 1
        m  = vG['lambda']/n
    
        alpha = np.ones(nz)*vG['alpha']
        theta_S = np.ones(nz)*vG['theta_S']
        theta_R = np.ones(nz)*vG['theta_R']
        lambdA = np.ones(nz)*vG['lambda']
        n = np.ones(nz)*n
        m = np.ones(nz)*m
        
        if 'KsV' in params:
          ksat = np.ones(nz)*params['KsV']/3600.
        else:
          ksat = np.ones(nz)*vG['Ks']/3600.

        phi_veg = {'alpha': alpha,
               'theta_R': theta_R,
               'theta_S': theta_S,
               'lambda': lambdA,
               'n': n, 
               'm': m, 
               'ksat': ksat,
              }
      
        # seal 
        si =  np.where(z == z[-1] - zseal)[0][0]+1  
        
        
        ksatB = ksat.copy()
        ksatB[si:] = params['KsB']/3600.
        
        phi_bare = {'alpha': alpha,
           'theta_R': theta_R,
           'theta_S': theta_S,
           'lambda': lambdA,
           'n': n, 
           'm': m, 
           'ksat': ksatB,
          }
         
        return phi_veg, phi_bare, nz

        
def write_phi(path,phi_veg,phi_bare, h_init,
              dz, zmax):
        """
        writes soil vanG parameters and initial conditions to 
        """
        nz = h_init.shape[0]            
        f = open('{0}/vanG.dat'.format(path), 'w')      
        f.write('dz zmax \n')  
        f.write('{0:<13} {1:<13}  \n'.format(dz, zmax))
        f.write('nz \n')  
        f.write('{0:<13}  \n'.format(nz))    
        f.write("alpha    theta_S   theta_R    lambda    Ksat                   h_init \n")        
        f.write(' vegetated cells  \n'.format(nz))  
        for ndum in range(nz):
            f.write('{0}   {1}     {2}     {3}      {4}      {5} \n'.format(phi_veg['alpha'][ndum],
                            phi_veg['theta_S'][ndum], phi_veg['theta_R'][ndum], 
                            phi_veg['lambda'][ndum], phi_veg['ksat'][ndum],
                            h_init[ndum]))
        
        f.write(' bare soil cells \n'.format(nz))
        for ndum in range(nz):
            f.write('{0}   {1}     {2}     {3}      {4}      {5} \n'.format(phi_bare['alpha'][ndum],
                            phi_bare['theta_S'][ndum], phi_bare['theta_R'][ndum], 
                            phi_bare['lambda'][ndum], phi_bare['ksat'][ndum],
                            h_init[ndum]))
        
        f.close()

    

def vanGenuchten(h,phi) :
    """
    """
    import numpy as np
    alpha   = phi['alpha']
    theta_S = phi['theta_S']
    theta_R = phi['theta_R']
    n       = phi['n']
    m       = phi['m']
    Ksat    = phi['ksat'] 
    L = phi['L']
    # Compute the volumetric moisture content
    theta = (theta_S - theta_R)/(1 + (alpha*abs(h))**n)**m + theta_R
    # Compute the effective saturation
    Se = ((theta - theta_R)/(theta_S - theta_R))
    # Compute the hydraulic conductivity
    K = Ksat*Se**(L)*(1 - (1 - Se**(1./m))**m)**2
    # Compute the specific moisture storage
    C =  -alpha*n*np.sign(h)*(1./n - 1)*(alpha*abs(h))**(n - 1)*(theta_R - 
         theta_S)*((alpha*abs(h))**n + 1)**(1/n - 2)

    if type(h) == float:
      h = [h]
    for i in range(len(h)):
        if h[i] > 0:
            K[i] = Ksat[i]
            C[i] = 0.
            theta[i] = theta_S[i]

    return [C,K,theta]


