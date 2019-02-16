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
        L  = vG['L']
    
        alpha = np.ones(nz)*vG['alpha']
        theta_s = np.ones(nz)*vG['theta_s']
        theta_r = np.ones(nz)*vG['theta_r']
        lambdA = np.ones(nz)*vG['lambda']
        n = np.ones(nz)*n
        m = np.ones(nz)*m
        L = np.ones(nz)*L
        
        if 'KsV' in params:
          ksat = np.ones(nz)*params['KsV']/3600.
        else:
          ksat = np.ones(nz)*vG['Ks']/3600.

        phi_veg = {'alpha': alpha,
               'theta_r': theta_r,
               'theta_s': theta_s,
               'lambda': lambdA,
               'n': n, 
               'm': m, 
               'ksat': ksat,
               'L' : L
              }
      
        # seal 
        si =  np.where(z == z[-1] - zseal)[0][0]+1  
        
        
        ksatB = ksat.copy()
        ksatB[si:] = params['KsB']/3600.
        
        phi_bare = {'alpha': alpha,
           'theta_r': theta_r,
           'theta_s': theta_s,
           'lambda': lambdA,
           'n': n, 
           'm': m, 
           'ksat': ksatB,
           'L' : L
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
        f.write('L \n')  
        f.write('{0:<13}  \n'.format(phi_veg['L'][0]))                
        f.write(' vegetated cells  \n'.format(nz))        
        for ndum in range(nz):
            f.write('{0} {1} {2} {3} {4} {5} \n'.format(phi_veg['alpha'][ndum],
                            phi_veg['theta_s'][ndum], phi_veg['theta_r'][ndum], 
                            phi_veg['lambda'][ndum], phi_veg['ksat'][ndum],
                            h_init[ndum]))
        
        f.write(' bare soil cells \n'.format(nz))
        for ndum in range(nz):
            f.write('{0} {1} {2} {3} {4} {5} \n'.format(phi_bare['alpha'][ndum],
                            phi_bare['theta_s'][ndum], phi_bare['theta_r'][ndum], 
                            phi_bare['lambda'][ndum], phi_bare['ksat'][ndum],
                            h_init[ndum]))
        
        f.close()

    

def vanGenuchten(h,phi) :
    """
    """
    import numpy as np
    alpha   = phi['alpha']
    theta_s = phi['theta_s']
    theta_r = phi['theta_r']
    n       = phi['n']
    m       = phi['m']
    Ksat    = phi['ksat'] 
    L = phi['L']
    # Compute the volumetric moisture content
    theta = (theta_s - theta_r)/(1 + (alpha*abs(h))**n)**m + theta_r
    # Compute the effective saturation
    Se = ((theta - theta_r)/(theta_s - theta_r))
    # Compute the hydraulic conductivity
    K = Ksat*Se**(L)*(1 - (1 - Se**(1./m))**m)**2
    # Compute the specific moisture storage
    C =  -alpha*n*np.sign(h)*(1./n - 1)*(alpha*abs(h))**(n - 1)*(theta_r - 
         theta_s)*((alpha*abs(h))**n + 1)**(1/n - 2)

    if type(h) == float:
      h = [h]
    for i in range(len(h)):
        if h[i] > 0:
            K[i] = Ksat[i]
            C[i] = 0.
            theta[i] = theta_s[i]

    return [C,K,theta]


