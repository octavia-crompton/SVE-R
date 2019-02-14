import numpy as np
import os, sys


def myfloat(b):
    
    try: 
        b = float(b)
    except ValueError:      
        if '-' in b:  
            b = [b for b in b.split('-') if b]
            b = float(b[0])*10**(-float(b[1]))
        elif '+' in b:
            b = 0.
    return b

def read_time(path = 'test') :  
    """
    read time.out 
    """
    t_p = []   
    cfl = []  #  print step
    f =  open("{0}/output/time.out".format(path), 'r')
    f.next()
    for line in f:
        a = (line.strip().split(" "))
        a = [b for b in a if b]
        t_p.append(float(a[0]))
        cfl.append(float(a[1]))

    t_p = np.array(t_p)
    t_p = np.round(t_p, 2)
    
    cfl = np.array(cfl)
    
    return t_p, cfl
    
def read_tp(path): 
    """
    tp : time of ponding 
    t_final :  final sim time 
    runtime : fortran runtime
    """  
    runtime = np.nan
    tp_for = np.inf
    t_final = np.nan 
    for line in open('{0}/summary.txt'.format(path), 'r'):
        a = line.strip().split(' ')
        a =  [b for b in a if b] 
        try:
          if a[0] == 'runtime':
            runtime = float(a[-1])            
          if a[0] == '1st':
            tp_for = float(a[-1])            
          if a[0] == 'final':                            
            t_final = float(a[-1])                          
        except:
          continue    
        
    return runtime, t_final, tp_for

def read_hydro(path ,  Lx=1, Ly=1):
    
    """ Reads hydrograph timeseries 
    hydro units are f*ds = m^3/s, converted to cm/s"""
    tdum = []
    hydro = []
    f =  open("{0}/output/hydro.out".format(path), 'r')
    f.next()
    for line in f:
        a = (line.strip().split(" "))
        a = [myfloat(b) for b in a if b]
        try:
            tdum.append(float(a[0]))
            hydro.append(float(a[1]))

        except IndexError:
            raise KeyboardInterrupt

    tdum = np.array(tdum)
    hydro = - np.array(hydro)/Lx/Ly*100  # now  cm/s   
    tdum = np.round(tdum, 2)
    
    return tdum, hydro # now cm/s


def read_ptsTheta(path, nz = 51 ):
    """
    reads H and Theta profiles for two selected points
    """  
    ThetaH = []
    
    ThetaHdum =  np.zeros([nz, 4])
    for line in open("{0}/output/ptsTheta.out".format(path), 'r'):
        a = (line.strip().split(" "))
        a = [myfloat(b) for b in a if b]
        
        if len(a) > 2:
            k = int(a[0])-1
            ThetaHdum[k, :] = a[2:]

        elif len(a)== 2:
            
            dumt = int(a[0])
            ThetaH.append(ThetaHdum.copy())    
    
    ThetaH = np.array(ThetaH)
    vegTheta = ThetaH[:, :, 0]
    vegH = ThetaH[:, :, 1]  
    bareTheta = ThetaH[:, :, 2]
    bareH = ThetaH[:, :, 3]  

    return vegTheta,vegH,bareTheta,bareH
    
    
def read_h(path, ncol , nrow , dx  ):
    """
    Reads the sve output variables from output/h.out
    Returns  h (m), u and v (m/s), infiltration (m2*cm) and intercell fluxes (m2*cm)
        with time resolution dt_p
      
     Fluxes are positive directed out of the domain
    """
    h = []
    hdum =  np.zeros([ncol, nrow])
    v = []
    vdum =  np.zeros([ncol, nrow])
    u = []
    udum =  np.zeros([ncol, nrow])

    zinflmap2 = []
    infldum =  np.zeros([ncol, nrow])
    
    xflux0 = []
    xfluxdum0 =  np.zeros([ncol, nrow])    
    
    yflux0 = []
    yfluxdum0 =  np.zeros([ncol, nrow])
    
    xflux1 = []
    xfluxdum1 =  np.zeros([ncol, nrow])    
    
    yflux1 = []
    yfluxdum1 =  np.zeros([ncol, nrow])    
    f = open("{0}/output/h.out".format(path), 'r')
    f.next()
    for line in f:
        a = (line.strip().split(" "))
        a = [myfloat(b) for b in a if b]
        try:
            j = int(a[0])-1
            k = int(a[1])-1
            hdum[j, k] = a[2]
            udum[j, k] = a[3]
            vdum[j, k] = a[4]
            infldum[j,k] = a[5]
            xfluxdum0[j,k] = a[6]
            yfluxdum0[j,k] = a[7] 
            xfluxdum1[j,k] = a[8]
            yfluxdum1[j,k] = a[9]                        
        except IndexError:
            dumt = int(a[0])
            h.append(hdum.copy())    
            u.append(udum.copy())
            v.append(vdum.copy())
            zinflmap2.append(infldum.copy())
            xflux0.append(xfluxdum0.copy())            
            yflux0.append(yfluxdum0.copy()) 

            xflux1.append(xfluxdum1.copy())            
            yflux1.append(yfluxdum1.copy())                       
    h = np.array(h)
    u = np.array(u)
    v = np.array(v)
    zinflmap2 = - np.array(zinflmap2)    
    xflux0 = np.array(xflux0)
    yflux0 = np.array(yflux0)  
    xflux1 = np.array(xflux1)
    yflux1 = np.array(yflux1)        
    # convert fortran zinflmap2 units of m3 to cm*m2
    inflVmap = zinflmap2*100 
    
    # convert fortran  m3 to cm*m2
    xflux0 = xflux0*100  # m3 --> cm*m2
    xflux1 = xflux1*100  # m3 --> cm*m2
    yflux0 = yflux0*100  # m3 --> cm*m2
    yflux1 = yflux1*100  # m3 --> cm*m2          
    
    return h,u,v,inflVmap,xflux0,yflux0,xflux1,yflux1

        
def read_dvol(path, Lx, Ly):
    """
    Reads fortran timeseries related to mass balance, saved with resolution dt_p. 
    The output is normalized by hillslope dimension, converting the units from m3 to cm.
    
    Returns:
        dvol (cm): change in surface volume 
        flux (cm): sum of lateral boundary fluxes 
        infl (cm): total infiltration
    """

    ta = []
    dvol = []
    infl = [] 
    flux = []

    f = open('{0}/output/dvol.out'.format(path), 'r'); 
    f.next()
    for line in f:
        a = (line.strip().split(" "))
        a = [myfloat(b) for b in a if b]
        ta.append(a[0])
        dvol.append(a[1])
        flux.append(a[2])
        infl.append(a[3])

    ta = np.array(ta)
    dvol = np.array(dvol)
    flux = np.array(flux)
    infl = np.array(infl)  

    dvol = dvol/Lx/Ly*100
    flux = flux/Lx/Ly*100
    infl = infl/Lx/Ly*100
    return dvol, flux, infl

def lateral_fluxes(path, params):
    """ 
    Read the boundary fluxes, and convert units from m^3 to cm/s 
    (normalizing by  the hillslope area and time step )

    """
    Ly =  params['nrow']*params['dx']
    Lx =  params['ncol']*params['dx']

    dt_p  = params['dt_p' ]
    flux1 = []
    flux2 = []
    flux3 = [] 
    flux4 = []
    f = open('{0}/output/fluxes1234.out'.format(path), 'r'); 
    f.next()
    for line in f:
        a = (line.strip().split(" "))
        a = [myfloat(b) for b in a if b]
        flux1.append(a[0])
        flux2.append(a[1])
        flux3.append(a[2])
        flux4.append(a[3])

    # convert units from m^3 to cm/s 
    flux1 = - np.array(flux1)/Lx/Ly*100/dt_p
    flux2 = - np.array(flux2)/Lx/Ly*100/dt_p
    flux3 = - np.array(flux3)/Lx/Ly*100/dt_p
    flux4 = - np.array(flux4)/Lx/Ly*100/dt_p
    
    return flux1, flux2, flux3, flux4
      
