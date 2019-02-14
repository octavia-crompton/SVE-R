#==============================================================================
# Main file for reading the fortran output and saving to sim.pklz
#==============================================================================
import pickle
import os, sys
from os.path import dirname
import multiprocessing as mp
import numpy as np
import json
import datetime
import shutil
import zipfile
import time
from glob import glob
import gzip

current_dir = os.getcwd()


mymodules = ['output_dry']
for mymod in mymodules:
    if mymod in sys.modules:
        del sys.modules[mymod]

from output_dry import *


def main(argv):

  sim_read()
  

def sim_read(path):
    """ 
    Reads the fortran output files and saves a pandas series to sim.pklz.
    sim.pklz includes the overland flow fields  (hc,uc,vc), infiltration (inflVmap),
    hydrograph, lateral fluxes and a mass balance check.
    """
    # load modified parameter file
    fname = '{0}/params_modified.json'.format(path)
    with open(fname) as f:
        params = json.load(f)
    
    for key,val in params.items():
            exec(key + '=val')              

    # list of variables that will be saved to sim.pklz
    fortan_outvars = ['t_p', 'cfl_p', 
                      't_h', 'hydro',
                      'runtime', 't_final', 't_pond',
                      'early_exit', 'no_flow',
                      'hmax', 'vmax', 'qmax',
                      'hc', 'uc', 'vc', 'inflVmap',
                      'dvol', 'flux', 'infl',
                      'volD', 'fluxD', 'rainD', 'inflD',
                      'flux1', 'flux2', 'flux3', 'flux4',
                      'mass_bal', 'infl_frac',   'zinflc', 
                      'xflux0','yflux0','xflux1','yflux1',
                      'vegTheta','vegH','bareTheta','bareH']     

    # read fortran output files
    fortan_outvars +=  params.keys()
                             
    t_p, cfl_p = read_time(path)        
    
    # read hydrograph, units = cm/s  (m3/s divided by hillslope area)
    t_h, hydro = read_hydro(path, Lx=Lx, Ly=Ly)  
    
    runtime, t_final, t_pond = read_tp(path)  
    
    early_exit = False  
    if  t_p[-1]/60. < tr:
      early_exit = True
    
    hc,uc,vc,inflVmap,xflux0,yflux0,xflux1,yflux1 = \
            read_h(path, ncol=ncol, nrow=nrow, dx= dx) 
        
    vmax = np.max(np.sqrt(uc**2 + vc**2), 0)*100    
    hmax = np.max(hc, 0)*100    
    qmax = hmax*vmax
    
    if np.max(np.abs(vc)) < 1e-3:
      no_flow = True
    else:
      no_flow = False
    
    #  mass balance file:  dvol.out
    dvol, flux, infl  = read_dvol(path, Lx = Lx, Ly = Lx) # cm                    
    
    # boundary fluxes. units: cm/s
    flux1, flux2, flux3, flux4 = lateral_fluxes(path, params) 
    
    # sample soil profiles
    vegTheta,vegH,bareTheta,bareH = read_ptsTheta(path , nz = int(zmax/dz)+ 1)      
      
    coord_path = '{0}/input/coords.pklz'.format(path)
    ff = gzip.open(coord_path ,'rb')
    coords = pickle.load(ff)
    ff.close()
    
    veg_path = '{0}/input/veg.pklz'.format(path)   
    ff = gzip.open(veg_path,'rb')
    veg = pickle.load(ff)
    ff.close()
    isvegc = veg['isveg'][:-1, :-1]  # trim because veg is defined at the nodes
    
    # infilration map (cm)
    zinflc = np.nansum(inflVmap, 0)/dx**2        
    
    # check mass conservation
    volD = np.nanmean(hc[-1])*100     # D for depth in cm
    fluxD = -np.nansum(flux)          #  units = cm
    rainD = rain*t_rain*100.            # units = cm
    inflD = - np.nansum(infl)         # units = cm, sim.zinflc.mean()
    mass_bal = rainD - volD - fluxD - inflD 
    
    if   rainD > 0:
        infl_frac = inflD/rainD  # infiltration fraction
    else:
        infl_frac = 0
        # need to compute volume from flux3s    
    
    maxQ= np.max(hydro)*3.6e4 # mm/hr
    maxV =  np.max(vmax) # cm/hr

        
    Q_eq = p - Ks*isvegc.mean() -  KsB*(1-isvegc.mean())  # Q = P - Ks, cm/hr
    
    ## save it again
    sim_dict = dict([(name,locals()[name]) for name in fortan_outvars])
    sim_dict.update(coords)  # update sim_dict with coords
    sim_dict.update({'isvegc' : isvegc})  # update sim_dict with veg
    
    fname = '{0}/sim.pklz'.format(path)
    ff = gzip.open(fname, 'wb')
    pickle.dump(sim_dict,ff)
    ff.close()  

  
if __name__ == '__main__':
    main(sys.argv)

