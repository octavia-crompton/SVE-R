import os
import sys
from os.path import dirname
import pickle
import shutil
import json

mymodules = ['input_boundary', 
             'input_coords',
             'input_phi',              
             'input_param',  
             'input_veg',                         
            ]
             
for mymod in mymodules:
    if mymod in sys.modules:  
        del sys.modules[mymod]
               

from input_coords import *
from input_boundary import *
from input_veg import *
from input_phi import *

   
def main(argv):
    """          
    writes the input files for dry.f:       
        coords.dat
        params.dat
        vanG.dat
        boundary.dat
    copies source code (dry.for, dry.inc) to simulation directory     
    """

    sim_input(sim_dir) 
      
    
def sim_input(path):
    """
    tasks:  
      write the Fortran input files.
        includes: 
          coords.dat  
          nodes.dat
          veg.dat
          boundary.dat
          vanG.dat
      
    copy source code (dry.for, dry.inc)    
    
    calls functions:
      build_coords: save coords to path + '/input/coords.dat'  
    
    """ 

    fname = '{0}/params.json'.format(path)
    params = json.load(open(fname))

    for key,val in params.items():
            exec(key + '=val')

    # make input andoutput subfolders
    input_path = '/'.join([path, 'input'])
    os.system('mkdir {0}'.format(input_path))  if os.path.isdir(input_path) == False else True

    output_path = '/'.join([path, 'output'])
    os.system('mkdir {0}'.format(output_path))  if os.path.isdir(output_path) == False else True

    #
    if vegtype == 'randv':
        
        isveg = wrap_veg(path + '/input', params)                          

    elif vegtype == 'stripes':

      isveg = wrap_stripe_veg(path + '/input', params)

    elif vegtype == 'image':
                 
      isveg = wrap_image_veg(path , params)
  
    else:
        isveg = wrap_veg(path + '/input', ncol, nrow, 0, fV, 0)
    
    
    params['Ly'] =  params['nrow']*params['dx']
    params['Lx'] =  params['ncol']*params['dx']
    
    track_points = json.load(open('{0}/track_points.json'.format(input_path)))
    params.update(track_points)    

    try:  
      plot_veg(isveg, path)
    except Exception as e:
      print e    

    wrap_coords(input_path, params)
    nop = write_nodes(input_path, ncol, nrow)  

    write_boundary(input_path, params)
    # sae track points
    for key,val in track_points.items():
            exec(key + '=val') 
    
    params['t_rain'] =  params['tr']*60    # storm duration in seconds
    params['rain'] =  params['p']/3.6e5    # in m/s
    if 'dt_p' not in params:
         params['dt_p'] = int(params['t_rain']/params['print_count'])
    
    if 'tmax_scale' in params:
        params['tmax'] =  params['t_rain']*params['tmax_scale']
    else:
        params['tmax'] =  params['t_rain']*4

    
    nprt = int(np.maximum( params['dt_p']/dt_sw, 1))  
    nt_sw = int(params['tmax']/dt_sw)
    params['nt'] = nt_sw
    params['nprt'] = nprt    
    

    params['iscale']  = int(params['dt_r']/params['dt_sw'])
    
    # roughness parameters are hard-coded to Manning's equation
    params['m'] = 2./3         
    params['alphaB'] = 0.03
    params['mB'] = 2/3.

    params['eta'] = 1/2. if params['m'] != 2 else 1
    params['etaB'] = 1/2. if params['mB'] != 2 else 1
          
    
    soil_params =  { "Ks" : Ks,
                    "L": 0.5,
                    "n": 1.47,
                    "theta_r": 0.0378,
                    "theta_s": 0.472,
                    "alpha": 0.0096,
                    "lambda": 0.47}
                    
    params['vG'] = soil_params
    
    write_param(path,  params)


    with open('{0}/params_modified.json'.format(path), 'w') as file:
      file.write(json.dumps(params)) 
            
    phi_veg, phi_bare, nz  = make_phi(params)
    

    H_init =  np.flip(np.arange( H_i, H_i + zmax+dz, dz), 0)
    H_init[H_init > -1] = -1  #  H_init cannot exceed -

    write_phi(input_path,phi_veg,phi_bare, H_init, dz, zmax)          

    # modify dry.inc and copy it to path as dry.inc
    write_inc(path, params)  

    # copy the model to the new folder
    shutil.copyfile('dry.for', '{0}/dry.for'.format(path))

def write_param(path, params):
    
    for key,val in params.items():
            exec(key + '=val')
    
    fname = '{0}/input/params.dat'.format(path)    
    f = open(fname, 'w')
    f.write('grav        dt           \n')
    f.write('9.806d0     {0}          \n'.format(dt_sw))
    f.write('tmax      t_rain    \n')
    f.write('{0}       {1}     \n'.format(tmax,t_rain)) 
    f.write('prate      nt  \n')
    f.write('{0}        {1}   \n'.format(rain, nt))        
    f.write('epsh      beta     \n')  
    f.write('{0}        1.0     \n'.format(epsh))
    f.write('nprt      \n')
    f.write('{0}       \n'.format(int(nprt)))   
    f.write('iscale   \n')
    f.write('{0}       \n'.format( iscale))  
    f.write('stop_tol   \n')
    f.write('{0}       \n'.format( stop_tol))      
    f.write('h0      u0    v0   \n ')
    f.write('0.0   0.0    0.0  \n ')
    f.write('mV      etaV     alphaV   \n')
    f.write('{0}     {1}      {2} \n'.format(m, eta, alpha))  
    f.write('mB      etaB     alphaB   \n')
    f.write('{0}     {1}      {2} \n'.format(mB, etaB, alphaB))        
    f.write('jveg     kveg       \n')
    f.write('{0}      {1} \n'.format( jveg+1, kveg+1))  
    f.write('jbare     kbare       \n')
    f.write('{0}       {1} \n'.format( jbare+1, kbare+1)) 
    f.write('tsat_min \n')
    f.write('{0} \n'.format(tsat_min))
    f.close()  

def write_inc(path, params):
    #  Avoid allocating to much space   
    nt = params['nt']
    zmax = params['zmax']
    dz = params['dz']    
    nz =  int(zmax/dz)+ 1
    ncol = params['ncol']    
    nrow = params['nrow']    
    with open('dry.inc', 'r') as input_file: 
      with open('{0}/dry.inc'.format(path), 'w') as output_file:
        for line in input_file:
            if line[6:15] == 'parameter':
                # print 'old line:', line
                a = (line.strip().split(" ")[2].split(","))
                nn0 = int(a[0].split('=')[-1])
                nt0 = int(a[1].split('=')[-1])
                nx0 = int(a[2].split('=')[-1])
                ny0 = int(a[3].split('=')[-1])
                nz0 = int(a[4].split('=')[-1])

                newline = '      parameter ( nn={0},ntp={1},nx={2},ny={3},nz={4} )'.format(
                     (ncol+1)*(nrow+1)+1, nt+1, ncol+2, nrow+2, nz)

                output_file.write(newline)
            else:
                output_file.write(line)
     
        
if __name__ == '__main__':
    main(sys.argv)
    
