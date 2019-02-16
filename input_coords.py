import numpy as np
import scipy as sp
import os, sys
import gzip, pickle

def wrap_coords(path,params):
    """
    inputs: 
      path: path to save input files
      params : dictionary with parameters, including: 
          topo: topography case
          ncol: across slope number of cells
          nrow: along slope number of cells
          dx: grid cell width
          So : slope, m/m (not in percent)
    
    1. call build_coords  --> x,y,z values at nodes
    2. call write_coords --> write to coords.dat  
    3. call cc_coords --> save cell-center coordinates to coords.pklz (for python)
    """   
    
    ncol = params['ncol']
    nrow = params['nrow']    
    dx = params['dx']  
    
    nop = write_nodes(path, ncol, nrow)  
    x, y, z  = build_coords(params)
    write_coords(path, ncol, nrow, dx, x, y, z)
    cc_coords(path, ncol, nrow, nop, x, y, z, dx)
  
    return x,y,z

def cc_coords(path, ncol, nrow, nop, x, y, z, dx):
  """
  save cell center coordinates
    
  input:
    path, ncol, nrow
    x, y, z : node coordinates

  1. interpolate nodes to cell centers
  2. save a pickle in the path directory
  """
  xc = interp2nodes(ncol,nrow, nop, x)    
  yc = interp2nodes(ncol,nrow, nop, y)    
  zc = interp2nodes(ncol,nrow, nop, z)   
  
  d2divide = nrow - yc/dx # number of grid cells (not m)
  coord_dict = {'xc' : xc, 'yc': yc, 'zc' :zc, 'd2divide' : d2divide}

  fname = '{0}/coords.pklz'.format(path)
  f = gzip.open(fname, 'wb')
  pickle.dump(coord_dict,f)
  f.close()   
          

def build_coords(params):
    """
    creates x,y,z coordinates, called by wrap_coords
                
    input: 
      params : dictionary with parameters, including:
          ncol, nrow, dx, slope, 
          topo (topography case)
    
    returns: 
      xdum: [ncol+1, nrow+1] , x at nodes
      ydum: [ncol+1, nrow+1] , y at nodes
      zdum: [ncol+1, nrow+1] , z at nodes
    """
    for key,val in params.items():
            exec(key + '=val')
        
    if 'plane' in topo:
        # construct grid
        x = np.arange(0, (ncol+1)*dx - 1e-10, dx )
        y = np.arange(0, (nrow+1)*dx - 1e-10, dx )
        y, x = np.meshgrid(y, x)

        zymax = So*(np.max(y) - np.min(y))
        zrow =  np.linspace(0,zymax, nrow+1)
        z = np.tile(zrow, [ncol+1]).reshape([ncol+1, nrow+1])
    
    elif 'nonplane' in topo:
        # example non-planar topography
        omega = 0.0001
        Ly = nrow*dx
        H = Ly*So
        y = 0.82*np.arange(nrow+1)
        x =  (-20 + 4/3 * ( np.arange(ncol+1)  ))
        y,x = np.meshgrid(y, x)
        
        a = -2*omega/So
        x =  x*np.exp(a*y) + 20

        z = H*( y/Ly ) + 2*(omega)*( x - 20 )
        
    return x.ravel(), y.ravel(), z.ravel()          
   
   
def write_coords(path, ncol, nrow, dx, x, y, z):        
    
    """
     writes x,y,z coordinates to coords.dat
    """
    npt = (ncol+1)*(nrow+1)  # number of points
    ne = nrow*ncol           # number of edges
    
    fname = '{0}/coords.dat'.format(path)
    f = open(fname, 'w')
    f.write('{0:<13}   {1:<13}\n'.format(npt, ne))
    
    # write x, y, z
    for n in range(npt):
        f.write('{0:<13.6f} {1:<13.6f} {2:<13.6f}  \n'.format(
                    x[n],y[n],z[n]))
    f.close()
    
def write_nodes(path, ncol, nrow):
    """
    write cell node indices to nodes.dat
    """
    npt = (ncol+1)*(nrow+1)  # number of points    
    nodes = np.arange(1, npt+1, dtype = int).reshape([ncol+1, nrow+1])

    nop = np.zeros([ncol, nrow, 4], dtype = int)
    for j in range(ncol):
        for k in range(nrow):
            nop[j, k] =  nodes[j,k], nodes[j+1, k], nodes[j+1,k+1], nodes[j,k+1]

    fname = '{0}/nodes.dat'.format(path)
    f = open(fname, 'w')
    f.write('{0:<10} {1:<10}  {2:<10} {3:<10}\n'.format("n1", "n2", "n3", "n4")) 
    # write node numbers  
    for j in range(ncol):
        for k in range(nrow):
            n1 = nop[j, k, 0] 
            n2 = nop[j, k, 1]       
            n3 = nop[j, k, 2]        
            n4 = nop[j, k, 3] 
            f.write('{0:<10} {1:<10}  {2:<10} {3:<10}\n'.format(n1, n2, n3, n4)) 
    f.close() 

    return nop

def interp2nodes(ncol,nrow, nop, x):
    """ 
    interpolates values of array x from nodes to cell centers
    """
    xcc  = np.zeros([ncol, nrow])    

    for j in range(ncol):
        for k in range( nrow):
            n1 = nop[j, k, 0] - 1
            n2 = nop[j, k, 1] - 1    
            n3 = nop[j, k, 2] - 1            
            n4 = nop[j, k, 3] - 1    
            xcc[j,k] = 0.25*(x[n1] + x[n2] + x[n3] + x[n4])                            
                
    return xcc

            
if __name__ == '__main__':
     main(sys.argv)
 
