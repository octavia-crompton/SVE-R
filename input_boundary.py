import numpy as np
import scipy as sp
import os, sys


def write_boundary(path, params):
    """
    Writes  boundary file, boundary.dat, for a rectangular domain where
    all boundary cells of a given orientation (i.e. south, north, east,or west) 
    are of the same type. 
    name conventions for the orientations are in grid/boundaries.png 
        ipos = 3 for the divide / top of hill
        ipos = 1 for the channel / bottom of hill.
    """
    ncol = params['ncol']
    nrow = params['nrow']
    itype3 = params['itype3']
    itype1 = params['itype1']
    itype2 = params['itype2']
    itype4 =  params['itype4']
                   
    inum = np.zeros([ncol+1, nrow+1], dtype = int)
    inum[1:, 1] = 1
    inum[1:, -1]= 1
    inum[1, 1:] = 1
    inum[-1, 1:] = 1
    inum[1, 1] = 2
    inum[1, -1] = 2
    inum[-1, -1] = 2
    inum[-1, 1] = 2
    
    ipos = np.zeros( [ncol+1, nrow+1, 2], dtype = int)
    
    # bottom boundary
    ipos[2:-1, 1,0] = 1
    ipos[1, 1,1] = 1
    ipos[-1, 1,1] = 1

    # right boundary
    ipos[-1, 1:-1, 0] = 2
    ipos[-1, -1,1] = 2

    # left boundary
    ipos[1, 1:, 0] = 4

    # top boundary
    ipos[2:, -1,0] = 3
    ipos[1, -1,1] = 3
    
    itype = np.zeros([ncol+1, nrow+1, 2], dtype = int)
    # bottom boundary
    itype[2:-1, 1,0] = itype1
    itype[1, 1,1] = itype1
    itype[-1, 1,1] = itype1

    # right boundary
    itype[-1, 1:-1, 0] = itype2
    itype[-1, -1,1] = itype2

    # left boundary
    itype[1, 1:,0] = itype4

    # top boundary
    itype[2:, -1,0] = itype3
    itype[1, -1,1] = itype3
    
    npt = (ncol+1)*(nrow+1)  # number of points
    nbcell = 2*ncol + 2*nrow - 4  # number of boundary cells
    
    fname = '{0}/boundary.dat'.format(path)    
    f = open(fname, 'w')
    f.write('number of boundary cell \n') 
    f.write('  {0} \n'.format(nbcell))
    f.write(' j    k          inum    itype             ipos \n')
    # f.write(' j \t k \tinum    itype \t\t ipos')
    j = 1
    for k in range(1, nrow+1):
        if inum[j, k] == 2:
            f.write( '{0:<5} {1:<13} {2:<7} {3:<8} {4:<9} {5:<8} {6:<6} \n'.format(
                        j, k, inum[j, k], itype[j, k, 0], itype[j, k, 1], 
                         ipos[j, k, 0], ipos[j, k, 1]))
        else:
            f.write( '{0:<5} {1:<13} {2:<7} {3:<18} {4:<10}   \n'.format(
                         j, k, inum[j, k],  itype[j, k, 0],  ipos[j, k, 0], ))

    for j in range(2, ncol+1):
        if inum[j, k] == 2:
            f.write( '{0:<5} {1:<13} {2:<7} {3:<8} {4:<9} {5:<8} {6:<6} \n'.format(
                        j, k, inum[j, k], itype[j, k, 0], itype[j, k, 1], 
                         ipos[j, k, 0], ipos[j, k, 1]))
        else:
            f.write( '{0:<5} {1:<13} {2:<7} {3:<18} {4:<10}   \n'.format(
                         j, k, inum[j, k],  itype[j, k, 0],  ipos[j, k, 0], ))

    for k in range(nrow-1,0,-1):
        if inum[j, k] == 2:
            f.write( '{0:<5} {1:<13} {2:<7} {3:<8} {4:<9} {5:<8} {6:<6} \n'.format(
                        j, k, inum[j, k], itype[j, k, 0], itype[j, k, 1], 
                         ipos[j, k, 0], ipos[j, k, 1]))
        else:
            f.write( '{0:<5} {1:<13} {2:<7} {3:<18} {4:<10}   \n'.format(
                         j, k, inum[j, k],  itype[j, k, 0],  ipos[j, k, 0], ))
            
    for j in range(ncol-1,1,-1):
        if inum[j, k] == 2:
            f.write( '{0:<5} {1:<13} {2:<7} {3:<8} {4:<9} {5:<8} {6:<6} \n'.format(
                        j, k, inum[j, k], itype[j, k, 0], itype[j, k, 1], 
                         ipos[j, k, 0], ipos[j, k, 1]))
        else:
            f.write( '{0:<5} {1:<13} {2:<7} {3:<18} {4:<10}   \n'.format(
                         j, k, inum[j, k],  itype[j, k, 0],  ipos[j, k, 0], ))

    kbeg = np.ones(ncol+1, dtype = int)
    kend = np.ones(ncol+1, dtype = int)*nrow
   
    f.write('ncol\n')
    f.write("{0}\n".format(ncol))
    f.write('nrow\n')
    f.write("{0}\n".format(nrow))    
    f.write('j     kbeg          kend \n')
    for j in range(1, ncol+1):
        f.write( '{0:>5}  {1:>5} {2:>13}   \n'.format(
                    j, kbeg[j],kend[k] ))
        
    f.close()
    boundary_fix(path , params )
    return inum, ipos, itype

def boundary_fix(path, params):  
    """
    write fixed flux boundary condition at top of hill, if if params['itype3']  == 4
     
    influx : unit influx q in m^2/s
        influx can be related to the equivalent p (in cm/hr) as influx*3.6e5/Ly,
        or p*Ly/3.6e5 = equiv influx
    
        if 'influx' is not provided, specify as equivalent 5cm/hr

    """
    fname = '{0}/boundary.dat'.format(path)    
    f = open(fname, 'a')
    ncol = params['ncol'] 
    nrow = params['nrow']     
    if params['itype3']  == 4:
        fixj = np.arange(1, ncol+1)
        fixk = np.ones(ncol, dtype = int)*nrow
        fixh = np.zeros(ncol, dtype = float )
        fixu = np.zeros(ncol, dtype = float)
        influx = params['influx'] if 'influx' in params else 5*params['Ly']/3.6e5
        fixv = - np.ones(ncol, dtype = float)*influx 
        ndir = len(fixj)
    else: 
        ndir = 0
    
    f.write('number of fixed bc cells, ndir    \n ')
    f.write('{0}  \n '.format(ndir))   
    f.write('j     k    fix h    fix u    fix v	\n')
    for i in range(ndir): 
        f.write('{0}     {1}     {2}     {3}     {4}     \n'.format(
                fixj[i],fixk[i],fixh[i],fixu[i],fixv[i]))
    f.close()
