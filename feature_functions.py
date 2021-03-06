import numpy as np
import pandas as pd
import scipy as sp
import scipy.ndimage
from scipy import ndimage
np.seterr(divide='ignore', invalid='ignore')

def RF_patterns(isveg, rvl_params):        
    """
    inputs:
      isveg: binary array of vegetation 
      rvl_params: dictionary of parameters specifying how features are computed:
         edge : "padding" added to all distances to boundaries, to reduce edge effects.
         saturate : maximum size of elements in the feature arrays 
         gsigma : smoothing lengthscales (in the Gaussian filter)

    output: 
      pattern_dict: a dictionary of feature maps, each with shape (ncol x nrow)

     does not include: d2divide (ncol x nrow), the distance to divide
    """
    isvegc = np.array(isveg, dtype = float) 
    ncol = isvegc.shape[0]
    nrow = isvegc.shape[1]

    edge = int(rvl_params['edge'])    
    saturate = int(rvl_params['saturate'])        
    gsigma = rvl_params['gsigma']
      
    if gsigma == 0:
        gsigma = []
    elif type(gsigma) == int:
        gsigma = [gsigma]
        
    d2uB = func_d2uB(isvegc, edge, saturate)
    d2dB = func_d2dB(isvegc, edge, saturate)     
    
    d2lB = func_d2lB(isvegc, edge, saturate)
    d2rB = func_d2rB(isvegc, edge, saturate)    
    
    d2xB = np.minimum(d2lB, d2rB)    
    d2xB[d2xB > saturate] = saturate
                  
    d2uV = func_d2uV(isvegc, edge, saturate)
    d2dV = func_d2dV(isvegc, edge, saturate)

    d2lV = func_d2lV(isvegc, edge, saturate)
    d2rV = func_d2rV(isvegc, edge, saturate)   
    d2xV = np.minimum(d2lV, d2rV)
                  
    patchL,patchLB = get_patchL(isvegc, saturate) 
    bareL, bareLV = get_bareL(isvegc, saturate) 

    # assemble base pattern_dict
    pattern_dict = {'isvegc' : isvegc,                                                
                    'd2uB' : d2uB, 
                    'd2dB' : d2dB, 
                    'd2xB' : d2xB,
                    'd2uV' : d2uV, 
                    'd2dV' : d2dV, 
                    'd2xV' : d2xV,
                    'patchLB' : patchLB,
                    'bareLV' : bareLV,
                  }

    pos_sigma = [gs for gs in gsigma if gs > 0]
    
    for gs in pos_sigma:           
        for key in ['d2uV','d2dV','d2xV','bareLV']:  
            # smooth bare features
            pattern_dict[key + '_s{0}'.format(gs)] =   smoothB(pattern_dict[key], isvegc, gs)        
        
        for key in ['d2uB', 'd2dB', 'd2xB','patchLB']:
            # smooth vegetion features
            pattern_dict[key + '_s{0}'.format(gs)] =   smoothV(pattern_dict[key], isvegc, gs)        
    
    # delete unsmoothed features if 0 not in gsigma
    if 0 not in gsigma:
        
        for key in ['d2uV','d2dV','d2xV','bareLV',  \
                    'd2uB', 'd2dB', 'd2xB','patchLB']: 
            del pattern_dict[key] 
        
    return pattern_dict                

def get_feature_matrix(sim, rvl_params, target = None):
    """
    Compute feature matrix for an input vegetation field (sim['isvegc'])
 
    input: 
         sim: dictionary containing vegetation map (isvegc) and grid resolution (dx)
         
     veg features, bare features
    
    tasks:
      get pattern_dict from RF_patterns  (in ravel_fxns_RF.py)
      
      for all the keys in features:
        if in pattern_dict, np.ravel and add to pattern_ravel
        if in sim dict (non-array), add as list with equiv size
        skip 'local_path' key
        add 'd2divide' manually (units : grid cell)
    
    """
    isvegc = sim['isvegc'].astype(float)
    dx = sim['dx']
    
    ncol = isvegc.shape[0]
    nrow = isvegc.shape[1]
    
    xc = np.arange(0, ncol*dx, dx)  + dx/2
    yc = np.arange(0, nrow*dx, dx)  + dx/2
    xc, yc = np.meshgrid(xc, yc)
    
    xc = xc.T
    yc = yc.T
            
    pattern_dict = RF_patterns(isvegc, rvl_params)    
    features =  pattern_dict.keys() + ["d2divide"]
     
    pattern_ravel = {}
    eqv_size =  int(np.size(isvegc.ravel()))

    for key in features:
        if key in pattern_dict:
            pattern_ravel[key] = np.ravel(pattern_dict[key])

        elif key in sim.keys() and type(sim[key]) != np.ndarray:            
            pattern_ravel[key] = np.array([sim[key]]*eqv_size)  

        elif key =='d2divide':
            pattern_ravel[key] = (nrow - yc/dx).ravel() # ( nrow - yc/sim.dx).ravel()  

        else:
            print key
    
    pattern_ravel['isvegc'] = np.ravel(isvegc)
    pattern_ravel['fV'] = np.array([np.mean(isvegc)]*eqv_size)
    
    if target:
        pattern_ravel[target] = np.ravel(sim[target])        
 
    pattern_ravel = pd.DataFrame(pattern_ravel)

    return pattern_ravel
    
def smoothB(U, isvegc, gsigma):
    """
    Smooths array U over bare soil cells with a Gaussian filter,
      ignoring vegetated cells 

    Inputs: 
        U : array to smooth over
        isvegc : vegetated cells
        gsigma : standard deviation of the Gaussian kernel, determining length scale of smoothing
    Outputs:
        Z : array containing U smoothed over the bare soil cells (ignoring vegetated)
    """   
    U = U.astype(float)
    U[isvegc == 1] = np.nan
    V=U.copy()
    V[U!=U]=0
    VV=sp.ndimage.gaussian_filter(V,gsigma)

    W=0*U.copy()+1
    W[U!=U]=0
    WW=sp.ndimage.gaussian_filter(W,gsigma)

    Z=VV/WW
    Z = Z.astype(int)
    Z[isvegc ==1] = 0
    return Z


def smoothV(U, isvegc, gsigma):
    """
    Smooths array U over vegetated cells with a Gaussian filter,
      ignoring bare soil cells 

    Inputs: 
        U : array to smooth over
        isvegc : vegetated cells
        gsigma : standard deviation of the Gaussian kernel, determining length scale of smoothing
    Outputs:
        Z : array containing U smoothed over the vegetated cells
    """
    U = U.astype(float)
    U[isvegc == 0] = np.nan
    V=U.copy()
    V[U!=U]=0
    VV=sp.ndimage.gaussian_filter(V,gsigma)

    W=0*U.copy()+1
    W[U!=U]=0
    WW=sp.ndimage.gaussian_filter(W,gsigma)

    Z=VV/WW
    Z = Z.astype(int)
    Z[isvegc ==0] = 0
    return Z


def func_d2uB(isvegc, edge, saturate):
    """
    Computes the distance to nearest upslope bare soil cell
    Input: 
        isvegc: binary array of vegetation locations 
    Returns:
        df1 : array containing the distance to the nearest upslope bare cell
             =  0 for bare cells
             =  1 for veg cells with a neighboring bare cell upslope
             >  1 for veg cells with a bare cell further upslope
    """
    ncol = isvegc.shape[0]
    nrow = isvegc.shape[1]
    arr = np.ones([ncol+ 2*edge, nrow + 2*edge]).T
    arr[edge:-edge, edge:-edge] = np.flipud(isvegc.T)
    a = pd.DataFrame(arr) != 0

    df1 = a.cumsum()-a.cumsum().where(~a).ffill().fillna(0).astype(int)
    df1 = np.array(df1)[edge:-edge, edge:-edge]
    df1 = np.flipud(df1).T
    df1[isvegc == 0] = 0
   
    df1[df1>saturate] = saturate
    
    return df1
      

def func_d2uV(isvegc, edge, saturate):
    """
    Computes the distance to nearest upslope vegetated cell
    Input: 
        isvegc: binary array of vegetation locations 
    Returns:
        df1 : array containing the distance to the nearest upslope vegetated cell
             =  0 for vegetated cells
             =  1 for bare cells with a neighboring veg cell upslope
             >  1 for bare cells with a veg cell further upslope
    """

    ncol = isvegc.shape[0]
    nrow = isvegc.shape[1]
    arr = np.ones([ncol+ 2*edge, nrow + 2*edge]).T
    arr[edge:-edge, edge:-edge] = 1 - np.flipud(isvegc.T)
    a = pd.DataFrame(arr) != 0

    df1 = a.cumsum()-a.cumsum().where(~a).ffill().fillna(0).astype(int)
    df1 = np.array(df1)[edge:-edge, edge:-edge]
    df1 = np.flipud(df1).T
    df1[isvegc == 1] = 0

    df1[df1>saturate] = saturate  
    return df1
    

def func_d2dB(isvegc, edge, saturate):
    """
    Computes the distance to nearest downslope bare cell
    Input: 
        isvegc: binary array of vegetation locations 
    Returns:
        df1 : array containing the distance to the nearest downslope bare cell
    """
    ncol = isvegc.shape[0]
    nrow = isvegc.shape[1]
    
    arr = np.ones([ncol+ 2*edge, nrow + 2*edge]).T
    arr[edge:-edge, edge:-edge] = isvegc.T
    a = pd.DataFrame(arr) != 0

    df1 = a.cumsum()-a.cumsum().where(~a).ffill().fillna(0).astype(int)
    df1 = np.array(df1)[edge:-edge, edge:-edge]
    df1 = df1.T
    df1[isvegc == 0] = 0    
     
    df1[df1>saturate] = saturate
 
    return df1

def func_d2dV(isvegc, edge, saturate):
    """
    Computes the distance to nearest downslope vegetated cell
    Input: 
        isvegc: binary array of vegetation locations 
    Returns:
        df1 : array containing the distance to the nearest downslope vegetated cell
    """
  

    ncol = isvegc.shape[0]
    nrow = isvegc.shape[1]
    arr = np.ones([ncol+ 2*edge, nrow + 2*edge]).T
    arr[edge:-edge, edge:-edge] = 1 - isvegc.T
    a = pd.DataFrame(arr) != 0

    df1 = a.cumsum()-a.cumsum().where(~a).ffill().fillna(0).astype(int)
    df1 = np.array(df1)[edge:-edge, edge:-edge]
    df1 =  df1.T
    df1[isvegc == 1] = 0
   
    df1[df1>saturate] = saturate
    
    return df1



def func_d2lB(isvegc, edge, saturate):
    """
    input:
      isvegc : [ncol x nrow] array of vegetation field

    output :
      d2lB : [ncol x nrow] array of distance to nearest left bare
        =  0 for bare cells
        =  1 for veg cells with a bare cell immediately left
    """

    ncol = isvegc.shape[0]
    nrow = isvegc.shape[1]
    arr = np.ones([ncol+ 2*edge, nrow + 2*edge])
    arr[edge:-edge, edge:-edge] = isvegc
    a = pd.DataFrame(arr) != 0

    df1 = a.cumsum()-a.cumsum().where(~a).ffill().fillna(0).astype(int)
    df1 = np.array(df1)[edge:-edge, edge:-edge]
    df1[isvegc == 0] = 0

    df1[df1>saturate] = saturate
    
    return df1


        
def func_d2lV(isvegc, edge, saturate):
    """
    input: 
      isvegc : [ncol x nrow] array of vegetation field

    output : 
      d2lV : [ncol x nrow] array, distance to nearest veg cell to the left
        =  0 for veg cells
        =  1 for bare cells with a veg cell immediately left  
    """

    ncol = isvegc.shape[0]
    nrow = isvegc.shape[1]
    arr = np.ones([ncol+ 2*edge, nrow + 2*edge])
    arr[edge:-edge, edge:-edge] = 1 - isvegc
    a = pd.DataFrame(arr) != 0

    df1 = a.cumsum()-a.cumsum().where(~a).ffill().fillna(0).astype(int)
    df1 = np.array(df1)[edge:-edge, edge:-edge]
    df1[isvegc == 1] = 0

    df1[df1>saturate] = saturate

    return df1


def func_d2rB(isvegc, edge, saturate):
    """
    input:
      isvegc : [ncol x nrow] array of vegetation field

    output :
      d2rB : [ncol x nrow] array;  distance to nearest bare cell to right
        =  0 for bare cells
        =  1 for veg cells with a veg cell immediately to right
    """

    ncol = isvegc.shape[0]
    nrow = isvegc.shape[1]
    arr = np.ones([ncol+ 2*edge, nrow + 2*edge])
    arr[edge:-edge, edge:-edge] = np.flipud(isvegc)
    a = pd.DataFrame(arr) != 0

    df1 = a.cumsum()-a.cumsum().where(~a).ffill().fillna(0).astype(int)
    df1 = np.array(df1)[edge:-edge, edge:-edge]
    df1 = np.flipud(df1)
    df1[isvegc == 0] = 0
      
    df1[df1>saturate] = saturate
    
    return df1
    
def func_d2rV(isvegc, edge, saturate):
    """
    input: 
      isvegc : [ncol x nrow] array of vegetation field

    output : 
      d2rV : [ncol x nrow] array;  distance to nearest veg cell to right
        =  0 for veg cells
        =  1 for bare cells with a veg cell immediately to right  
    """
    ncol = isvegc.shape[0]
    nrow = isvegc.shape[1]
    arr = np.ones([ncol+ 2*edge, nrow + 2*edge])
    arr[edge:-edge, edge:-edge] = 1- np.flipud(isvegc)
    a = pd.DataFrame(arr) != 0

    df1 = a.cumsum()-a.cumsum().where(~a).ffill().fillna(0).astype(int)
    df1 = np.array(df1)[edge:-edge, edge:-edge]
    df1 = np.flipud(df1)
    df1[isvegc == 1] = 0

    df1[df1>saturate] = saturate

    return df1
 

def get_patchL(isvegc, saturate):
    """
    input : isvegc from get_source(df)  
    
    output : 
  
      patchLv:  vegetated patch length
      patchLB:  upslope interspace patch length (paired to veg patch)
      patchLc:  charcteristic length  Lv/(Lv + Lb)
      
      Ldict:  dictionary of veg patch lengths. 
         Ldict key :  downslope patch coordinate 
      Bdict:  dictionary of paired upslope interspace lengths. 
  
    usage: 
      patchLv,patchLB,patchLc,Ldict,Bdict = get_patchL(isvegc)
    """
    ncol = isvegc.shape[0]
    nrow = isvegc.shape[1]
    patchLv = np.zeros(isvegc.shape, dtype = float)  # veg patch length
    patchLB = np.zeros(isvegc.shape, dtype = float)  # upslope interspace patch length (paired to veg patch)
    
    for i in range(ncol):  # loop over across-slope direction first
        count = 0           
        for j in range(nrow):    
            if isvegc[i, j] == 1:    #  if veg patch, add 1
                if j >= (nrow -1):  # if we're at the top of the hill                  
                  patchLv[i, j-count:] = count  # record veg patch length                  
                count += 1  
                                                        
            # if [i,j] is bare and the slope cell is vegetated, record.
            # each patch starts at [i,j-count] and ends at [i,j-1]
            elif isvegc[i,j] == 0 and isvegc[i, j-1] == 1:   
                if j > 0:
                  # veg patch starts at j-count and ends at j
                  patchLv[i, j-count:j] = count
                  try:
                      # find the nearest upslope veg cell
                      Lb = np.where(isvegc[i,j:] == 1)[0][0]                               
                      patchLB[i,j-count:j] = Lb
                  except IndexError:  # bare patch extends to top of hill
                      patchLB[i,j-count:j] = nrow - j
                  count = 0 
    patchLv[patchLv > saturate] = saturate
    patchLB[patchLB > saturate] = saturate
        
    return  patchLv, patchLB


def get_bareL(isvegc, saturate, skipflag = 0):
    """
    input : isvegc from get_source(df)  
    
    output : 
  
      bareL:  bare patch length
      bareLV:  upslope vegeted patch length (paired to bare patch)
      
    """
    ncol = isvegc.shape[0]
    nrow = isvegc.shape[1]
    bareLV = np.zeros(isvegc.shape, dtype = float)  # veg patch length
    bareL = np.zeros(isvegc.shape, dtype = float)  # upslope interspace patch length (paired to veg patch)

    for i in range(ncol):  # loop over across-slope direction first
        count = 0           
        for j in range(nrow):    
            if isvegc[i, j] == 0:    #  if bare, add 1
                if j >= (nrow -1):  # if we're at the top of the hill                  
                  bareL[i, j-count:] = count  # record bare length                  
                count += 1  
 
            # if [i,j] is veg and the slope cell is bare, record.
            # each patch starts at [i,j-count] and ends at [i,j-1]
            elif isvegc[i,j] == 1 and isvegc[i, j-1] == 0:   
                if j > 0:
                  # veg patch starts at j-count and ends at j
                  bareL[i, j-count:j] = count
                  if skipflag == 1 and j-count ==0:
                      bareL[i, j-count:j] = 0
                  try:
                      # find the nearest upslope bare cell
                      Lb = np.where(isvegc[i,j:] == 0)[0][0]                               
                      bareLV[i,j-count:j] = Lb
                  except IndexError:  # bare patch extends to top of hill
                      bareLV[i,j-count:j] = nrow - j
                  count = 0 
    bareLV[bareLV > saturate] = saturate
    bareL[bareL > saturate] = saturate

    return  bareL, bareLV
