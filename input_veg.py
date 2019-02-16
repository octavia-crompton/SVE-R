#==============================================================================
# functions for creating vegetation pattern and writing input/veg.dat
#==============================================================================
import numpy as np
import scipy as sp
import os, sys
import json
import gzip, pickle
from scipy.ndimage.filters import gaussian_filter
from skimage import io, color
import imageio
from scipy.interpolate import RegularGridInterpolator    
import matplotlib.pylab as plt

def wrap_veg(path, params):
    """
    This function is called if vegtype==randv.  It calls make_randv() to 
    create a random vegetation field, write_veg() to write input/veg.dat, 
    and get_points() to select two points for which the soil moisture profiles will be saved.
    
    args:
      path: path to save veg.dat
      params: parameter dictionary
    """   
   
    ncol = params['ncol']
    nrow = params['nrow']    
    fV = params['fV']        
    sigma_x = params['sigma_x']    
    sigma_y = params['sigma_y']        
    seed = params['seed']

    isveg = make_randv(ncol, nrow, fV, sigma_x, sigma_y, seed)
    isveg = isveg.astype(int)
    
    # try to plot it (but skip if matplotlib throws an error)
    try:  
      plot_veg(isveg, path)
    except: 
      pass

    write_veg(path, ncol, nrow, isveg )
    get_points(path, isveg, nrow, ncol)   
    
    return isveg
    
def make_randv(ncol, nrow, fV, sigma_x, sigma_y, seed):
    """
    Makes a random vegetation field with specified fractional cover (fV) and 
    lengthscales (sigma_x, sigma_y)
    """
    np.random.seed(seed) 

    # makes a binary vegetation array with density fV
    isveg = sp.rand(ncol+1, nrow+1) >=  1 - fV      
    isveg = isveg.astype(float)

    # apply gaussian filter if sigma > 0
    if (sigma_x > 0 or sigma_y > 0) and fV < 1.0:
      blurred = gaussian_filter(isveg.astype(float), sigma=(sigma_x, sigma_y))
      isveg = (blurred> np.percentile(blurred, 100*(1-fV))).astype(int)

    isveg = isveg.astype(float)
    return  isveg


def wrap_image_veg(path, params):
    """
    Creates vegetation map from an image file, specified by "image_name" in params.json
    This function is only called if vegtype == image
    """
    filepath = '/'.join([path,  params['image_name']])

    trim = 1
    threshold = 0.5
        
    image = imageio.imread(filepath)
    image = color.rgb2gray(image[trim:-trim, trim:-trim])
    image = np.fliplr(image.T)

    image = 1.0*(image > float(threshold)) # binarize
    scale = image.shape[1]/1.0/params['nrow'] # might be better to get scale from Ly

    y = np.arange(0, image.shape[0])
    x = np.arange(0, image.shape[1])
    interpolating_function = RegularGridInterpolator((y, x), image)

    xv = np.arange(0, image.shape[0], scale)
    yv = np.arange(0, image.shape[1], scale)

    xv, yv = np.meshgrid(xv,yv)
    image2 = interpolating_function((xv, yv)).T   
    image2 = 1.0*(image2 > float(threshold)) # binarize

    ncol, nrow = params['ncol'], params['nrow']
    isveg = np.zeros([ncol+1, nrow+1])

    nx = min(ncol, image2.shape[0])
    ny = min(nrow, image2.shape[1])

    isveg[:nx, :ny] = image2[:nx, :ny]

    write_veg(path + '/input', ncol, nrow, isveg )
    get_points(path+  '/input', isveg, nrow, ncol) 
    
    return isveg

def wrap_stripe_veg(path, params):
    """
    Creates a vegetation field of horizonal stripes with randomly-selected
     widths and spacings, determined by fractional cover and sigma_x

    args: 
      path: path to save coords.dat
      ncol: across-slope number of cells
      nrow: along-slope number of cells
      seed : random seed
    """
    ncol = params['ncol']
    nrow = params['nrow']    
    fV = params['fV']        
    sigma = params['sigma_x']    
    seed = params['seed']
    
    np.random.seed(seed) 

    isveg = (sp.rand( nrow+1) >=  1 - fV ).astype(float)                
    
    if sigma > 0:
        blurred = gaussian_filter(isveg, sigma=sigma)
        isveg = (blurred> np.percentile(blurred, 100*(1-fV))).astype(int)
    isveg = np.tile(isveg, (ncol + 1, 1))
    
    try:  
      plot_veg(isveg, path)
    except: 
      pass

    write_veg(path, ncol, nrow, isveg )
    get_points(path, isveg, nrow, ncol)   
    return isveg

def plot_veg(isveg, path):
    """
    saves a png of the vegetation field 
    """
    fig = plt.figure(figsize = (4,5))

    plt.pcolormesh(isveg[:-1, :-1].T, cmap = 'Greens')  
    name = path.split('/')[-1]  
    plt.savefig('{0}/veg_{1}.png'.format(path, name),bbox_inches='tight')

    plt.close(fig)

   
def write_veg(path, ncol, nrow, isveg):        
    """
    writes the input vegetation field to input/veg.dat 
    
    args:
    isveg
    2D isveg --> 1D isveg = isveg.ravel()
    1D isveg --> 2D isveg = isveg.reshape([ncol+1, nrow+1])
        
        call write_veg -->
         write isveg (1D) to veg.dat 
         write isvegc ([ncol, nrow]) to veg.pklz   
    """
    npt = (ncol+1)*(nrow+1)  # number of points
    isveg = isveg.ravel()
    
    fname = '{0}/veg.dat'.format(path)
    f = open(fname, 'w')
    for n in range(npt):
      f.write('{0}  \n'.format(isveg[n])) 
    
    f.close()    
    
    isveg = isveg.reshape([ncol+1, nrow+1])#[:-1, :-1]
    
    fname = '{0}/veg.pklz'.format(path)
    f = gzip.open(fname, 'wb')
    pickle.dump({'isveg' : isveg},f)
    f.close()

def load_veg(path, filename):
  
  f = gzip.open('/'.join([path, filename]),'rb')
  isveg = pickle.load(f)
  f.close()
  
  return isveg
  
def get_points(path, isveg, nrow, ncol):
    """
    jveg, kveg : find vegetated cell nearest to the center of the bottom of the hillslope
    jbare, kbare : find bare cell nearest to the center of the bottom of the hillslope    
    
    note: jveg, kveg,  jbare, kbare are in python coordinates
    """
    try:
      jveg = np.where(isveg == 1)[0]
      kveg = np.where(isveg == 1)[1]
  
      vegind = np.where((abs(jveg - ncol/2)**2 +  kveg**2) == np.min(abs(jveg -ncol/2)**2 +  kveg**2))[0][0]
      jveg = jveg[vegind]
      kveg =  kveg[vegind]  

    except:
      jveg = 1
      kveg =  1
    
    try:      
      jbare = np.where(isveg == 0)[0]
      kbare = np.where(isveg == 0)[1]

      bareind = np.where((abs(jbare - ncol/2)**2 +  kbare**2) == np.min(abs(jbare-ncol/2)**2 +  kbare**2))[0][0]      

      jbare = jbare[bareind]
      kbare = kbare[bareind]
    except:
        
      jbare = 0
      kbare = 0
      
    track_points = {'jveg': jveg, 'kveg': kveg, 'jbare': jbare, 'kbare' : kbare}
    
    with open('{0}/track_points.json'.format(path), 'w') as file:
      file.write(json.dumps(track_points))  
