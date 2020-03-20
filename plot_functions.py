import numpy as np
import pickle
import matplotlib.pylab as plt
import matplotlib.gridspec as gridspec
from scipy import ndimage
from scipy.ndimage import convolve
import pandas as pd
import matplotlib.colors as colors
import cmocean

import warnings
warnings.filterwarnings("ignore")

plt.style.use('ggplot')
import seaborn as sns
sns.set(font_scale = 1.3)


plt.rcParams['axes.facecolor']='w'
plt.rcParams['grid.color']= 'grey'
plt.rcParams['grid.alpha']=0.0
plt.rcParams['axes.linewidth']=0.5
plt.rc('axes',edgecolor='grey')

plt.rcParams['axes.spines.top']= 0
plt.rcParams['axes.spines.right']= 0
plt.rcParams['axes.spines.left']= 1
plt.rcParams['axes.spines.bottom']= 1   

import itertools
from itertools import cycle
lines = ["s","o","d","<"]
linecycler = cycle(lines)

def veg_points(isvegc, dx = 1.0, veg_size = 10, ax = '', c = 'g'):
    """
    input:
        isvegc: binary array of vegetation field
    output:
        scatter plot of vegetated points
    """
    ncol = isvegc.shape[0]
    nrow = isvegc.shape[1]
    xc = np.arange(0, ncol*dx, dx)  + dx/2
    yc = np.arange(0, nrow*dx, dx)  + dx/2
    xc, yc = np.meshgrid(xc, yc)
    
    isvegc = (isvegc*veg_size).astype(float)

    isvegc[isvegc == 0] = np.nan
    if ax == '':      
        fig = plt.figure(figsize = (4,6))
        ax = fig.add_axes()
    else:
        fig = plt.gcf()        
           
    vegplot = ax.scatter(xc+dx/2, yc+dx/2.,
                        s = isvegc.T,
                        c = c,  marker='o', alpha = 0.75)
    
    ax.set_xlim(xc.min(), xc.max())
    ax.set_ylim(yc.min(), yc.max())
    ax.set_xticks([], []);
    ax.set_yticks([], []);

    return vegplot
    
        
    
def color_topo(zc, dx = 1.0, ax = ''):
    
    ncol = zc.shape[0]
    nrow = zc.shape[1]
    xc = np.arange(0, ncol*dx, dx)  + dx/2
    yc = np.arange(0, nrow*dx, dx)  + dx/2
    xc, yc = np.meshgrid(xc, yc)
    
    if ax == '':      
        fig = plt.figure(figsize = (4,6))
        ax = fig.add_subplot(111)
    else:
        fig = plt.gcf()
            
    topo = plt.pcolormesh(xc, yc, zc.T , cmap = 'Greys_r',
                          alpha = 0.5, 
                          shading ='gouraud')
    
    plt.xlim(xc.min(), xc.max())
    plt.ylim(yc.min(), yc.max())
    plt.xticks([], []);
    plt.yticks([], []);
    
    
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='5%', pad=0.05)

    cbh = fig.colorbar(topo, cax = cax, shrink=1)
    # cbh.set_label(clabel, fontsize = 16)
    # cbh = plt.colorbar(topo,shrink=1)
    cbh.set_label('elevation (m)', fontsize = 16)
    
    return topo

def plot_c(r, dx = 1, ax  = ''):

    import cv2
    thresh = r.astype(np.uint8)
    contours, hierarchy = cv2.findContours(thresh,cv2.RETR_TREE,
        cv2.CHAIN_APPROX_NONE)
    if ax == '':     
        fig, ax = plt.subplots()

    for n, contour in enumerate(contours):
        contour = np.squeeze((contour))

        contours[n] = contour
        if len(contour) >2:
            ax.plot(contour[:, 1]*dx + dx, contour[:, 0]*dx + dx, 
                linewidth=0.7, c= 'k')
            
def plot_c_inv(r, dx = 1, ax  = ''):

    import cv2
    thresh = 1- r.astype(np.uint8)
    image, contours, hierarchy = cv2.findContours(thresh,cv2.RETR_TREE,
        cv2.CHAIN_APPROX_NONE)
    if ax == '':     
        fig, ax = plt.subplots()

    for n, contour in enumerate(contours):
        contour = np.squeeze((contour))

        contours[n] = contour
        if len(contour) >2:
            ax.plot(contour[:, 1]*dx + dx, contour[:, 0]*dx + dx/1.5, 
                linewidth=0.7, c= 'k')

def colormap(df, array,  ax = '',
             colorbar = True, veg_scale = False,
             bounds = '', clabel = '',
             cmin = False, cmax = False,           
             cround = '',cfontsize = 14,
             cmap = cmocean.cm.deep):


    isvegc = df['isvegc'].astype(float)
    dx = df['dx']
    
    ncol = isvegc.shape[0]
    nrow = isvegc.shape[1]
    
    xc = np.arange(0, ncol*dx, dx)  + dx/2
    yc = np.arange(0, nrow*dx, dx)  + dx/2
    xc, yc = np.meshgrid(xc, yc)
    
    xc = xc.T
    yc = yc.T
    
    
    if isvegc.sum() == 0:
      veg_scale = False  
      
    if ax == '':      
        fig = plt.figure(figsize = (5,5))
        ax = fig.add_subplot(111)
    else:
        fig = plt.gcf()
  
    if bounds == '':
    
        if veg_scale == True:
            scale_vals = array[isvegc == True].ravel()
        else: 
            scale_vals =  array.ravel()
        
        if type(cmin) == bool:            
            cmin = np.nanmin(scale_vals)

        if type(cmax) == bool:            
            cmax = np.nanmax(scale_vals)

    
        bounds = np.linspace(cmin, cmax, 100)

        if cround != '':
            cmin = np.round(cmin, cround)
            bounds = np.arange(cmin, cmax, 1/10.**cround)
            
        if np.sum(array.astype(int) - array) == 0:
             bounds = np.arange(cmin, cmax+1.1,1)

            
        if np.nanstd(array) < 1e-5:
            bounds = np.linspace(cmin-1, cmax+1, 10)

    norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)    

    zinflplot = ax.pcolormesh(xc.T,
                           yc.T,
                           array.T,
                           norm = norm,
                           cmap=cmap, alpha= 1);
                                               
    if colorbar == True:
      from mpl_toolkits.axes_grid1 import make_axes_locatable
    
      divider = make_axes_locatable(ax)
      cax = divider.append_axes('right', size='5%', pad=0.05)

      cbh = fig.colorbar(zinflplot,cax = cax,shrink=1)
      cbh.set_label(clabel, fontsize = cfontsize)
      cbh.ax.tick_params(labelsize=cfontsize) 
      


    ax.set_ylim(yc.min(), yc.max())
    ax.set_xlim(xc.min(), xc.max())
    ax.set_xticks([], []);
    ax.set_yticks([], []);

    return zinflplot