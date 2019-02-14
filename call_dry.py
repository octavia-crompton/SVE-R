import pickle
import os, sys
from os.path import dirname
import argparse
from commands import getoutput as cmd
import multiprocessing as mp
import numpy as np
import json
import datetime
import shutil
import zipfile
import time

from glob import glob
current_dir = os.getcwd()

modules = [ 'sim_read', 'sim_input']
for mod in modules:
    if mod in sys.modules:
        del sys.modules[mod]

from sim_read import *
from sim_input import *

def main(argv):
    """
    This script writes the fortan input files, compiles and executes the fortran code dry.for,
    and read the output files into a pandas dataframe, which is saved as a pklz file.

    input:
         sim_name: name of the simulation folder, which must contain  a parameter file (params.json) 
    """
    parser = argparse.ArgumentParser()
  
    parser.add_argument("sim_name", type=str, help="name of simulation folder")
    args = parser.parse_args()        
    
    sim_path =  '/'.join([current_dir, args.sim_name])    
    
    sim_input(sim_path)    # write the fortran input files
    runmodel(sim_path)     # compile and run fortran code
    sim_read(sim_path)

def runmodel(sim_path):   

    """ 
    compile dry.for source code
    execute ./sw
    """         
    c = os.system("gfortran -o {0}/sw  -framework accelerate {0}/dry.for".format(sim_path))

    a = os.system("cd {0} \n ./sw  \n cd {1}".format(sim_path, current_dir))
    
    return   a



if __name__ == '__main__':
    main(sys.argv)
