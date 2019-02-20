__Brief summary__

The core of the SVE-R model is in dry.f, and a number of python wrapper functions are provided to facilitate its use.  Please refer to documentation (model.pdf) for further details on the SVE-R model and python scripts.

__Quick start__

From a bash terminal: 

1. Create a subfolder `test`, containing a params.json file (see example subfolders for examples).  
2. Open a terminal and nagivate to the model directory.
3. Call the command: `python call_dry.py test`

If the above lines throw an error, it is likely due a difference in Python or Fortran versions.


The Fortran code was developed using a gfortran compiler, and has not been tested on other compilers.  Additionally, it requires the [accelarate](https://sites.ualberta.ca/~kbeach/lapack.html "Title") framework:

__Python dependencies__

The Python wrapper scripts and Jupyter notebooks use Python 2.7 and the libraries listed in `requirements.txt`.  
The easiest way to match environments is to create a virtual environment and install the required Python libraries with pip:

  `pip install -r requirements.txt`

To create a virtual environment that is compatible with the Jupyter notebook examples:

`conda create -n o_env  python=2.7 ipykernel` 
`source activate ipykernel_py2`    # On Windows, remove the word 'source'
`python -m ipykernel install --user`

Please refer to the Conda documentation for [more details](https://ipython.readthedocs.io/en/stable/install/kernel_install.html "Title")


__List of files:__

dry.for  : SVE-R Fortran core model
dry.inc  : specifies common variables for dry.for


Python files:
wrapper scripts:
* call\_dry.py
   * sim\_input.py
      * input\_phi.py
      * input\_veg.py
      * input\_coords.py
      * input\_boundary.py
   * sim\_read.py


Other files:
- plot\_functions.py  : plotting functions used by ipython example notebooks
- feature\_functions.py : python functions to convert a binary array (e.g. a vegetation map or urban landscape) into an array of features on which the random forest can operate.
- feature\_example.ipynb : contains visualizations of the features.

Examples: 
* example\_2by2:  illustrates the grid contruction for a 2x2 rectangular grid (described in more detail in the 
 [documentation](https://github.com/octavia-crompton/SVE-R/blob/master/doc/model.pdf "Title")
* example\_image: rain driven overland flow on a planar hillslope with the vegetation field read from an image file.
* example\_inflow: overland flow on a planar hillslope supplied by a fixed-flux boundary condition at the top of the hillslope.   
* example\_rain:  rain driven overland flow on a planar hillslope with the randomly-generated vegetation field.

Notes:
 
This code is provided as a supplement to the paper: "Name of Paper" (available 
 [here](http://example.com "Title"), referred to here as the emulator paper).
 The terminology used in the code vs. the companion paper differ in several ways:
Firstly, the model is referred to in the documentation as the SVE-R model (Saint Venant Equation - Richards equation), and in the paper as the "SVE" model for simplicity. 
Second, the emulator paper describe impermeable and permeable areas, which are labeled here as bare soil and vegetated patches.  
The model simulations described in the emulator paper were run in parallel, whereas the code provided here is is for a single simulation.

