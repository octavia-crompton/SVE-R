__Brief summary__

The core of the SVE-R model is in dry.f, and a number of python wrapper functions are provided to facilitate its use.  Please refer to documentation (model.pdf) for further details on the SVE-R model and python scripts.

__Quick start__

From a bash terminal: 

1. Create a subfolder `test`, containing a params.json file (see example subfolders for examples).  
2. Open a terminal and nagivate to the model directory.
3. Call the command: `python call_dry.py test`

If the above lines throw an error, it is likeley because the Python things differed. Try creating a virtual environmentand installing the same python library versions:

__Installing Python dependencies__

The required Python libraries can be installed from the requirements file `requirements.txt` with pip:
`pip install -r requirements.txt`

To create a virtual environment that is compatible with the jupyter notebook files:

`conda create -n o_env  python=2.7 ipykernel`
`source activate ipykernel_py2`    # On Windows, remove the word 'source'
`python -m ipykernel install --user

More details in the documentation:  https://ipython.readthedocs.io/en/stable/install/kernel_install.html


Notes:
- The Fortran code requires the accelerate framework.



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

