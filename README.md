Brief summary:

The core of the SVE-R model is in dry.f, and a number of python wrapper functions are provided to facilitate its use.  The SVE-R model and python scripts are documented in detail in doc/model.pdf.

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
* example\_2by2:
* example\_image:
* example\_infl:
* example\_rain:

Notes:
 
This code is provided as a supplement to the paper: "Name of Paper" (available 
 [here](http://example.com "Title")).
 The terminology used in the code vs. the companion paper differ in several ways:
First, impermeable and permeable areas are labeled here as bare soil and vegetated patches, respectively.  
The model simulations described in the emulator paper were run in parallel, whereas the code provided here is is for a single 

