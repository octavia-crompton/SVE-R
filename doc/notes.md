## coupling the model
### time grids
The structure of the SVE solver is presented here in pseudocode and as a schematic.  

The coupled SVE-R operates on several time grids.  The SVE solver requires the smallest timesteps for CFL stability (as small as 0.5 ms for a steep slope).  The Richards equation solver, by contrast, does not requires small time steps.  To limit the size of the SVE output files, the variables are saved every \code{dt\_p} seconds (or every \code{nprt}  SVE timesteps).

The time variables are specified/related as:

\begin{itemize}
	\item \code{dt\_sw} : SVE solver timestep.  
	\subitem  \code{it = t/dt\_sw} : iteration number of the SVE solver
	\item \code{dt\_r} : Richards equation timestep 
	\item \code{iscale = dt\_r/dt\_sw} : ratio of the SVE to the Richards equation timesteps.
	\subitem  \code{itp = t/dt\_p} : print iteration (number of times the output has already been printed to file.
	\subitem  \code{nprt } : number of SVE timesteps between saving the output
\end{itemize}



### running the model


\code{  flux1, flux2, flux3, flux4} are the total fluxes out of boundaries 1, 2, 3, 4, where bounbdary 1 refers to all the cells with boundary type = 1. 

\code{xflux0, yflux0, xflux1, xflux2 } are the fluxes between cells
