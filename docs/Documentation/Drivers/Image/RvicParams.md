# RvicParams

The routing process is done in two steps: Generation unit hydrograph (ref) and the convolution (ref)
*   The unit hydrograph is generated separately with the RVIC streamflow routing model. This RVIC is a stand alone python program see:  http://rvic.readthedocs.io/en/latest/. 
*   Convolution is calculated in VIC itself (if routing is enabled in the makefile) and uses the parameter file generated in step 1

VIC needs the parameter file generated in step 1. In formation how to generate this file can be found at http://rvic.readthedocs.io/en/latest/user-guide/parameters/.
