# Routing

There are two routing options available: rout_stub (default) and rout_rvic.

## Option: rout_stub

Option rout_stub is the default option. During compilation the routing module is left empty: VIC can thus not be used for routing and consequently no discharge will be calculated.

## Option: rout_rvic

Option rout_rvic is optional. 

The routing process is done in two steps: Generation unit hydrograph and the convolution
*   The unit hydrograph is generated separately with the RVIC streamflow routing model. This RVIC is a stand alone python program see:  http://rvic.readthedocs.io/en/latest/. VIC is tested with RVIC version 1.1.0.
*   Convolution is calculated in VIC itself (if routing is enabled in the makefile) and uses the RVIC-parameter file generated in the step mentioned above.

### Running: rout_rvic

This option can be enabled by adding the following commandline option to the `make`:
`make ROUT=rout_rvic`

VIC needs the RVIC-parameter file. In formation how to generate this file can be found at http://rvic.readthedocs.io/en/latest/user-guide/parameters/. Note that the domain (lat-lon information and mask) of the RVIC-parameter file must match with the one of the VIC domain.

The RIVC-parameter file must be set in the [Global Parameter File](GlobalParam.md) e.g.:
`ROUT_PARAM                     RVIC_params.nc                   # Routing parameter path/file`

An extra output has to be set in the [Global Parameter File](GlobalParam.md):
`OUTVAR                  OUT_DISCHARGE`
