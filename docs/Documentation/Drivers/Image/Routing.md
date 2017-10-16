# Routing

By default, river routing is not included in the VIC model. But the routing option can be turned on such that river routing is performed in real time when the VIC model is running, outputing routed river discharge. The routing scheme is identical to the RVIC routing model (http://rvic.readthedocs.io/en/latest/) and is tested with RVIC version 1.1.0. The steps to turn on the routing option are as follows.

## 1. Prepare routing-specific parameter file
To run routing in VIC, an extra parameter file that defines unit hydrograph information is needed. This parameter file contains information about unit hydrographs, which quantifies the reponse of a river outlet (whose dicharge is to be calculated) to an upstream runoff input. This parameter file can typically be generated using the RVIC model. Information on how to run RVIC and generate this file can be found at http://rvic.readthedocs.io/en/latest/. Note that the domain file (lat-lon information and mask) required by the RVIC model must match the one used for the VIC model. With this extra parameter file, VIC is able to convolve local runoff to calculate river discharge at downstream outlets.

## 2. Specify the routing-specific information in the VIC global file
The path of the routing-specific parameter file must be set in the [Global Parameter File](GlobalParam.md). For example:
```
ROUT_PARAM                     ./RVIC_params.nc                   # Routing parameter path/file
```

To output routed discharge results, an extra output variable has to be set in the [Global Parameter File](GlobalParam.md):
```

OUTVAR                  OUT_DISCHARGE
```

## 3. Turn on the routing option when compiling VIC
You must turn on the routing option when you compile VIC. There are two ways to turn on the routing option when compiling VIC:

1) Add the following command line option to the `make` command when compiling VIC:
`make ROUT=rout_rvic`

2) Change the value of `ROUT` in the Makefile to `rout_rvic`, i.e., `ROUT=rout_rvic`. This variable is by default `rout_stub`, which means running VIC without routing.

