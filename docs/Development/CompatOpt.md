# Backwards Compatibility of 4.1.2 with Earlier Versions of VIC

VIC versions 4.1.1 and later should be able to run using the same input files used for versions 4.0.6 and later. However, VIC versions 4.1.1 and later will not necessarily reproduce the results of earlier versions, due to their new features and bug fixes. To reproduce the results of earlier versions, certain options (found either in the [global parameter file](../Documentation/GlobalParam.md)) must be set to the values listed below. This still will not exactly reproduce the results of earlier versions, but will be fairly close.

## Options for compatibility with 4.1.1

To reproduce results from VIC 4.1.1, use the following option settings:

### Global Parameter File
| Name                      | Value                     |
|-----------------------    |-------------------------- |
| MTCLIM_SWE_CORR           | FALSE                     |
| VP_ITER                   | VP_ITER_ALWAYS (or omit)  |
| VP_INTERP                 | FALSE                     |
| LW_TYPE                   | LW_TVA (or omit)          |
| LW_CLOUD                  | LW_CLOUD_BRAS             |
| SW_PREC_THRESH            | 0 (or omit)               |

To reproduce results from VIC 4.1.0_r5, use the following option settings:

### Global Parameter File
| Name                  	| Value                    	|
|-----------------------	|--------------------------	|
| MTCLIM_SWE_CORR       	| FALSE                    	|
| VP_ITER               	| VP_ITER_ALWAYS (or omit) 	|
| VP_INTERP             	| FALSE (or omit)          	|
| LW_TYPE               	| LW_TVA (or omit)         	|
| LW_CLOUD              	| LW_CLOUD_BRAS            	|
| SW_PREC_THRESH        	| 0 (or omit)              	|
| AERO_RESIST_CANSNOW   	| AR_410                   	|
| GRND_FLUX_TYPE        	| GF_410 (or omit)         	|
| IMPLICIT              	| FALSE                    	|
| COMPUTE_TREELINE      	| FALSE (or omit)          	|
| JULY_TAVG_SUPPLIED    	| FALSE (or omit)          	|
| PLAPSE                	| FALSE                    	|

## Options for compatibility with 4.0.6

To reproduce results from VIC 4.0.6, use the following option settings:

### Global Parameter File
| Name                  	| Value                    	|
|-----------------------	|--------------------------	|
| MTCLIM_SWE_CORR       	| FALSE                    	|
| VP_ITER               	| VP_ITER_ALWAYS (or omit) 	|
| VP_INTERP             	| FALSE (or omit)          	|
| LW_TYPE               	| LW_TVA (or omit)         	|
| LW_CLOUD              	| LW_CLOUD_BRAS            	|
| SW_PREC_THRESH        	| 0 (or omit)              	|
| AERO_RESIST_CANSNOW   	| AR_406                   	|
| GRND_FLUX_TYPE        	| GF_406                   	|
| IMPLICIT              	| FALSE                    	|
| EXP_TRANS             	| FALSE (or omit)          	|
| EQUAL_AREA            	| FALSE (or omit)          	|
| LAKES                 	| FALSE (or omit)          	|
| LAKE_PROFILE          	| FALSE (or omit)          	|
| BLOWING               	| FALSE (or omit)          	|
| PLAPSE                	| FALSE                    	|
| SNOW_ALBEDO           	| USACE (or omit)          	|
| SNOW_DENSITY          	| DENS_BRAS (or omit)      	|

### user_def.h File
| Name                  	| Value                    	|
|-----------------------	|--------------------------	|
| CLOSE_ENERGY          	| FALSE                    	|
| SPATIAL_FROST         	| FALSE                    	|
| SPATIAL_SNOW          	| FALSE                    	|
| EXCESS_ICE            	| FALSE                    	|
