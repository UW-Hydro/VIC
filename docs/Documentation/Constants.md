# Model Constants and Miscellaneous Parameters

The VIC model requires the specification of many model parameters to run. While some parameters are spatially distributed (e.g. _soil_, _vegetation_, or _snowbands_), others are more appropriately held constant across the model domain.

## Physical Constants

Physical constants are set in the `vic_physical_constants.h` header file. These constants rarely need to be changed and the model must be recompiled after any modifications.

## Model Constants and Parameters

Some model constants and model parameters need to be changed more frequently. VIC supports this through the use of the optional `CONSTANTS` file. The path to this file should be specified in the global parameter file.

There are over 150 constants that can be changed through this file. An example `CONSTANTS` file is shown below.

```
# Set HUGE_RESIST to smaller value
HUGE_RESIST 1e5

# Set EMISS_GRND to unrealistic value
EMISS_GRND 0.3
```

The table below lists the constants available for manipulation via the `CONSTANTS` file. Descriptions will be added later or as needed. Default values can be found in `initialize_parameters.c`

| Constant                     | Description |
|------------------------------|-------------|
| LAPSE_RATE                   |             |
| GAUGE_HEIGHT                 |             |
| HUGE_RESIST                  |             |
| ALBEDO_BARE_SOIL             |             |
| ALBEDO_H20_SURF              |             |
| EMISS_GRND                   |             |
| EMISS_ICE                    |             |
| EMISS_VEG                    |             |
| EMISS_SNOW                   |             |
| EMISS_H2O                    |             |
| SOIL_RESID_MOIST             |             |
| SOIL_SLAB_MOIST_FRACT        |             |
| VEG_LAI_SNOW_MULTIPLIER      |             |
| VEG_MIN_INTERCEPTION_STORAGE |             |
| VEG_LAI_WATER_FACTOR         |             |
| CANOPY_CLOSURE               |             |
| CANOPY_RSMAX                 |             |
| CANOPY_VPDMINFACTOR          |             |
| LAKE_TMELT                   |             |
| LAKE_MAX_SURFACE             |             |
| LAKE_BETA                    |             |
| LAKE_FRACMIN                 |             |
| LAKE_FRACLIM                 |             |
| LAKE_DM                      |             |
| LAKE_SNOWCRIT                |             |
| LAKE_ZWATER                  |             |
| LAKE_ZSNOW                   |             |
| LAKE_RHOSNOW                 |             |
| LAKE_CONDI                   |             |
| LAKE_CONDS                   |             |
| LAKE_LAMISW                  |             |
| LAKE_LAMILW                  |             |
| LAKE_LAMSSW                  |             |
| LAKE_LAMSLW                  |             |
| LAKE_LAMWSW                  |             |
| LAKE_LAMWLW                  |             |
| LAKE_A1                      |             |
| LAKE_A2                      |             |
| LAKE_QWTAU                   |             |
| LAKE_MAX_ITER                |             |
| SVP_A                        |             |
| SVP_B                        |             |
| SVP_C                        |             |
| PHOTO_OMEGA                  |             |
| PHOTO_LAIMAX                 |             |
| PHOTO_LAILIMIT               |             |
| PHOTO_LAIMIN                 |             |
| PHOTO_EPAR                   |             |
| PHOTO_FCMAX                  |             |
| PHOTO_FCMIN                  |             |
| PHOTO_ZENITHMIN              |             |
| PHOTO_ZENITHMINPAR           |             |
| PHOTO_ALBSOIPARMIN           |             |
| PHOTO_MINMAXETRANS           |             |
| PHOTO_MINSTOMCOND            |             |
| PHOTO_FCI1C3                 |             |
| PHOTO_FCI1C4                 |             |
| PHOTO_OX                     |             |
| PHOTO_KC                     |             |
| PHOTO_KO                     |             |
| PHOTO_EC                     |             |
| PHOTO_EO                     |             |
| PHOTO_EV                     |             |
| PHOTO_ER                     |             |
| PHOTO_ALC3                   |             |
| PHOTO_FRDC3                  |             |
| PHOTO_EK                     |             |
| PHOTO_ALC4                   |             |
| PHOTO_FRDC4                  |             |
| PHOTO_THETA                  |             |
| PHOTO_FRLEAF                 |             |
| PHOTO_FRGROWTH               |             |
| SRESP_E0_LT                  |             |
| SRESP_T0_LT                  |             |
| SRESP_WMINFM                 |             |
| SRESP_WMAXFM                 |             |
| SRESP_WOPTFM                 |             |
| SRESP_RHSAT                  |             |
| SRESP_RFACTOR                |             |
| SRESP_TAULITTER              |             |
| SRESP_TAUINTER               |             |
| SRESP_TAUSLOW                |             |
| SRESP_FAIR                   |             |
| SRESP_FINTER                 |             |
| SNOW_MAX_SURFACE_SWE         |             |
| SNOW_LIQUID_WATER_CAPACITY   |             |
| SNOW_NEW_SNOW_DENSITY        |             |
| SNOW_NEW_SNOW_DENS_MAX       |             |
| SNOW_DEPTH_THRES             |             |
| SNOW_DENS_DMLIMIT            |             |
| SNOW_DENS_DMLIMIT_FACTOR     |             |
| SNOW_DENS_MAX_CHANGE         |             |
| SNOW_DENS_ETA0               |             |
| SNOW_DENS_C1                 |             |
| SNOW_DENS_C2                 |             |
| SNOW_DENS_C3                 |             |
| SNOW_DENS_C3_CONST           |             |
| SNOW_DENS_C4                 |             |
| SNOW_DENS_C4WET              |             |
| SNOW_DENS_C5                 |             |
| SNOW_DENS_C6                 |             |
| SNOW_DENS_F                  |             |
| SNOW_DENS_EXP                |             |
| SNOW_DENS_DENOM              |             |
| SNOW_NEW_SNT_C1              |             |
| SNOW_NEW_SNT_C2              |             |
| SNOW_NEW_SNT_C3              |             |
| SNOW_NEW_BRAS_DENOM          |             |
| SNOW_MIN_SWQ_EB_THRES        |             |
| SNOW_A1                      |             |
| SNOW_A2                      |             |
| SNOW_L1                      |             |
| SNOW_L2                      |             |
| SNOW_NEW_SNOW_ALB            |             |
| SNOW_ALB_ACCUM_A             |             |
| SNOW_ALB_ACCUM_B             |             |
| SNOW_ALB_THAW_A              |             |
| SNOW_ALB_THAW_B              |             |
| SNOW_TRACESNOW               |             |
| SNOW_CONDUCT                 |             |
| SNOW_MAX_SNOW_TEMP           |             |
| SNOW_MIN_RAIN_TEMP           |             |
| BLOWING_KA                   |             |
| BLOWING_CSALT                |             |
| BLOWING_UTHRESH              |             |
| BLOWING_KIN_VIS              |             |
| BLOWING_MAX_ITER             |             |
| BLOWING_K                    |             |
| BLOWING_SETTLING             |             |
| BLOWING_NUMINCS              |             |
| TREELINE_TEMPERATURE         |             |
| SNOW_DT                      |             |
| SURF_DT                      |             |
| SOIL_DT                      |             |
| CANOPY_DT                    |             |
| CANOPY_VP                    |             |
| TOL_GRND                     |             |
| TOL_OVER                     |             |
| FROZEN_MAXITER               |             |
| MAX_ITER_GRND_CANOPY         |             |
| NEWT_RAPH_MAXTRIAL           |             |
| NEWT_RAPH_TOLX               |             |
| NEWT_RAPH_TOLF               |             |
| NEWT_RAPH_R_MAX              |             |
| NEWT_RAPH_R_MIN              |             |
| NEWT_RAPH_RELAX1             |             |
| NEWT_RAPH_RELAX2             |             |
| NEWT_RAPH_RELAX3             |             |
| NEWT_RAPH_EPS2               |             |
| ROOT_BRENT_MAXTRIES          |             |
| ROOT_BRENT_MAXITER           |             |
| ROOT_BRENT_TSTEP             |             |
| ROOT_BRENT_T                 |             |
