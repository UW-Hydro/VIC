#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
#include <string.h>

static char vcid[] = "$Id$";

out_data_struct *create_output_list() {
/*************************************************************
  create_output_list()      Ted Bohn     September 08, 2006

  This routine creates the list of output variables.

*************************************************************/

  extern option_struct options;
  int v;
  out_data_struct *out_data;

  out_data = (out_data_struct *)calloc(N_OUTVAR_TYPES,sizeof(out_data_struct));

  // Build the list of supported output variables

  // Water Balance Terms - state variables
  strcpy(out_data[OUT_SNOW_CANOPY].varname,"OUT_SNOW_CANOPY");         /* snow interception storage in canopy [mm] */
  strcpy(out_data[OUT_SNOW_COVER].varname,"OUT_SNOW_COVER");           /* fractional area of snow cover [fraction] */
  strcpy(out_data[OUT_SNOW_DEPTH].varname,"OUT_SNOW_DEPTH");           /* depth of snow pack [cm] */
  strcpy(out_data[OUT_SOIL_ICE].varname,"OUT_SOIL_ICE");               /* soil ice content [mm] for each soil layer */
  strcpy(out_data[OUT_SOIL_LIQ].varname,"OUT_SOIL_LIQ");               /* soil liquid moisture content [mm] for each soil layer */
  strcpy(out_data[OUT_SOIL_MOIST].varname,"OUT_SOIL_MOIST");           /* soil total moisture content [mm] for each soil layer */
  strcpy(out_data[OUT_SWE].varname,"OUT_SWE");                         /* snow water equivalent in snow pack [mm] */
  strcpy(out_data[OUT_WDEW].varname,"OUT_WDEW");                       /* total moisture interception storage in canopy [mm] */

  // Water Balance Terms - fluxes
  strcpy(out_data[OUT_BASEFLOW].varname,"OUT_BASEFLOW");               /* baseflow out of the bottom layer [mm] */
  strcpy(out_data[OUT_EVAP].varname,"OUT_EVAP");                       /* total net evaporation [mm] */
  strcpy(out_data[OUT_EVAP_BARE].varname,"OUT_EVAP_BARE");             /* net evaporation from bare soil [mm] */
  strcpy(out_data[OUT_EVAP_CANOP].varname,"OUT_EVAP_CANOP");           /* net evaporation from canopy interception [mm] */
  strcpy(out_data[OUT_TRANSP_VEG].varname,"OUT_TRANSP_VEG");           /* net transpiration from vegetation [mm] */
  strcpy(out_data[OUT_INFLOW].varname,"OUT_INFLOW");                   /* moisture that reaches top of soil column [mm] */
  strcpy(out_data[OUT_PREC].varname,"OUT_PREC");                       /* incoming precipitation [mm] */
  strcpy(out_data[OUT_RUNOFF].varname,"OUT_RUNOFF");                   /* surface runoff [mm] */
  strcpy(out_data[OUT_SUB_CANOP].varname,"OUT_SUB_CANOP");             /* net sublimation from snow stored in canopy [mm] */
  strcpy(out_data[OUT_SUB_SNOW].varname,"OUT_SUB_SNOW");               /* net sublimation from snow pack [mm] */

  // Energy Balance Terms - state variables
  strcpy(out_data[OUT_ALBEDO].varname,"OUT_ALBEDO");                   /* albedo [fraction] */
  strcpy(out_data[OUT_FDEPTH].varname,"OUT_FDEPTH");                   /* depth of freezing fronts [m] for each freezing front */
  strcpy(out_data[OUT_RAD_TEMP].varname,"OUT_RAD_TEMP");               /* average radiative surface temperature [K] */
  strcpy(out_data[OUT_SOIL_TEMP].varname,"OUT_SOIL_TEMP");             /* soil temperature [C] for each soil layer */
  strcpy(out_data[OUT_SURF_TEMP].varname,"OUT_SURF_TEMP");             /* average surface temperature [C] */
  strcpy(out_data[OUT_TDEPTH].varname,"OUT_TDEPTH");                   /* depth of thawing fronts [m] for each thawing front */

  // Energy Balance Terms - fluxes
  strcpy(out_data[OUT_ADVECTION].varname,"OUT_ADVECTION");             /* advected energy [W/m2] */
  strcpy(out_data[OUT_DELTACC].varname,"OUT_DELTACC");                 /* rate of change in cold content in snow pack [W/m2] */
  strcpy(out_data[OUT_DELTAH].varname,"OUT_DELTAH");                   /* rate of change in heat storage [W/m2] */
  strcpy(out_data[OUT_ENERGY_ERROR].varname,"OUT_ENERGY_ERROR");       /* energy budget error [W/m2] */
  strcpy(out_data[OUT_GRND_FLUX].varname,"OUT_GRND_FLUX");             /* net heat flux into ground [W/m2] */
  strcpy(out_data[OUT_LATENT].varname,"OUT_LATENT");                   /* net upward latent heat flux [W/m2] */
  strcpy(out_data[OUT_LONGWAVE].varname,"OUT_LONGWAVE");               /* incoming longwave radiation [W/m2] */
  strcpy(out_data[OUT_NET_LONG].varname,"OUT_NET_LONG");               /* net downward longwave flux [W/m2] */
  strcpy(out_data[OUT_NET_SHORT].varname,"OUT_NET_SHORT");             /* net downward shortwave flux [W/m2] */
  strcpy(out_data[OUT_R_NET].varname,"OUT_R_NET");                     /* net downward radiation flux [W/m2] */
  strcpy(out_data[OUT_REFREEZE_ENERGY].varname,"OUT_REFREEZE_ENERGY"); /* net energy used to refreeze liquid water in snowpack [W/m2] */
  strcpy(out_data[OUT_SENSIBLE].varname,"OUT_SENSIBLE");               /* net upward sensible heat flux [W/m2] */
  strcpy(out_data[OUT_SHORTWAVE].varname,"OUT_SHORTWAVE");             /* incoming shortwave radiation [W/m2] */
  strcpy(out_data[OUT_SNOW_FLUX].varname,"OUT_SNOW_FLUX");             /* energy flux through snow pack [W/m2] */

  // Miscellaneous Terms
  strcpy(out_data[OUT_AERO_RESIST].varname,"OUT_AERO_RESIST");         /* canopy aerodynamic resistance [s/m] */
  strcpy(out_data[OUT_AERO_COND].varname,"OUT_AERO_COND");             /* canopy aerodynamic conductance [m/s] */
  strcpy(out_data[OUT_AIR_TEMP].varname,"OUT_AIR_TEMP");               /* air temperature [C] */
  strcpy(out_data[OUT_DENSITY].varname,"OUT_DENSITY");                 /* near-surface atmospheric density [kg/m3] */
  strcpy(out_data[OUT_PRESSURE].varname,"OUT_PRESSURE");               /* near surface atmospheric pressure [kPa] */
  strcpy(out_data[OUT_REL_HUMID].varname,"OUT_REL_HUMID");             /* relative humidity [fraction]*/
  strcpy(out_data[OUT_SURF_COND].varname,"OUT_SURF_COND");             /* surface conductance [m/s] */
  strcpy(out_data[OUT_VP].varname,"OUT_VP");                           /* near surface vapor pressure [kPa] */
  strcpy(out_data[OUT_WIND].varname,"OUT_WIND");                       /* near surface wind speed [m/s] */

  // Band-specific quantities
  strcpy(out_data[OUT_ADVECTION_BAND].varname,"OUT_ADVECTION_BAND");             /* advected energy [W/m2] */
  strcpy(out_data[OUT_ALBEDO_BAND].varname,"OUT_ALBEDO_BAND");                   /* albedo [fraction] */
  strcpy(out_data[OUT_DELTACC_BAND].varname,"OUT_DELTACC_BAND");                 /* change in cold content in snow pack [W/m2] */
  strcpy(out_data[OUT_GRND_FLUX_BAND].varname,"OUT_GRND_FLUX_BAND");             /* net heat flux into ground [W/m2] */
  strcpy(out_data[OUT_LATENT_BAND].varname,"OUT_LATENT_BAND");                   /* net upward latent heat flux [W/m2] */
  strcpy(out_data[OUT_NET_LONG_BAND].varname,"OUT_NET_LONG_BAND");               /* net downward longwave flux [W/m2] */
  strcpy(out_data[OUT_NET_SHORT_BAND].varname,"OUT_NET_SHORT_BAND");             /* net downward shortwave flux [W/m2] */
  strcpy(out_data[OUT_REFREEZE_ENERGY_BAND].varname,"OUT_REFREEZE_ENERGY_BAND"); /* net energy used to refreeze liquid water in snowpack [W/m2] */
  strcpy(out_data[OUT_SENSIBLE_BAND].varname,"OUT_SENSIBLE_BAND");               /* net upward sensible heat flux [W/m2] */
  strcpy(out_data[OUT_SNOW_CANOPY_BAND].varname,"OUT_SNOW_CANOPY_BAND");         /* snow interception storage in canopy [mm] */
  strcpy(out_data[OUT_SNOW_COVER_BAND].varname,"OUT_SNOW_COVER_BAND");           /* fractional area of snow cover [fraction] */
  strcpy(out_data[OUT_SNOW_DEPTH_BAND].varname,"OUT_SNOW_DEPTH_BAND");           /* depth of snow pack [cm] */
  strcpy(out_data[OUT_SNOW_FLUX_BAND].varname,"OUT_SNOW_FLUX_BAND");             /* energy flux through snow pack [W/m2] */
  strcpy(out_data[OUT_SWE_BAND].varname,"OUT_SWE_BAND");                         /* snow water equivalent in snow pack [mm] */

  // Set number of elements
  for (v=0; v<N_OUTVAR_TYPES; v++) {
    out_data[v].nelem = 1;
  }
  if (options.FROZEN_SOIL) {
    out_data[OUT_FDEPTH].nelem = MAX_FRONTS;
    out_data[OUT_TDEPTH].nelem = MAX_FRONTS;
  }
  out_data[OUT_SOIL_ICE].nelem = options.Nlayer;
  out_data[OUT_SOIL_LIQ].nelem = options.Nlayer;
  out_data[OUT_SOIL_MOIST].nelem = options.Nlayer;
  out_data[OUT_SOIL_TEMP].nelem = options.Nlayer;
  out_data[OUT_ADVECTION_BAND].nelem = options.SNOW_BAND;
  out_data[OUT_ALBEDO_BAND].nelem = options.SNOW_BAND;
  out_data[OUT_DELTACC_BAND].nelem = options.SNOW_BAND;
  out_data[OUT_GRND_FLUX_BAND].nelem = options.SNOW_BAND;
  out_data[OUT_LATENT_BAND].nelem = options.SNOW_BAND;
  out_data[OUT_NET_LONG_BAND].nelem = options.SNOW_BAND;
  out_data[OUT_NET_SHORT_BAND].nelem = options.SNOW_BAND;
  out_data[OUT_REFREEZE_ENERGY_BAND].nelem = options.SNOW_BAND;
  out_data[OUT_SENSIBLE_BAND].nelem = options.SNOW_BAND;
  out_data[OUT_SNOW_CANOPY_BAND].nelem = options.SNOW_BAND;
  out_data[OUT_SNOW_COVER_BAND].nelem = options.SNOW_BAND;
  out_data[OUT_SNOW_DEPTH_BAND].nelem = options.SNOW_BAND;
  out_data[OUT_SNOW_FLUX_BAND].nelem = options.SNOW_BAND;
  out_data[OUT_SWE_BAND].nelem = options.SNOW_BAND;

  // Allocate space for data
  for (v=0; v<N_OUTVAR_TYPES; v++) {
    out_data[v].data = (double *)calloc(out_data[v].nelem, sizeof(double));
  }

  // Initialize data values
  init_output_list(out_data, FALSE, "%.1f", OUT_TYPE_FLOAT, 1);

  return out_data;

}


void init_output_list(out_data_struct *out_data, int write, char *format, int type, float mult) {
/*************************************************************
  init_output_list()      Ted Bohn     September 08, 2006

  This routine initializes the output information for all output variables.

*************************************************************/
  int varid, i;

  for (varid=0; varid<N_OUTVAR_TYPES; varid++) {
    out_data[varid].write = write;
    strcpy(out_data[varid].format,format);
    out_data[varid].type = type;
    out_data[varid].mult = mult;
    for(i=0; i<out_data[varid].nelem; i++) {
      out_data[varid].data[i] = 0;
    }
  }

}


int set_output_var(out_data_file_struct *out_data_files,
                    int write,
                    int filenum,
                    out_data_struct *out_data,
                    char *varname,
                    int varnum,
                    char *format,
                    int type,
                    float mult) {
/*************************************************************
  set_output_var()      Ted Bohn     September 08, 2006

  This routine updates the output information for a given output variable.

*************************************************************/
  int varid;
  int found=FALSE;
  int status=0;

  for (varid=0; varid<N_OUTVAR_TYPES; varid++) {
    if (strcmp(out_data[varid].varname,varname) == 0) {
      found = TRUE;
      out_data[varid].write = write;
      if (strcmp(format,"*") != 0)
        strcpy(out_data[varid].format,format);
      if (type != 0)
        out_data[varid].type = type;
      if (mult != 0)
        out_data[varid].mult = mult;
      out_data_files[filenum].varid[varnum] = varid;
    }
  }
  if (!found) {
    status = -1;
    fprintf(stderr, "Error: set_output_var: \"%s\" was not found in the list of supported output variable names.  Please use the exact name listed in vicNl_def.h.\n", varname);
  }
  return status;

}


void zero_output_list(out_data_struct *out_data) {
/*************************************************************
  zero_output_list()      Ted Bohn     September 08, 2006

  This routine resets the values of all output variables to 0.

*************************************************************/
  int varid, i;

  for (varid=0; varid<N_OUTVAR_TYPES; varid++) {
    for(i=0; i<out_data[varid].nelem; i++) {
      out_data[varid].data[i] = 0;
    }
  }

}

void free_out_data_files(out_data_file_struct **out_data_files) {
/*************************************************************
  free_out_data_files()      Ted Bohn     September 08, 2006

  This routine frees the memory in the out_data_files array.

*************************************************************/
  extern option_struct options;
  int filenum;

  for (filenum=0; filenum<options.Noutfiles; filenum++) {
    free((char*)(*out_data_files)[filenum].varid);
  }
  free((char*)(*out_data_files));

}

