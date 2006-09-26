/* RCS Id String
 * $Id$
 */
/********************************************************************
  Modifications:
  2006-Sep-11 Implemented flexible output configuration. Added output
              variable types.  Added binary output format types. 
              Removed all output files except the state file from
              the outfiles_struct and the filenames_struct.  Added
              Noutfiles to the option_struct.  Created new
              out_data_struct and out_data_files_struct.  TJB

********************************************************************/
/***** Version Information *****/
#define VERSION		"VIC Release 4.0.6"

/***** Model Constants *****/
#define MAXSTRING    2048
#define MINSTRING    20
#define HUGE_RESIST  1.e20	/* largest allowable double number */
#define SMALL        1.e-12	/* smallest allowable double number */
#define MISSING      -99999.	/* missing value for multipliers in
				   BINARY format */
#define LITTLE 1		/* little-endian flag */
#define BIG 2			/* big-endian flag */
#define BROOKS 1		/* Brooks-Corey parameters for unsaturated flow
				 */ 
#define ASCII 1			/* met file format flag */
#define BINARY 2		/* met file format flag */

/***** Forcing Variable Types *****/
#define N_FORCING_TYPES 13
#define AIR_TEMP  0  /* air temperature per time step [C] */
#define ALBEDO    1  /* surface albedo [fraction] */
#define DENSITY   2  /* atmospheric density [kg/m3] */
#define LONGWAVE  3  /* incoming longwave radiation [W/m2] */
#define PREC      4  /* precipitation [mm] */
#define PRESSURE  5  /* atmospheric pressure [kPa] */
#define SHORTWAVE 6  /* incoming shortwave [W/m2] */
#define TMAX      7  /* maximum daily temperature [C] */
#define TMIN      8  /* minimum daily temperature [C] */
#define TSKC      9  /* cloud cover [fraction] */
#define VP        10 /* vapor pressure [kPa] */
#define WIND      11 /* wind speed [m/s] */
#define SKIP      12 /* place holder for unused data columns */

/***** Output Variable Types *****/
#define N_OUTVAR_TYPES 100
// Water Balance Terms - state variables
#define OUT_SNOW_CANOPY      0  /* snow interception storage in canopy [mm] */
#define OUT_SNOW_COVER       1  /* fractional area of snow cover [fraction] */
#define OUT_SNOW_DEPTH       2  /* depth of snow pack [cm] */
#define OUT_SOIL_ICE         3  /* soil ice content [mm] for each soil layer */
#define OUT_SOIL_LIQ         4  /* soil liquid content [mm] for each soil layer */
#define OUT_SOIL_MOIST       5  /* soil total moisture content [mm] for each soil layer */
#define OUT_SWE              6  /* snow water equivalent in snow pack [mm] */
#define OUT_WDEW             7  /* total moisture interception storage in canopy [mm] */
// Water Balance Terms - fluxes
#define OUT_BASEFLOW         8  /* baseflow out of the bottom layer [mm] */
#define OUT_EVAP             9  /* total net evaporation [mm] */
#define OUT_EVAP_BARE       10  /* net evaporation from bare soil [mm] */
#define OUT_EVAP_CANOP      11  /* net evaporation from canopy interception [mm] */
#define OUT_TRANSP_VEG      12  /* net transpiration from vegetation [mm] */
#define OUT_INFLOW          13  /* moisture that reaches top of soil column [mm] */
#define OUT_PREC            14  /* incoming precipitation [mm] */
#define OUT_RUNOFF          15  /* surface runoff [mm] */
#define OUT_SUB_CANOP       16  /* net sublimation from snow stored in canopy [mm] */
#define OUT_SUB_SNOW        17  /* net sublimation from snow pack [mm] */
// Energy Balance Terms - state variables
#define OUT_ALBEDO          18  /* albedo [fraction] */
#define OUT_FDEPTH          19  /* depth of freezing fronts [m] for each freezing front */
#define OUT_RAD_TEMP        20  /* average radiative surface temperature [K] */
#define OUT_SOIL_TEMP       21  /* soil temperature [C] for each soil layer */
#define OUT_SURF_TEMP       22  /* average surface temperature [C] */
#define OUT_TDEPTH          23  /* depth of thawing fronts [m] for each thawing front */
// Energy Balance Terms - fluxes
#define OUT_ADVECTION       24  /* advected energy [W/m2] */
#define OUT_DELTACC         25  /* rate of change in cold content in snow pack [W/m2] */
#define OUT_DELTAH          26  /* rate of change in heat storage [W/m2] */
#define OUT_ENERGY_ERROR    27  /* energy budget error [W/m2] */
#define OUT_GRND_FLUX       28  /* net heat flux into ground [W/m2] */
#define OUT_LATENT          29  /* net upward latent heat flux [W/m2] */
#define OUT_LONGWAVE        30  /* incoming longwave radiation [W/m2] */
#define OUT_NET_LONG        31  /* net downward longwave flux [W/m2] */
#define OUT_NET_SHORT       32  /* net downward shortwave flux [W/m2] */
#define OUT_R_NET           33  /* net downward radiation flux [W/m2] */
#define OUT_REFREEZE_ENERGY 34  /* net energy used to refreeze liquid water in snowpack [W/m2] */
#define OUT_SENSIBLE        35  /* net upward sensible heat flux [W/m2] */
#define OUT_SHORTWAVE       36  /* incoming shortwave radiation [W/m2] */
#define OUT_SNOW_FLUX       37  /* energy flux through snow pack [W/m2] */
// Miscellaneous Terms
#define OUT_AERO_RESIST     38  /* canopy aerodynamic resistance [s/m] */
#define OUT_AERO_COND       39  /* canopy aerodynamic conductance [m/s] */
#define OUT_AIR_TEMP        40  /* air temperature [C] */
#define OUT_DENSITY         41  /* near-surface atmospheric density [kg/m3]*/
#define OUT_PRESSURE        42  /* near surface atmospheric pressure [kPa] */
#define OUT_REL_HUMID       43  /* relative humidity [fraction]*/
#define OUT_SURF_COND       44  /* surface conductance [m/s] */
#define OUT_VP              45  /* near surface vapor pressure [kPa] */
#define OUT_WIND            46  /* near surface wind speed [m/s] */
// Band-specific quantities
#define OUT_ADVECTION_BAND       47  /* advected energy [W/m2] */
#define OUT_ALBEDO_BAND          48  /* albedo [fraction] */
#define OUT_DELTACC_BAND         49  /* change in cold content in snow pack [W/m2] */
#define OUT_GRND_FLUX_BAND       50  /* net heat flux into ground [W/m2] */
#define OUT_LATENT_BAND          51  /* net upward latent heat flux [W/m2] */
#define OUT_NET_LONG_BAND        52  /* net downward longwave flux [W/m2] */
#define OUT_NET_SHORT_BAND       53  /* net downward shortwave flux [W/m2] */
#define OUT_REFREEZE_ENERGY_BAND 54  /* net energy used to refreeze liquid water in snowpack [W/m2] */
#define OUT_SENSIBLE_BAND        55  /* net upward sensible heat flux [W/m2] */
#define OUT_SNOW_CANOPY_BAND     56  /* snow interception storage in canopy [mm] */
#define OUT_SNOW_COVER_BAND      57  /* fractional area of snow cover [fraction] */
#define OUT_SNOW_DEPTH_BAND      58  /* depth of snow pack [cm] */
#define OUT_SNOW_FLUX_BAND       59  /* energy flux through snow pack [W/m2] */
#define OUT_SWE_BAND             60  /* snow water equivalent in snow pack [mm] */

/***** Physical Constants *****/
#define BARE_SOIL_ALBEDO 0.2	    /* albedo for bare soil */
#define RESID_MOIST      0.0        /* define residual moisture content 
				       of soil column */
#define ice_density      917.	    /* density of ice (kg/m^3) */
#define T_lapse          6.5        /* tempreature lapse rate of US Std 
				       Atmos in C/km */
#define von_K        0.40	/* Von Karmin constant for evapotranspiration */
#define KELVIN       273.15	/* conversion factor C to K */
#define STEFAN_B     5.6696e-8	/* stefan-boltzmann const in unit W/m^2/K^4 */
#define Lf           3.337e5	/* Latent heat of freezing (J/kg) at 0C */
#define RHO_W        1000.0	/* Density of water (kg/m^3) at 0C */
#define Cp           1004.0	/* Specific heat at constant pressure of air 
				   (J/deg/K) */
#define CH_ICE       2100.0e3	/* Volumetric heat capacity (J/(m3*C)) of ice */

#define SECPHOUR     3600	/* seconds per hour */
#define SNOW_DT       5.0	/* Used to bracket snow surface temperatures
				   while computing the snow surface energy 
				   balance (C) */
#define SURF_DT       1.0	/* Used to bracket soil surface temperatures 
                                   while computing energy balance (C) */
#define SOIL_DT       0.25      /* Used to bracket soil temperatures while
                                   solving the soil thermal flux (C) */
#define HOURSPERDAY   24        /* number of hours per day */
#define HOURSPERYEAR  24*365    /* number of hours per year */

/***** Physical Constraints *****/
#define MINSOILDEPTH 0.001	/* minimum layer depth with which model can
					work (m) */
#define STORM_THRES  0.001      /* thresehold at which a new storm is 
				   decalred */
/***** Define Boolean Values *****/
#ifndef FALSE
#define FALSE 0
#define TRUE !FALSE
#endif

#ifndef WET
#define WET 0
#define DRY 1
#endif

#ifndef SNOW
#define RAIN 0
#define SNOW 1
#endif

/***** Baseflow parametrizations *****/
#define ARNO        0
#define NIJSSEN2001 1

/***** Output BINARY format types *****/
#define OUT_TYPE_DEFAULT 0 /* Default data type */
#define OUT_TYPE_CHAR    1 /* char */
#define OUT_TYPE_SINT    2 /* short int */
#define OUT_TYPE_USINT   3 /* unsigned short int */
#define OUT_TYPE_INT     4 /* int */
#define OUT_TYPE_FLOAT   5 /* single-precision floating point */
#define OUT_TYPE_DOUBLE  6 /* double-precision floating point */

#define min(a,b) (a < b) ? a : b
#define max(a,b) (a > b) ? a : b

#include <user_def.h>
#include <snow.h>

#define DAYS_PER_YEAR 365.
#define DtoR 0.017453293	/* degrees to radians */
#ifndef PI
#define PI 3.1415927
#endif
#define STEFAN 5.6696e-8	/* Stefan boltzmann constant */
#define SOLAR_CONSTANT 1400.0	/* Solar constant in W/m^2 */
#define SEC_PER_DAY 86400.	/* seconds per day */

/* define constants for saturated vapor pressure curve (kPa) */
#define A_SVP 0.61078
#define B_SVP 17.269
#define C_SVP 237.3

/* define constants for penman evaporation */
#define CP_PM 1013		/* specific heat of moist air J/kg/C 
				   (Handbook of Hydrology) */
#define PS_PM 101300		/* sea level air pressure in Pa */
#define LAPSE_PM -0.006		/* environmental lapse rate in C/m */

/***** Codes for displaying version information *****/
#define DISP_VERSION 1
#define DISP_COMPILE_TIME 2
#define DISP_ALL 3

/***** VIC model version *****/
extern char *version;

/* global variables */
extern int NR;			/* array index for atmos struct that indicates
				   the model step avarage or sum */
extern int NF;			/* array index loop counter limit for atmos
				   struct that indicates the SNOW_STEP values */

/***** Data Structures *****/

/** file structures **/
typedef struct {
  FILE *forcing[2];     /* atmospheric forcing data files */
  FILE *globalparam;    /* global parameters file */
  FILE *init_snow;      /* snowpack initialization file */
  FILE *snowband;       /* snow elevation band data file */
  FILE *soilparam;      /* soil parameters for all grid cells */
  FILE *veglib;         /* vegetation parameters for all vege types */
  FILE *vegparam;       /* fractional coverage info for grid cell */
  FILE *statefile;      /* initial model state file */
} infiles_struct;

typedef struct {
  FILE *statefile;      /* output model state file */
} outfiles_struct;

typedef struct {
  char  forcing[2][MAXSTRING];  /* atmospheric forcing data file names */
  char  global[MAXSTRING];      /* global control file name */
  char  init_state[MAXSTRING];  /* initial model state file name */
  char  result_dir[MAXSTRING];  /* directory where results will be written */
  char  snow_band[MAXSTRING];   /* snow band parameter file name */
  char  soil[MAXSTRING];        /* soil parameter file name, or name of 
				   file that has a list of all soil 
				   ARC/INFO files */
  char  soil_dir[MAXSTRING];    /* directory from which to read ARC/INFO 
				   soil files */
  char  veg[MAXSTRING];         /* vegetation grid coverage file */
  char  veglib[MAXSTRING];      /* vegetation parameter library file */
} filenames_struct;

typedef struct {

  // simulation modes
  int    AboveTreelineVeg; // Default veg type to use above treeline
                           // Negative number indicates bare soil.
  char   COMPUTE_TREELINE; // TRUE = Determine treeline and exclude overstory
                           // vegetation from higher elevations
  char   CORRPREC;       /* TRUE = correct precipitation for gage undercatch */
  char   DIST_PRCP;      /* TRUE = Use distributed precipitation model */
  char   FROZEN_SOIL;    /* TRUE = Use frozen soils code */
  char   FULL_ENERGY;    /* TRUE = Use full energy code */
  char   GRND_FLUX;      /* TRUE = compute ground heat flux and energy
                            balance */
  float  MIN_WIND_SPEED; /* Minimum wind speed in m/s that can be used by
                            the model. **/
  char   MOISTFRACT;     /* TRUE = output soil moisture as moisture content */
  int    Nlayer;         /* Number of layers in model */
  int    Nnode;          /* Number of soil thermal nodes in the model */
  char   NOFLUX;         /* TRUE = Use no flux lower bondary when computing
                            soil thermal fluxes */
  char   QUICK_FLUX;     /* TRUE = Use Liang et al., 1999 formulation for
                            ground heat flux, if FALSE use explicit finite
                            difference method */
  int    SNOW_BAND;      /* Number of elevation bands over which to solve the
                            snow model */
  int    SNOW_STEP;      /* Time step in hours to use when solving the
                            snow model */

  // input options
  char   ARC_SOIL;       /* TRUE = use ARC/INFO gridded ASCII files for soil
                            parameters*/
  char   BASEFLOW;       /* ARNO: read Ds, Dm, Ws, c; NIJSSEN2001: read d1, d2, d3, d4 */
  int    GRID_DECIMAL;   /* Number of decimal places in grid file extensions */
  char   GLOBAL_LAI;     /* TRUE = read LAI values for each vegetation type
                            from the veg param file */
  char   JULY_TAVG_SUPPLIED; /* If TRUE and COMPUTE_TREELINE is also true,
                                then average July air temperature will be read
                                from soil file and used in calculating treeline */
  int    ROOT_ZONES;     /* Number of root zones used in simulation */

  // state options
  char   BINARY_STATE_FILE; /* TRUE = model state file is binary (default) */
  char   INIT_STATE;     /* TRUE = initialize model state from file */
  char   SAVE_STATE;      /* TRUE = save state file */

  // output options
  char   ALMA_OUTPUT;    /* TRUE = output variables are in ALMA-compliant units */
  char   BINARY_OUTPUT;  /* TRUE = output files are in binary, not ASCII */
  char   COMPRESS;       /* TRUE = Compress all output files */
  int    Noutfiles;      /* Number of output files (not including state files) */
  float  PREC_EXPT;      /* Exponential that controls the fraction of a
                            grid cell that receives rain during a storm
                            of given intensity */
  char   PRT_SNOW_BAND;  /* TRUE = print snow parameters for each snow band */

} option_struct;

#if LINK_DEBUG

typedef struct {
  FILE    *fg_balance;
  FILE    *fg_energy;
  FILE    *fg_grid;
  FILE    *fg_kappa;
  FILE    *fg_modelstep_atmos;
  FILE    *fg_moist;
  FILE    *fg_snow;
  FILE    *fg_snowstep_atmos;
  FILE    *fg_temp;
  char     DEBUG;
  char     PRT_ATMOS;
  char     PRT_BALANCE;
  char     PRT_FLUX;
  char     PRT_GLOBAL;
  char     PRT_GRID;
  char     PRT_KAPPA;
  char     PRT_MOIST;
  char     PRT_SNOW;
  char     PRT_SOIL;
  char     PRT_TEMP;
  char     PRT_VAR;
  char     PRT_VEGE;
  char     debug_dir[512];
  double **inflow[2];
  double **outflow[2];
  double **store_moist[2];
} debug_struct;

#endif

/*******************************************************
  Stores forcing file input information.
*******************************************************/
typedef struct {
  char    SIGNED;
  int     SUPPLIED;
  double  multiplier;
} force_type_struct;

/******************************************************************
  This structure records the parameters set by the forcing file
  input routines.  Those filled, are used to estimate the paramters
  needed for the model run in initialize_atmos.c.
  ******************************************************************/
typedef struct {
  force_type_struct TYPE[N_FORCING_TYPES];
  int  FORCE_DT[2];     /* forcing file time step */
  int  FORCE_ENDIAN[2]; /* endian-ness of input file, used for
			   DAILY_BINARY format */
  int  FORCE_FORMAT[2]; /* ASCII or BINARY */
  int  FORCE_INDEX[2][N_FORCING_TYPES];
  int  N_TYPES[2];
} param_set_struct;

/*******************************************************
  This structure stores all model run global parameters.
  *******************************************************/
typedef struct {
  char   statename[MAXSTRING];  /* name of file in which to store model state */
  double MAX_SNOW_TEMP; /* maximum temperature at which snow can fall (C) */
  double MIN_RAIN_TEMP; /* minimum temperature at which rain can fall (C) */
  double measure_h;  /* height of measurements (m) */
  double wind_h;     /* height of wind measurements (m) */ 
  float  resolution; /* Model resolution (degrees) */
  int    dt;         /* Time step in hours (24/dt must be an integer) */
  int    endday;     /* Last day of model simulation */
  int    endmonth;   /* Last month of model simulation */
  int    endyear;    /* Last year of model simulation */
  int    forceday[2];   /* day forcing files starts */
  int    forcehour[2];  /* hour forcing files starts */
  int    forcemonth[2]; /* month forcing files starts */
  int    forceskip[2];  /* number of model time steps to skip at the start of
			the forcing file */
  int    forceyear[2];  /* year forcing files start */
  int    nrecs;      /* Number of time steps simulated */
  int    skipyear;   /* Number of years to skip before writing output data */
  int    startday;   /* Starting day of the simulation */
  int    starthour;  /* Starting hour of the simulation */
  int    startmonth; /* Starting month of the simulation */
  int    startyear;  /* Starting year of the simulation */
  int    stateday;   /* Day of the simulation at which to save model state */
  int    statemonth; /* Month of the simulation at which to save model state */
  int    stateyear;  /* Year of the simulation at which to save model state */
} global_param_struct;

/***********************************************************
  This structure stores the soil parameters for a grid cell.
  ***********************************************************/
typedef struct {
  int      FS_ACTIVE;                 /* if TRUE frozen soil algorithm is 
					 active in current grid cell */
  double   Ds;                        /* fraction of maximum subsurface flow 
					 rate */
  double   Dsmax;                     /* maximum subsurface flow rate 
					 (mm/day) */
  double   Ksat[MAX_LAYERS];          /* saturated hydraulic  conductivity 
					 (mm/day) */
  double   Wcr[MAX_LAYERS];           /* critical moisture level for soil 
					 layer, evaporation is no longer 
					 affected moisture stress in the 
					 soil (mm) */
  double   Wpwp[MAX_LAYERS];          /* soil moisture content at permanent 
					 wilting point (mm) */
  double   Ws;                        /* fraction of maximum soil moisture */
  double   alpha[MAX_NODES];          /* thermal solution constant */
  double   annual_prec;               /* annual average precipitation (mm) */
  double   avg_temp;                  /* average soil temperature (C) */
  double   avgJulyAirTemp;            /* Average July air temperature (C) */
  double   b_infilt;                  /* infiltration parameter */
  double   beta[MAX_NODES];           /* thermal solution constant */
  double   bubble[MAX_LAYERS];        /* bubbling pressure, HBH 5.15 (cm) */
  double   bubble_node[MAX_NODES];    /* bubbling pressure (cm) */
  double   bulk_density[MAX_LAYERS];  /* soil bulk density (kg/m^3) */
  double   c;                         /* exponent */
  double   depth[MAX_LAYERS];         /* thickness of each soil moisture 
					 layer (m) */
  double   dp;                        /* soil thermal damping depth (m) */
  double   dz_node[MAX_NODES];        /* thermal node thickness (m) */
  double   expt[MAX_LAYERS];          /* pore-size distribution per layer, 
					 HBH 5.15 */
  double   expt_node[MAX_NODES];      /* pore-size distribution per node */
  double   gamma[MAX_NODES];          /* thermal solution constant */
  double   init_moist[MAX_LAYERS];    /* initial layer moisture level (mm) */
  double   max_infil;                 /* maximum infiltration rate */
  double   max_moist[MAX_LAYERS];     /* maximum moisture content (mm) per 
					 layer */
  double   max_moist_node[MAX_NODES]; /* maximum moisture content (mm/mm) per 
					 node */
  double   phi_s[MAX_LAYERS];         /* soil moisture diffusion parameter 
					 (mm/mm) */
  double   porosity[MAX_LAYERS];      /* porosity (fraction) */
  double   quartz[MAX_LAYERS];        /* quartz content of soil (fraction) */
  double   resid_moist[MAX_LAYERS];   /* residual moisture content of soil 
					 layer */
  double   rough;                     /* soil surface roughness (m) */
  double   snow_rough;                /* snow surface roughness (m) */
  double   soil_density[MAX_LAYERS];  /* soil partical density (kg/m^3) */
  double  *AreaFract;                 /* Fraction of grid cell included in 
					 each elevation band */
  double  *Pfactor;                   /* Change in Precipitation due to 
					 elevation (fract) */
  double  *Tfactor;                   /* Change in temperature due to 
					 elevation (C) */
  char    *AboveTreeLine;             // Flag to indicate if band is above 
                                      // the treeline
#if QUICK_FS
  double **ufwc_table_layer[MAX_LAYERS];
  double **ufwc_table_node[MAX_NODES]; 
#endif
  float    elevation;                 /* grid cell elevation (m) */
  float    lat;                       /* grid cell central latitude */
  float    lng;                       /* grid cell central longitude */
  float    time_zone_lng;             /* central meridian of the time zone */
  float  **layer_node_fract;          /* fraction of all nodes within each 
					 layer */
  int      gridcel;                   /* grid cell number */
} soil_con_struct;

/*******************************************************************
  This structure stores information about the vegetation coverage of
  the current grid cell.
  *******************************************************************/
typedef struct {
  double  Cv;               /* fraction of vegetation coverage */ 
  double  Cv_sum;           /* total fraction of vegetation coverage */
  float   root[MAX_LAYERS]; /* percent of roots in each soil layer (fraction) */
  float  *zone_depth;       /* depth of root zone */
  float  *zone_fract;       /* fraction of roots within root zone */
  int     veg_class;        /* vegetation class reference number */
  int     vegetat_type_num; /* number of vegetation types in the grid cell */
} veg_con_struct;

/******************************************************************
  This structure stores parameters for individual vegetation types.
  ******************************************************************/
typedef struct {
  char   overstory;        /* TRUE = overstory present, important for snow 
			      accumulation in canopy */
  double LAI[12];          /* monthly leaf area index */
  double Wdmax[12];        /* maximum monthly dew holding capacity (mm) */
  double albedo[12];       /* vegetation albedo (added for full energy) 
			      (fraction) */
  double displacement[12]; /* vegetation displacement height (m) */
  double emissivity[12];   /* vegetation emissivity (fraction) */
  double rad_atten;        /* radiation attenuation due to canopy, 
			      default = 0.5 (N/A) */
  double rarc;             /* architectural resistance (s/m) */
  double rmin;             /* minimum stomatal resistance (s/m) */
  double roughness[12];    /* vegetation roughness length (m) */
  double trunk_ratio;      /* ratio of trunk height to tree height, 
			      default = 0.2 (fraction) */
  double wind_atten;       /* wind attenuation through canopy, 
			      default = 0.5 (N/A) */
  double wind_h;           /* height at which wind is measured (m) */
  float  RGL;              /* Value of solar radiation below which there 
			      will be no transpiration (ranges from 
			      ~30 W/m^2 for trees to ~100 W/m^2 for crops) */
  int    veg_class;        /* vegetation class reference number */
} veg_lib_struct;

/***************************************************************************
   This structure stores the atmospheric forcing data for each model time 
   step for a single grid cell.  Each array stores the values for the 
   SNOW_STEPs during the current model step and the value for the entire model
   step.  The latter is referred to by array[NR].  Looping over the SNOW_STEPs
   is done by for (i = 0; i < NF; i++) 
***************************************************************************/
typedef struct {
  char   snowflag[25];  /* TRUE if there is snowfall in any of the snow 
			   bands during the timestep, FALSE otherwise*/
  double air_temp[25];  /* air temperature (C) */
  double density[25];   /* atmospheric density (kg/m^3) */
  double longwave[25];  /* incoming longwave radiation (W/m^2) (net incoming 
			   longwave for water balance model) */
  double out_prec;      /* Total precipitation for time step - accounts
			   for corrected precipitation totals */
  double prec[25];      /* average precipitation in grid cell (mm) */
  double pressure[25];  /* atmospheric pressure (kPa) */
  double shortwave[25]; /* incoming shortwave radiation (W/m^2) */
  double vp[25];        /* atmospheric vapor pressure (kPa) */
  double vpd[25];       /* atmospheric vapor pressure deficit (kPa) */
  double wind[25];      /* wind speed (m/s) */
} atmos_data_struct;

/*************************************************************************
  This structure stores information about the time and date of the current
  time step.
  *************************************************************************/
typedef struct {
  int day;                      /* current day */
  int day_in_year;              /* julian day in year */
  int hour;                     /* beginning of current hour */
  int month;                    /* current month */
  int year;                     /* current year */
} dmy_struct;			/* array of length nrec created */

/***************************************************************
  This structure stores all soil variables for each layer in the
  soil column.
  ***************************************************************/
typedef struct {
  double Cs;                /* average volumetric heat capacity of the 
			       current layer (J/m^3/K) */
  double T;                 /* temperature of the unfrozen sublayer (C) */
  double evap;              /* evapotranspiration from soil layer (mm) */
  double ice;               /* ice content of the frozen sublayer (mm) */
  double kappa;             /* average thermal conductivity of the current 
			       layer (W/m/K) */
  double moist;             /* moisture content of the unfrozen sublayer 
			       (mm) */
  double phi;               /* moisture diffusion parameter */
} layer_data_struct;

/******************************************************************
  This structure stores soil variables for the complete soil column 
  for each grid cell.
  ******************************************************************/
typedef struct {
  double aero_resist[3];               /* aerodynamic resistance (s/m) 
					  [0] = over vegetation or bare soil 
					  [1] = over snow-filled overstory
					  [2] = over snow */
  double aero_resist_used;             /* The (stability-corrected) aerodynamic
                                          resistance (s/m) that was actually used
                                          in flux calculations.  For cases in which
                                          a cell uses 2 different resistances for
                                          flux computations in the same time step
                                          (i.e. cell contains overstory and snow
                                          is present on the ground), aero_resist_used
                                          will contain the snow pack's resistance. */
  double baseflow;                     /* baseflow from current cell (mm/TS) */
  double inflow;                       /* moisture that reaches the top of 
					  the soil column (mm) */
  double runoff;                       /* runoff from current cell (mm/TS) */
  layer_data_struct layer[MAX_LAYERS]; /* structure containing soil variables 
					  for each layer (see above) */
} cell_data_struct;

/***********************************************************************
  This structure stores energy balance components, and variables used to
  solve the thermal fluxes through the soil column.
  ***********************************************************************/
typedef struct {
  char    frozen;                /* TRUE = frozen soil present */
  double  Cs[2];                 /* heat capacity for top two layers 
				    (J/m^3/K) */
  double  Cs_node[MAX_NODES];    /* heat capacity of the soil thermal nodes 
				    (J/m^3/K) */
  double  T[MAX_NODES];          /* thermal node temperatures (C) */
  double  Trad[2];               /* surface temperature of energy balance 
				    (C) */
  double  advection;             /* advective flux (Wm-2) */
  double  albedo;                /* surface albedo (fraction) */
  double  deltaCC;               /* change in snow heat storage (Wm-2) */
  double  deltaH;                /* change in soil heat storage (Wm-2) */
  double  error;                 /* energy balance error (W/m^2) */
  double  fdepth[MAX_FRONTS];    /* all simulated freezing front depths */
  double  grnd_flux;             /* ground heat flux (Wm-2) */
  double  ice[MAX_NODES];        /* thermal node ice content */
  double  kappa[2];              /* soil thermal conductivity for top two 
				    layers (W/m/K) */
  double  kappa_node[MAX_NODES]; /* thermal conductivity of the soil thermal 
				    nodes (W/m/K) */
  double  latent;                /* net latent heat flux (Wm-2) */
  double  longwave;              /* net longwave flux (Wm-2) */
  double  moist[MAX_NODES];      /* thermal node moisture content */
  double  refreeze_energy;       /* energy used to refreeze the snowpack 
				    (Wm-2) */
  double  sensible;              /* net sensible heat flux (Wm-2) */
  double  shortwave;             /* incoming shortwave heat (Wm-2) */
  double  snow_flux;             /* thermal flux through the snow pack 
				    (Wm-2) */
  double  tdepth[MAX_FRONTS];    /* all simulated thawing front depths */
  double  unfrozen;              /* frozen layer water content that is 
				    unfrozen */
  int     Nfrost;                /* number of simulated freezing fronts */
  int     Nthaw;                 /* number of simulated thawing fronts */
  int     T1_index;              /* soil node at the bottom of the top layer */
} energy_bal_struct;

/***********************************************************************
  This structure stores vegetation variables for each vegetation type in 
  a grid cell.
  ***********************************************************************/
typedef struct {
  double canopyevap;		/* evaporation from canopy (mm/TS) */
  double throughfall;		/* water that reaches the ground through 
                                   the canopy (mm/TS) */
  double Wdew;			/* dew trapped on vegetation (mm) */
} veg_var_struct;

/************************************************************************
  This structure stores snow pack variables needed to run the snow model.
  ************************************************************************/
typedef struct {
  char   MELTING;           /* flag indicating that snowpack melted 
			       previously */
  int    snow;              /* TRUE = snow, FALSE = no snow */
  double Qnet;              /* New energy at snowpack surface */
  double albedo;            /* snow surface albedo (fraction) */
  double canopy_vapor_flux; /* depth of water evaporation, sublimation, or 
			       condensation from intercepted snow (m) */
  double coldcontent;       /* cold content of snow pack */
  double coverage;          /* fraction of snow band that is covered with 
			       snow */
  double density;           /* snow density (kg/m^3) */
  double depth;             /* snow depth (m) */
  double mass_error;        /* snow mass balance error */
  double pack_temp;         /* depth averaged temperature of the snowpack 
			       (C) */
  double pack_water;        /* liquid water content of the snow pack (m) */
  double snow_canopy;       /* amount of snow on canopy (m) */
  double surf_temp;         /* depth averaged temperature of the snow pack 
			       surface layer (C) */
  double surf_water;        /* liquid water content of the surface layer (m) */
  double swq;               /* snow water equivalent of the entire pack (m) */
  double tmp_int_storage;   /* temporary canopy storage, used in snow_canopy */
  double vapor_flux;        /* depth of water evaporation, sublimation, or 
			       condensation from snow pack (m) */
  int    last_snow;         /* time steps since last snowfall */
} snow_data_struct;	    /* an array of size Nrec */

/*****************************************************************
  This structure stores all variables needed to solve, or save 
  solututions for all versions of this model.  Vegetation and soil
  variables are created for both wet and dry fractions of the grid
  cell (for use with the distributed precipitation model).
*****************************************************************/
typedef struct {
  cell_data_struct  **cell[2];    /* Stores soil layer variables (wet and 
				     dry) */
  double             *mu;         /* fraction of grid cell that receives 
				     precipitation */
  energy_bal_struct **energy;     /* Stores energy balance variables */
  snow_data_struct  **snow;       /* Stores snow variables */
  veg_var_struct    **veg_var[2]; /* Stores vegetation variables (wet and 
				     dry) */
} dist_prcp_struct;

/*******************************************************
  This structure stores output information for one variable.
  *******************************************************/
typedef struct {
  char		varname[20]; /* name of variable */
  int		write;       /* FALSE = don't write; TRUE = write */
  char		format[10];  /* format, when written to an ascii file;
		                should match the desired fprintf format specifier, e.g. %.4f */
  int		type;        /* type, when written to a binary file;
		                OUT_TYPE_USINT  = unsigned short int
		                OUT_TYPE_SINT   = short int
		                OUT_TYPE_FLOAT  = single precision floating point
		                OUT_TYPE_DOUBLE = double precision floating point */
  float		mult;        /* multiplier, when written to a binary file */
  int		nelem;       /* number of data values */
  double	*data;       /* array of data values */
} out_data_struct;

/*******************************************************
  This structure stores output information for one output file.
  *******************************************************/
typedef struct {
  char		prefix[20];  /* prefix of the file name, e.g. "fluxes" */
  char		filename[MAXSTRING]; /* complete file name */
  FILE		*fh;         /* filehandle */
  int		nvars;       /* number of variables to store in the file */
  int		*varid;      /* id numbers of the variables to store in the file
		                (a variable's id number is its index in the out_data array).
		                The order of the id numbers in the varid array
		                is the order in which the variables will be written. */
} out_data_file_struct;

/********************************************************
  This structure holds all variables needed for the error
  handling routines.
  ********************************************************/
typedef struct {
  atmos_data_struct *atmos;
  double             dt;
  energy_bal_struct *energy;
  infiles_struct     infp;
  int                rec;
  out_data_struct   *out_data;
  outfiles_struct    outfp;
  out_data_file_struct    *out_data_files;
  snow_data_struct  *snow;
  soil_con_struct    soil_con;
  veg_con_struct    *veg_con;
  veg_var_struct    *veg_var;
} Error_struct;

