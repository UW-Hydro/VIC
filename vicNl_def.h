/***** Model Constants *****/
#define MAXSTRING    512
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
#define N_FORCING_TYPES 14
#define AIR_TEMP  0  /* air temperature per time step (C) */
#define ALBEDO    1  /* surface albedo (fraction) */
#define DENSITY   2  /* atmospheric density (kg/m^3) */
#define LONGWAVE  3  /* incoming longwave radiation (W/m^2) */
#define PREC      4  /* precipitation (mm) */
#define PRESSURE  5  /* atmospheric pressure (kPa) */
#define SHORTWAVE 6  /* incoming shortwave (W/m^2) */
#define TMAX      7  /* maximum daily temperature (C) */
#define TMIN      8  /* minimum daily temperature (C) */
#define TSKC      9  /* cloud cover (fraction) */
#define VP        10 /* vapor pressure (kPa) */
#define WIND      11 /* wind speed (m/s) */
#define WINDE     12
#define SKIP      13 /* place holder for unused data columns */
		 
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
#define Cp           1010.0	/* Specific heat at constant pressure of air 
				   (J/deg/K) (H.B.H. p.4.7)*/
#define CH_ICE       2100.0e3	/* Volumetric heat capacity (J/(m3*C)) of ice */
#define CH_WATER     4186.8e3   // volumetric heat capacity of water
#define GRAMSPKG     1000.      // convert grams to kilograms
#define K_SNOW       2.9302e-6  // conductivity of snow (W/mK)

#define EPS          0.62196351 // (MW/MA) unitless ratio of molec weights */
#define JOULESPCAL   4.1868     /* Joules per calorie */
#define SECPHOUR     3600	/* seconds per hour */
#define CANOPY_DT    1.0	/* Used to bracket canopy air temperatures 
                                   while computing energy balance (C) */
#define CANOPY_VP    25.0	/* Used to bracket canopy vapor pressures 
                                   while computing moisture balance (Pa) */
#define SNOW_DT      5.0	/* Used to bracket snow surface temperatures
				   while computing the snow surface energy 
				   balance (C) */
#define SURF_DT      1.0	/* Used to bracket soil surface temperatures 
                                   while computing energy balance (C) */
#define SOIL_DT      0.25       /* Used to bracket soil temperatures while
                                   solving the soil thermal flux (C) */
#define HOURSPERDAY  24         /* number of hours per day */
#define HOURSPERYEAR 24*365     /* number of hours per year */
#define G            9.81       // gravity

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
#define kPa2Pa 1000.            /* converts kPa to Pa */

/* define constants for saturated vapor pressure curve (kPa) */
#define A_SVP 0.61078
#define B_SVP 17.269
#define C_SVP 237.3

/* define constants for penman evaporation */
#define CP_PM 1013		/* specific heat of moist air J/kg/C 
				   (Handbook of Hydrology) */
#define PS_PM 101300		/* sea level air pressure in Pa */
#define LAPSE_PM -0.006		/* environmental lapse rate in C/m */

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
  FILE *init_soil;      /* soil temp and mosit initialization file */
#if LAKE_MODEL
  FILE *lakeparam;      /* lake parameter file */
#endif // LAKE_MODEL
  FILE *snowband;       /* snow elevation band data file */
  FILE *soilparam;      /* soil parameters for all grid cells */
  FILE *veglib;         /* vegetation parameters for all vege types */
  FILE *vegparam;       /* fractional coverage info for grid cell */
  FILE *statefile;      /* initial model state file */
} infiles_struct;

typedef struct {
  FILE *fdepth;
  FILE *fluxes;
#if LAKE_MODEL
  FILE *lake;
#endif // LAKE_MODEL
  FILE *snow;
  FILE *snowband;
#if SAVE_STATE
  FILE *statefile;
#endif // SAVE_STATE
} outfiles_struct;

typedef struct {
  char  fdepth[MAXSTRING];      /* frozen soils depth (output) */
  char  fluxes[MAXSTRING];      /* grid cell surface fluxes (output) */
  char  forcing[2][MAXSTRING];  /* atmospheric forcing data file names */
  char  global[MAXSTRING];      /* global control file name */
  char  init_state[MAXSTRING];  /* initial model state file name */
#if LAKE_MODEL
  char  lakeparam[MAXSTRING];   /* lake model constants file */
#endif // LAKE_MODEL
  char  result_dir[MAXSTRING];  /* directory where results will be written */
  char  snow[MAXSTRING];        /* snow pack depth and swq (output) */
  char  snow_band[MAXSTRING];   /* snow band parameter file name */
  char  snowband[MAXSTRING];    /* snow band pack depth and swq (output) */
  char  soil[MAXSTRING];        /* soil parameter file name, or name of 
				   file that has a list of all aoil 
				   ARC/INFO files */
  char  soil_dir[MAXSTRING];    /* directory from which to read ARC/INFO 
				   soil files */
  char  veg[MAXSTRING];         /* vegetation grid coverage file */
  char  veglib[MAXSTRING];      /* vegetation parameter library file */
} filenames_struct;

typedef struct {
  char   ARC_SOIL;       /* TRUE = use ARC/INFO gridded ASCII files for soil 
			    parameters*/
  char   BINARY_OUTPUT;  /* TRUE = output files are in binary, not ASCII */
  char   COMPRESS;       /* TRUE = Compress all output files */
  char   CORRPREC;       /* TRUE = correct precipitation for gage undercatch */
  char   DIST_PRCP;      /* TRUE = Use distributed precipitation model */
  char   FROZEN_SOIL;    /* TRUE = Use frozen soils code */
  char   FULL_ENERGY;    /* TRUE = Use full energy code */
  char   GLOBAL_LAI;     /* TRUE = read LAI values for each vegetation type
			    from the veg param file */
  char   GRND_FLUX;      /* TRUE = compute ground heat flux and energy 
			    balance */
  char   INIT_STATE;     /* TRUE = initialize model state from file */
  char   NEW_ARNO_TYPE;  /* FALSE = Ds, Dm, Ws, c  TRUE = d1, d2, d3, d4 */
  char   LAKE_PROFILE;   /* TRUE = user-specified lake/area profile */
  char   LAKES;          /* TRUE = use lake energy code */
  char   MOISTFRACT;     /* TRUE = output soil moisture as moisture content */
  char   NOFLUX;         /* TRUE = Use no flux lower bondary when computing 
			    soil thermal fluxes */
  char   PRT_SNOW_BAND;  /* TRUE = print snow parameters for each snow band */
  char   QUICK_FLUX;     /* TRUE = Use Liang et al., 1999 formulation for
			    ground heat flux, if FALSE use explicit finite
			    difference method */
  char   QUICK_SOLVE;    /* TRUE = Use Liang et al., 1999 formulation for 
			    iteration, but explicit finite difference
			    method for final step. */
  char   BLOWING;        /* TRUE = calculate sublimation from blowing snow */
  float  MIN_WIND_SPEED; /* Minimum wind speed in m/s that can be used by 
			    the model. **/
  float  PREC_EXPT;      /* Exponential that controls the fraction of a
			    grid cell that receives rain during a storm
			    of given intensity */
  int    GRID_DECIMAL;   /* Number of decimal places in grid file extensions */
#if LAKE_MODEL
  int    Nlakenode;      /* Number of lake thermal nodes in the model. */
#endif // LAKE_MODEL
  int    Nlayer;         /* Number of layers in model */
  int    Nnode;          /* Number of soil thermal nodes in the model */
  int    ROOT_ZONES;     /* Number of root zones used in simulation */
  int    SNOW_BAND;      /* Number of elevation bands over which to solve the 
			    snow model */
  int    SNOW_STEP;      /* Time step in hours to use when solving the 
			    snow model */
} option_struct;

#if LINK_DEBUG

typedef struct {
  FILE    *fg_balance;
  FILE    *fg_energy;
  FILE    *fg_grid;
  FILE    *fg_kappa;
  FILE    *fg_lake;
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
  char     PRT_LAKE;
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
#if SAVE_STATE
  char   statename[MAXSTRING];  /* name of file in which to store model state */
#endif
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
#if SAVE_STATE
  int    stateday;   /* Day of the simulation at which to save model state */
  int    statemonth; /* Month of the simulation at which to save model state */
  int    stateyear;  /* Year of the simulation at which to save model state */
#endif
} global_param_struct;

#if LAKE_MODEL
/******************************************************************
  This structure stores the lake/wetland parameters for a grid cell
  ******************************************************************/
typedef struct {
  double Cl[MAX_LAKE_NODES];      /* fractional lake coverage area */
  double z[MAX_LAKE_NODES];       /* Fixed elevation from bottom of each Cl. */  
  double b;                       /* Exponent controlling lake depth y=Ax^b. */
  double basin[MAX_LAKE_NODES];   /* Area of the basin at each node. */
  double cell_area;               /* area of grid cell */
  double depth_in;                /* initial lake depth */
  double eta_a;                   /* Decline of solar rad w/ depth. */
  double maxdepth;                /* Maximum lake depth. */
  double maxrate;                 
  double mindepth;                /* Minimum lake depth. */
  double maxvolume;
  float  bpercent;
  float  rpercent;
  int wetland_veg_class;
  int    gridcel;
  int    numnod;                  /* Maximum number of solution nodes. */
} lake_con_struct;

/*****************************************************************
  This structure stores the lake/wetland variables for a grid cell
  *****************************************************************/
typedef struct {
  /** Use MAX_LAKE_NODES **/
  double aero_resist;
  double baseflow_in;
  double baseflow_out;
  double density[MAX_LAKE_NODES];
  double evapw;
  double fraci;                   /* Fractional coverage of ice. */
  double hice;                    /* Height of lake ice. */ 
  double ldepth;
  double runoff_in;
  double runoff_out;
  double sarea;
  double sdepth;                  /* Depth of snow on top of ice. */
  double snowmlt;
  double surface[MAX_LAKE_NODES];
  double swe;
  double temp[MAX_LAKE_NODES];    /* Lake water temp. at each node (C). */
  double tempavg;
  double tempi;                   /* Lake ice temp (C). */
  double tp_in;                   /* Lake skin temperature (C). */
  double volume;
  double dz;                      /* Distance between each water layer. */
  int    activenod;
  int    mixmax;                  /* top depth (node #) of local 
                                     instability. */
} lake_var_struct;
#endif // LAKE_MODEL 

/***********************************************************
  This structure stores the soil parameters for a grid cell.
  ***********************************************************/
typedef struct {
  char     FS_ACTIVE;                 /* if TRUE frozen soil algorithm is 
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
  double   b_infilt;                  /* infiltration parameter */
  double   beta[MAX_NODES];           /* thermal solution constant */
  double   bubble[MAX_LAYERS];        /* bubbling pressure, HBH 5.15 (cm) */
  double   bubble_node[MAX_NODES];    /* bubbling pressure (cm) */
  double   bulk_density[MAX_LAYERS];  /* soil bulk density (kg/m^3) */
  double   c;                         /* exponent */
  double   depth[MAX_LAYERS];         /* thickness of each soil moisture 
					 layer (m) */
#if SPATIAL_SNOW
  double   depth_full_snow_cover;     // minimum depth for full snow cover
#endif // SPATIAL_SNOW
  double   dp;                        /* soil thermal damping depth (m) */
  double   dz_node[MAX_NODES];        /* thermal node thickness (m) */
  double   expt[MAX_LAYERS];          /* pore-size distribution per layer, 
					 HBH 5.15 */
  double   expt_node[MAX_NODES];      /* pore-size distribution per node */
#if SPATIAL_FROST
  double   frost_fract[FROST_SUBAREAS]; /* spatially distributed frost 
					   coverage fractions */
  double   frost_slope;               // slope of frost distribution
#endif // SPATIAL_FROST
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
  float   sigma_slope;      /* Std. deviation of terrain slope for each vegetation class. */
  float   lag_one;          /* Lag one gradient autocorrelation of terrain slope */
  float   fetch;            /* Average fetch length for each vegetation class. */
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
#if OPTIMIZE || LINK_DEBUG
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
#else
  char   *snowflag;  /* TRUE if there is snowfall in any of the snow 
			   bands during the timestep, FALSE otherwise*/
  double *air_temp;  /* air temperature (C) */
  double *density;   /* atmospheric density (kg/m^3) */
  double *longwave;  /* incoming longwave radiation (W/m^2) (net incoming 
			   longwave for water balance model) */
  double out_prec;      /* Total precipitation for time step - accounts
			   for corrected precipitation totals */
  double *prec;      /* average precipitation in grid cell (mm) */
  double *pressure;  /* atmospheric pressure (kPa) */
  double *shortwave; /* incoming shortwave radiation (W/m^2) */
  double *vp;        /* atmospheric vapor pressure (kPa) */
  double *vpd;       /* atmospheric vapor pressure deficit (kPa) */
  double *wind;      /* wind speed (m/s) */
#endif // OPTIMIZE or LINK_DEBUG
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
#if SPATIAL_FROST
  double ice[FROST_SUBAREAS]; /* ice content of the frozen sublayer (mm) */
#else
  double ice;               /* ice content of the frozen sublayer (mm) */
#endif
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
  double aero_resist[3];               /* aerodynamic resistane (s/m) 
					  [0] = over bare vegetation or soil 
					  [2] = over snow */
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
  double  AlbedoLake;            /* albedo of lake surface (fract) */
  double  AlbedoOver;            /* albedo of intercepted snow (fract) */
  double  AlbedoUnder;           /* surface albedo (fraction) */
  double  AtmosError;
  double  AtmosLatent;           /* latent heat exchange with atmosphere */
  double  AtmosLatentSub;        /* latent sub heat exchange with atmosphere */
  double  AtmosSensible;         /* sensible heat exchange with atmosphere */
  double  Cs[2];                 /* heat capacity for top two layers 
				    (J/m^3/K) */
  double  Cs_node[MAX_NODES];    /* heat capacity of the soil thermal nodes 
				    (J/m^3/K) */
  double  LongOverIn;            /* incoming longwave to overstory */
  double  LongUnderIn;           /* incoming longwave to understory */
  double  LongUnderOut;          /* outgoing longwave to understory */
  double  NetLongAtmos;          /* net longwave radiation to the atmosphere 
				    (W/m^2) */
  double  NetLongOver;           /* net longwave radiation from the canopy 
				    (W/m^2) */
  double  NetLongUnder;          /* net longwave radiation from the canopy 
				    (W/m^2) */
  double  NetShortAtmos;         /* net shortwave to the atmosphere */
  double  NetShortGrnd;          /* net shortwave penetrating snowpack */
  double  NetShortOver;          /* net shortwave radiation from the canopy 
				    (W/m^2) */
  double  NetShortUnder;         /* net shortwave radiation from the canopy 
				    (W/m^2) */
  double  ShortOverIn;           /* incoming shortwave to overstory */
  double  ShortUnderIn;          /* incoming shortwave to understory */
  double  T[MAX_NODES];          /* thermal node temperatures (C) */
  double  Tcanopy;               /* temperature of the canopy air */
  double  Tfoliage;              /* temperature of the overstory vegetation */
  double  Tsurf;                 /* temperature of the understory */
  double  Trad[2];               /* surface temperature of energy balance 
				    (C) */
  double  advected_sensible;     /* net sensible heat flux advected to 
				    snowpack (Wm-2) */
  double  advection;             /* advective flux (Wm-2) */
  double  canopy_advection;      /* advection heat flux from the canopy 
				    (W/m^2) */
  double  canopy_latent;         /* latent heat flux from the canopy (W/m^2) */
  double  canopy_latent_sub;     /* latent heat flux of sublimation 
				    from the canopy (W/m^2) */
  double  canopy_refreeze;       /* energy used to refreeze/melt canopy 
				    intercepted snow (W/m^2) */
  double  canopy_sensible;       /* sensible heat flux from canopy 
                                    interception (W/m^2) */
  double  deltaCC;               /* change in snow heat storage (Wm-2) */
  double  deltaH;                /* change in soil heat storage (Wm-2) */
  double  error;                 /* energy balance error (W/m^2) */
  double  fdepth[MAX_FRONTS];    /* all simulated freezing front depths */
  double  fusion;                /* energy used to freeze/thaw soil water */
  double  grnd_flux;             /* ground heat flux (Wm-2) */
  double  ice[MAX_NODES];        /* thermal node ice content */
  double  kappa[2];              /* soil thermal conductivity for top two 
				    layers (W/m/K) */
  double  kappa_node[MAX_NODES]; /* thermal conductivity of the soil thermal 
				    nodes (W/m/K) */
  double  latent;                /* net latent heat flux (Wm-2) */
  double  latent_sub;            /* net latent heat flux from snow (Wm-2) */
  double  longwave;              /* net longwave flux (Wm-2) */
  double  melt_energy;           /* energy used to reduce snow cover 
				    fraction (Wm-2) */
  double  moist[MAX_NODES];      /* thermal node moisture content */
  double  out_long_canopy;       /* outgoing longwave to canopy */
  double  out_long_surface;      /* outgoing longwave to surface */
  double  refreeze_energy;       /* energy used to refreeze the snowpack 
				    (Wm-2) */
  double  sensible;              /* net sensible heat flux (Wm-2) */
  double  shortwave;             /* net shortwave radiation (Wm-2) */
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
  char   snow;              /* TRUE = snow, FALSE = no snow */
  double Qnet;              /* New energy at snowpack surface */
  double albedo;            /* snow surface albedo (fraction) */
  double canopy_albedo;     /* albedo of the canopy (fract) */
  double canopy_vapor_flux; /* depth of water evaporation, sublimation, or 
			       condensation from intercepted snow (m) */
  double coldcontent;       /* cold content of snow pack */
  double coverage;          /* fraction of snow band that is covered with 
			       snow */
  double density;           /* snow density (kg/m^3) */
  double depth;             /* snow depth (m) */
  double mass_error;        /* snow mass balance error */
  double max_swq;           /* last maximum swq - used to determine coverage
			       fraction during current melt period (m) */
  double melt;              /* snowpack melt (mm) */
  double pack_temp;         /* depth averaged temperature of the snowpack 
			       (C) */
  double pack_water;        /* liquid water content of the snow pack (m) */
  double snow_canopy;       /* amount of snow on canopy (m) */
  double store_coverage;    /* stores coverage fraction covered by new 
			       snow (m) */
  double store_swq;         /* stores newly accumulated snow over an 
			       established snowpack melt distribution (m) */
  double surf_temp;         /* depth averaged temperature of the snow pack 
			       surface layer (C) */
  double surf_water;        /* liquid water content of the surface layer (m) */
  double swq;               /* snow water equivalent of the entire pack (m) */
  double swq_slope;         /* slope of uniform snow distribution (m/fract) */
  double tmp_int_storage;   /* temporary canopy storage, used in snow_canopy */
  double vapor_flux;        /* depth of water evaporation, sublimation, or 
			       condensation from snow pack (m) */
  double blowing_flux;      /* depth of sublimation from blowing snow (m) */
  double surface_flux;      /* depth of sublimation from blowing snow (m) */
  int    last_snow;         /* time steps since last snowfall */
  int    store_snow;        /* flag indicating whether or not new accumulation
			       is stored on top of an existing distribution */
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
#if LAKE_MODEL
  lake_var_struct     lake_var;   /* Stores lake/wetland variables */
#endif // LAKE_MODEL
  snow_data_struct  **snow;       /* Stores snow variables */
  veg_var_struct    **veg_var[2]; /* Stores vegetation variables (wet and 
				     dry) */
} dist_prcp_struct;

/*******************************************************
  This structure stores all variables needed for output.
  *******************************************************/
typedef struct {
  double Wdew;                         /* canopy interception of moisture */
  double advected_sensible[MAX_BANDS+1]; /* grid cell net sensible heat flux 
					  advected to snowpack (Wm-2) */
  double advection[MAX_BANDS+1];       /* grid cell advection (snow only) 
					  (Wm-2) */
  double aero_resist;                  /* grid cell mean aerodynamic 
					  resistence  [s/m] */
  double air_temp;                     /* grid cell air temperature */
  double albedo;                       /* grid cell mean albedo */ 
  double baseflow;                     /* baseflow out of the bottom layer */
  double bot_energy_error[2];
  double coverage[MAX_BANDS+1];        /* fractional coverage of grid cell 
					  with snow */
  double deltaCC[MAX_BANDS+1];         /* change of cold content in the 
					  snowpack [Wm-2] */
  double deltaH;                       /* grid cell change in heat storage 
					  (snow only) */
  double energy_error;                 /* energy balance error */
  double evap;                         /* grid cell evaporation */
  double evap_bare;                    /* grid cell net evaporation from 
					  bare ground */
  double evap_canop;                   /* grid cell net evaporation from 
					  canopy interception */
  double evap_lake;                    // evaporation from lake surface
  double evap_veg;                     /* grid cell net evapotraspiration 
					  from vegetation */
  double fdepth[MAX_FRONTS];           /* depth of all freezing fronts */
  double fusion;                       /* energy used to freeze/thaw soil 
					  water */
  double grnd_flux;                    /* grid cell ground flux */
  double ice[MAX_LAYERS];              /* frozen layer ice content */
  double in_long;                      /* grid cell net incoming longwave 
					  flux */
  double inflow;                       /* moisture that reaches the top of 
					  the soil column */
#if LAKE_MODEL
  double lake_depth;                   /* lake depth */
  double lake_ice;                     /* moisture stored as lake ice */
  double lake_ice_fract;               /* fractional coverage of lake ice */
  double lake_ice_height;              /* thickness of lake ice */
  double lake_ice_temp;                /* temperature of lake ice */
  double lake_moist;                   /* liquid water stored in lake */
  double lake_surf_area;               /* lake surface area */
  double lake_surf_temp;               /* lake surface temperature */
  double lake_volume;                  /* lake volume */
#endif // LAKE_MODEL
  double latent;                       /* grid cell net latent heat flux */
  double latent_sub[MAX_BANDS+1];      /* grid cell net latent heat flux */
  double melt[MAX_BANDS+1];            /* snowpack melt (mm) */
  double melt_energy[MAX_BANDS+1];     /* energy used to reduce snow 
					  coverage (Wm-2) */
  double moist[MAX_LAYERS];            /* current moisture in each layer */
  double net_long;                     /* grid cell net longwave flux */
  double net_short;                    /* grid cell net shortwave flux */
  double prec;                         /* incoming precipitation */
  double r_net;                        /* grid cell net radiation W/m^2 */
  double rad_temp;                     /* grid cell average radiative surface 
					  temperature */
  double refreeze_energy[MAX_BANDS+1]; /* energy used to refreeze snowpack 
					  [Wm-2] */
  double rel_humid;
  double runoff;                       /* runoff from the surface */
  double sensible;                     /* grid cell net sensible heat flux */
  double shortwave;                    /* grid cell incoming shortwave flux */
  double snow_canopy[MAX_BANDS+1];     /* snow captured by canopy (mm) */
  double snow_depth[MAX_BANDS+1];      /* snow depth (cm) */
  double snow_flux[MAX_BANDS+1];       /* energy flux through the snowpack 
					  [Wm-2] */
  double snow_pack_temp[MAX_BANDS+1];  /* snow pack surface layer temp (C) */
  double snow_surf_temp[MAX_BANDS+1];  /* snow pack surface layer temp (C) */
  double sub_canop;                    /* grid cell net sublimation from 
					  canopy interception */
  double sub_snow[MAX_VEG];            /* grid cell net sublimation from 
					  bare ground from snow pack */
  double sub_total;
  double sub_surface;
  double sub_blowing;
  double surf_cond;                    /* grid cell mean surface conductance 
					  [m/s] */
  double surf_frost_fract;             /* fraction of surface frozen (fract) */
  double surf_temp;                    /* grid cell average daily surface 
					  temperature */
  double swq[MAX_BANDS+1];             /* snow water equivalent (mm) */
  double tdepth[MAX_FRONTS];           /* depth of all thawing fronts */
  double wind;                         /* grid cell wind speed */
} out_data_struct;

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
  snow_data_struct  *snow;
  soil_con_struct    soil_con;
  veg_con_struct    *veg_con;
  veg_var_struct    *veg_var;
} Error_struct;

