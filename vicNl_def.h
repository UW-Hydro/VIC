/***** Physical Constants *****/
#define RESID_MOIST  0.0        /* define residual moisture content of soil column */
#define ice_density  917.	/* density of ice (kg/m^3) */
#define T_lapse      6.5        /* tempreature lapse rate of US Std Atmos in C/km */
#define von_K        0.40	/* Von Karmin constant for evapotranspiration */
#define KELVIN       273.15	/* conversion factor C to K */
#define STEFAN_B     5.6696e-8	/* stefan-boltzmann const in unit W/m^2/K^4 */
#define Lf           3.337e5	/* Latent heat of freezing (J/kg) at 0C */
#define RHO_W        1000.0	/* Density of water (kg/m^3) at 0C */
#define Cp           1004.0	/* Specific heat at constant pressure of air 
				   (J/deg/K) */
#define CH_ICE       2100.0e3	/* Volumetric heat capacity (J/(m3*C)) of ice */

#define SECPHOUR     3600	/* seconds per hour */
#define DELTAT       25.	/* Used in SensibleHeatFlux to bracket the 
                                   effective surface temperature (C) */
#define HOURSPERDAY  24         /* number of hours per day */

/***** Physical Constraints *****/
#define MINSOILDEPTH 0.001	/* minimum layer depth with which model can
					work (m) */
#define STORM_THRES  0.001      /* thresehold at which a new storm is 
				   decalred */
/***** Model Constants *****/
#define MAXSTRING    512
#define MINSTRING    20
#define HUGE_RESIST  1.e20	/* largest allowable double number */
#define SMALL        1.e-15	/* smallest allowable double number */

/***** Define Boolean Values *****/
#ifndef FALSE
#define FALSE 0
#define TRUE !FALSE
#endif

#ifndef WET
#define WET 0
#define DRY 1
#endif

#define min(a,b) (a < b) ? a : b
#define max(a,b) (a > b) ? a : b

#include <user_def.h>
#include <rad_and_vpd.h>
#include <snow.h>

/***** Data Structures *****/

/** file structures **/
typedef struct {
  FILE *forcing[2];     /* atmospheric forcing data files */
  FILE *soilparam;	/* soil parameters for all grid cells */
  FILE *veglib;		/* vegetation parameters for all vege types */
  FILE *vegparam;	/* fractional coverage info for grid cell */
  FILE *snowband;       /* snow elevation band data file */
  FILE *globalparam;	/* global parameters file */
  FILE *init_soil;	/* soil temp and mosit initialization file */
  FILE *init_snow;	/* snowpack initialization file */
} infiles_struct;

typedef struct {
  FILE *fdepth;
  FILE *fluxes;
  FILE *snow;
  FILE *snowband;
} outfiles_struct;

typedef struct {
  char  forcing[2][MAXSTRING];	/* atmospheric forcing data file names */
  char  global[MAXSTRING];      /* global control file name */
  char  soil[MAXSTRING];        /* soil parameter file name, or name of 
				   file that has a list of all aoil 
				   ARC/INFO files */
  char  soil_dir[MAXSTRING];    /* directory from which to read ARC/INFO 
				   soil files */
  char  veglib[MAXSTRING];	/* vegetation parameter library file */
  char  veg[MAXSTRING];		/* vegetation grid coverage file */
  char  snow_band[MAXSTRING];   /* snow band parameter file name */
  char  result_dir[MAXSTRING];  /* directory where results will be written */
  char  fluxes[MAXSTRING];	/* grid cell surface fluxes (output) */
  char  fdepth[MAXSTRING];	/* frozen soils depth (output) */
  char  snow[MAXSTRING];        /* snow pack depth and swq (output) */
  char  snowband[MAXSTRING];    /* snow band pack depth and swq (output) */
  char  init_soil[MAXSTRING];	/* soil temp and moist initialization file */
  char  init_snow[MAXSTRING];	/* snowpack initialization file */
} filenames_struct;

typedef struct {
  char   FULL_ENERGY;    /* TRUE = Use full energy code */
  char   FROZEN_SOIL;	 /* TRUE = Use frozen soils code */
  char   SNOW_MODEL;	 /* TRUE = use internal snow model */
  char   CALC_SNOW_FLUX; /* TRUE = compute ground flux interaction with snow */
  char   DIST_PRCP;	 /* TRUE = Use distributed precipitation model */
  char   RADAR;		 /* TRUE = Use radar precipitation fields */
  char   INIT_SOIL;	 /* TRUE = Initial Soil Layers with Named File */
  char   INIT_SNOW;	 /* TRUE = Use Named File to Initialize Snowpack */
  char   FORCE_TYPE[20]; /* type of forcing files being provided */
  char   HP;		 /* TRUE = Use hourly precip in SAWD file */
  char   COMPRESS;       /* TRUE = Compress all output files */
  char   CORRPREC;       /* TRUE = correct precipitation for gage undercatch */
  char   MOISTFRACT;	 /* TRUE = output soil moisture as moisture content */
  char   BINARY_OUTPUT;  /* TRUE = output files are in binary, not ASCII */
  char   ARC_SOIL;       /* TRUE = use ARC/INFO gridded ASCII files for soil 
			  parameters*/
  char   PRT_SNOW_BAND;  /* TRUE = print snow parameters for each snow band */
  int    GRID_DECIMAL;   /* Number of decimal places in grid file extensions */
  int    Nlayer;	 /* Number of layers in model (4 for frozen
			    soils code) */
  int    SNOW_BAND;      /* Number of elevation bands over which to solve the 
			    snow model */
} option_struct;

typedef struct {
  char    debug_dir[512];
  char    DEBUG;
  char    PRT_SOIL;
  char    PRT_VEGE;
  char    PRT_GLOBAL;
  char    PRT_ATMOS;
  char    PRT_SNOW;
  char    PRT_FLUX;
  char    PRT_VAR;
  char    PRT_TEMP;
  char    PRT_MOIST;
  char    PRT_KAPPA;
  char    PRT_BALANCE;
  char    PRT_GRID;
  FILE   *fg_temp;
  FILE   *fg_moist;
  FILE   *fg_kappa;
  FILE   *fg_balance;
  FILE   *fg_energy;
  FILE   *fg_snow;
  FILE   *fg_grid;
  double **inflow[2];
  double **outflow[2];
  double **store_moist[2];
} debug_struct;

/******************************************************************
  This structure records the parameters set by the forcing file
  input routines.  Those filled, are used to estimate the paramters
  needed for the model run in initialize_atmos.c.
  ******************************************************************/
typedef struct {
  char SHORTWAVE;  /* incoming shortwave (W/m^2) */
  char LONGWAVE;   /* incoming longwave (W/m^2) */
  char PRESSURE;   /* atmospheric pressure (kPa) */
  char TSKC;       /* cloud cover (fraction) */
  char SVP;        /* saturated vapor pressure (kPa) */
  char VP;         /* vapor pressure (kPa) */
  char VPD;        /* vapor pressure deficit (kPa) */
  char REL_HUMID;  /* relative humidity (%) */
  char SPEC_HUMID; /* specific humidity (fraction) */
  char ALBEDO;     /* surface albedo (fraction) */
  char AIR_TEMP;   /* air temperature per time step (C) */
  char TMAX;       /* maximum daily temperature (C) */
  char TMIN;       /* minimum daily temperature (C) */
  char PREC;       /* precipitation (mm) */
  char WIND;       /* wind speed (m/s) */
  char DENSITY;    /* atmospheric density (kg/m^3) */
} param_set_struct;

/*******************************************************
  This structure stores all model run global parameters.
  *******************************************************/
typedef struct {
  float  resolution; /* Model resolution (degrees) */
  int	 dt;	     /* Time step in hours (24/dt must be an integer) */
  int    startyear;  /* Starting year of the simulation */
  int    startmonth; /* Starting month of the simulation */
  int    startday;   /* Starting day of the simulation */
  int    starthour;  /* Starting hour of the simulation */
  int    endyear;
  int    endmonth;
  int    endday;
  int    nrecs;      /* Number of time steps simulated */
  int    Nnodes;     /* Number of soil thermal nodes for soil column thermal
		        fluxes calculations */
  double wind_h;     /* height of wind measurements (m) */ 
  double measure_h;  /* height of measurements (m) */
  double MAX_SNOW_TEMP; /* maximum temperature at which snow can fall (C) */
  double MIN_RAIN_TEMP; /* minimum temperature at which rain can fall (C) */
} global_param_struct;

/***********************************************************
  This structure stores the soil parameters for a grid cell.
  ***********************************************************/
typedef struct {
  int     gridcel;	         /* grid cell number */
  float   lat;		         /* grid cell central latitude */
  float   lng;		         /* grid cell central longitude */
  double  b_infilt;  	         /* infiltration parameter */
  double  Ds;		         /* fraction of maximum subsurface flow rate */
  double  Dsmax;  	         /* maximum subsurface flow rate (mm/day) */
  double  Ws;		         /* fraction of maximum soil moisture */
  double  c;                     /* exponent */
  double  expt[MAXlayer];        /* pore-size distribution, HBH 5.15 */
  double  Ksat[MAXlayer];        /* saturated hydraulic  conductivity
				    (mm/day) */
  double  phi_s[MAXlayer];       /* soil moisture diffusion parameter (mm/mm) 
				  */
  double  init_moist[MAXlayer];  /* initial layer moisture level (mm) */
  float   elevation;	         /* grid cell elevation (m) */
  double  depth[MAXlayer];       /* dthickness of each layer (m) */
  double  avg_temp;	         /* average soil temperature (C) */
  double  dp;			 /* soil thermal damping depth (m) */
  double  bubble;	         /* bubbling pressure, HBH 5.15 (cm) */
  double  quartz;		 /* quartz content of soil (fraction) */
  double  resid_moist[MAXlayer]; /* residual moisture content of soil layer */
  double  bulk_density[MAXlayer];/* soil bulk density (kg/m^3) */
  double  soil_density;		 /* soil partical density (kg/m^3) */
  double  rough;		 /* soil surface roughness (m) */
  double  snow_rough;            /* snow surface roughness (m) */
  double  max_moist[MAXlayer];   /* maximum moisture content (mm) per layer */
  double  max_infil;	         /* maximum infiltration rate */
  double  Wcr[MAXlayer];	 /* critical moisture level for soil layer,
			            evaporation is no longer affected moisture
			            stress in the soil (mm) */
  double  Wpwp[MAXlayer];        /* soil moisture content at permanent wilting
			            point (mm) */
  float   time_zone_lng;	 /* central meridian of the time zone */
  double *Tfactor;               /* Change in temperature due to elevation (C) */
  double *Pfactor;               /* Change in Precipitation due to elevation 
				    (fract) */
  double *AreaFract;             /* Fraction of grid cell included in each
				    elevation band */
  double *dz_node;		/* thermal node thickness (m) */
  double *expt_node;		/* soil infiltration parameter */
  double *max_moist_node;	/* maximum soil moisture (mm/mm) */
  double *alpha;		/* thermal solution constant */
  double *beta;			/* thermal solution constant */
  double *gamma;		/* thermal solution constant */
} soil_con_struct;

/*******************************************************************
  This structure stores information about the vegetation coverage of
  the current grid cell.
  *******************************************************************/
typedef struct {
  int    vegetat_type_num;
  int    veg_class;		/* vegetation class reference number */
  double Cv;			/* fraction of vegetation coverage */ 
  double Cv_sum;		/* total fraction of vegetation coverage */
} veg_con_struct;

/******************************************************************
  This structure stores parameters for individual vegetation types.
  ******************************************************************/
typedef struct {
  int    veg_class;		/* vegetation class reference number */
  char   overstory;		/* TRUE = overstory present,
				   important for snow accumulation
				   in canopy */
  double rarc;			/* architectural resistance (s/m) */
  double rmin;			/* minimum stomatal resistance (s/m) */
  double wind_h;		/* height at which wind is measured (m) */
  double root[MAXlayer];	/* percent of roots in each soil layer
                                   (fraction) */
  double LAI[12];		/* monthly leaf area index */
  double displacement[12];	/* vegetation displacement height (m) */
  double roughness[12];		/* vegetation roughness length (m) */
  double Wdmax[12];		/* maximum monthly dew holding capacity (mm) */
  double albedo[12];		/* vegetation albedo (added for full energy)
                                   (fraction) */
  double emissivity[12];	/* vegetation emissivity
                                   (fraction) */
} veg_lib_struct;

/**********************************************************************
  This structure stores the atmospheric forcing data for each time step
  for a single grid cell.
  **********************************************************************/
typedef struct {
  double melt;			/* snow melt (mm) */
  double prec;			/* average precipitation in grid cell (mm) */
  double air_temp;		/* air temperature (C) */
  double rainonly;		/* amount of precip that is rain, not snow
				   (mm) */
  double wind;			/* wind speed (m/s) */
  double rad;			/* net radiation (W/m^2) */
  double vpd;			/* atmospheric vapor pressure deficit (kPa) */
  double vp;			/* atmospheric vapor pressure (kPa) */
  double pressure;		/* atmospheric pressure (kPa) */
  double density;		/* atmospheric density (kg/m^3) */
  double rel_humid;		/* relative humidity (%) */
  double spec_humid;		/* specific humidity (fraction)*/
  double tmin;			/* minimum daily temperature (C) */
  double tmax;			/* maximum daily temperature (C) */
  double priest;		/* preistly evaporation (mm/day) */
  double penman_temp;		/* penman temperature */
  double tskc;			/* total sky cover (fraction) */
  double trans;			/* atmospheric transmissivity */
  double shortwave;		/* incoming shortwave radiation (W/m^2) */
  double net_short;		/* net shortwave at surface (W/m^2) */
  double longwave;		/* incoming longwave radiation (W/m^2) */
  double albedo;		/* bare soil albedo (fraction) */
  double gamma;			/* pshychometric constant (used for snow
                                   met calcs) */
} atmos_data_struct;

/*************************************************************************
  This structure stores information about the time and date of the current
  time step.
  *************************************************************************/
typedef struct {
  int hour;                     /* current hour */
  int day;                      /* current day */
  int month;                    /* current month */
  int year;                     /* current year */
  int day_in_year;		/* julian day in year */
  int day_count;		/* total number of days */
} dmy_struct;			/* array of length nrec created */

/***************************************************************
  This structure stores all soil variables for each layer in the
  soil column.
  ***************************************************************/
typedef struct {
  double moist_thaw;        /* moisture content of the thawed sublayer (mm) */
  double moist_froz;        /* unfrozen moisture content of the frozen
			       sublayer (mm) */
  double moist;             /* moisture content of the unfrozen sublayer
			       (mm) */
  double ice;               /* ice content of the frozen sublayer (mm) */
  double unfrozen;          /* maximum unfrozen water content in frozen
			       layer (mm) */
  double T_thaw;            /* temperature of the thawed sublayer (C) */
  double T_froz;            /* temperature of the frozen sublayer (C) */
  double T;                 /* temperature of the unfrozen sublayer (C) */
  double kappa;             /* average thermal conductivity of the current
			       layer (W/m/K) */
  double Cs;                /* average volumetric heat capacity of the
			       current layer (J/m^3/K) */
  double evap;              /* evapotranspiration from soil layer (mm) */
  double tdepth;            /* thawing front depth in layer (m) */
  double fdepth;            /* freezing front depth in layer (m) */
  double phi;               /* moisture diffusion parameter */
} layer_data_struct;

/******************************************************************
  This structure stores soil variables for the complete soil column 
  for each grid cell.
  ******************************************************************/
typedef struct {
  double aero_resist[3];    /* aerodynamic resistane (s/m)
			       [0] = over bare vegetation or soil
			       [2] = over snow */
  double runoff;            /* runoff from current cell (mm/TS) */
  double baseflow;          /* baseflow from current cell (mm/TS) */
  double inflow;            /* moisture that reaches the top of the soil
			       column (mm) */
  layer_data_struct *layer; /* structure containing soil variables for
			       each layer (see above) */
} cell_data_struct;

/***********************************************************************
  This structure stores energy balance components, and variables used to
  solve the thermal fluxes through the soil column.
  ***********************************************************************/
typedef struct {
  double  shortwave;	        /* incoming shortwave heat (Wm-2) */
  double  longwave;	        /* net longwave flux (Wm-2) */
  double  latent;	        /* net latent heat flux (Wm-2) */
  double  sensible;	        /* net sensible heat flux (Wm-2) */
  double  grnd_flux;	        /* ground heat flux (Wm-2) */
  double  advection;	        /* advective flux (Wm-2) */
  double  deltaH;	        /* change in soil heat storage (Wm-2) */
  double  deltaCC;	        /* change in snow heat storage (Wm-2) */
  double  snow_flux;            /* thermal flux through the snow pack (Wm-2) */
  double  refreeze_energy;      /* energy used to refreeze the snow pack 
				   (Wm-2) */
  double  albedo;	        /* surface albedo (fraction) */
  double  error;	        /* energy balance error (W/m^2) */
  double  Trad[2];              /* surface temperature of energy balance (C) */
  char    frozen;	       	/* TRUE = frozen soil present */
  double *T;			/* thermal node temperatures (C) */
  double *dz;			/* thermal node thickness (m) */
  int     T1_index;		/* soil node at the bottom of the top layer */
  double  kappa[2];		/* soil thermal conductivity for top two
				   layers (W/m/K) */
  double  Cs[2];		/* heat capacity for top two layers
				   (J/m^3/K) */
  double  fdepth[2];		/* [0] freezing front depth, [1] thawing */
  double  unfrozen;		/* frozen layer water content that is 
                                   unfrozen */
  double *ice;			/* frozen layer ice content */
} energy_bal_struct;

/***********************************************************************
  This structure stores vegetation variables for each vegetation type in 
  a grid cell.
  ***********************************************************************/
typedef struct {
  double canopyevap;		/* evaporation from canopy (mm/TS) */
  double Wdew;			/* dew trapped on vegetation (mm) */
  double throughfall;		/* water that reaches the ground through 
                                   the canopy (mm/TS) */
} veg_var_struct;

/************************************************************************
  This structure stores snow pack variables needed to run the snow model.
  ************************************************************************/
typedef struct {
  char   snow;		    /* TRUE = snow, FALSE = no snow */
  int    last_snow;	    /* time steps since last snowfall */
  double snow_canopy;	    /* amount of snow on canopy (m) */
  double swq;               /* snow water equivalent of the entire pack (m) */
  double surf_water;        /* liquid water content of the surface layer (m) */
  double pack_water;        /* liquid water content of the snow pack (m) */
  double surf_temp;         /* depth averaged temperature of the snow pack
                               surface layer (C) */
  double pack_temp;         /* depth averaged temperature of the snow pack
                               (C) */
  double vapor_flux;        /* depth of water evaporation, sublimation, or 
                               condensation from snow pack (m) */
  double canopy_vapor_flux; /* depth of water evaporation, sublimation, or 
                               condensation from intercepted snow (m) */
  double albedo;            /* snow surface albedo (fraction) */
  double coldcontent;       /* cold content of snow pack */
  double mass_error;	    /* snow mass balance error */
  double density;	    /* snow density (kg/m^3) */
  double depth;		    /* snow depth (m) */
  double tmp_int_storage;   /* temporary canopy storage, used in snow_canopy */
  double Qnet;
  double coverage;          /* fraction of snow band that is covered with 
			       snow */
} snow_data_struct;	    /* an array of size Nrec */

/*****************************************************************
  This structure stores all variables needed to solve, or save 
  solututions for all versions of this model.  Vegetation and soil
  variables are created for both wet and dry fractions of the grid
  cell (for use with the distributed precipitation model).
*****************************************************************/
typedef struct {
  double *mu;		          /* fraction of grid cell that receives 
				     precipitation */
  energy_bal_struct **energy;     /* Stores energy balance variables */
  snow_data_struct  **snow;       /* Stores snow variables */
  cell_data_struct  **cell[2];    /* Stores soil layer variables 
				     (wet and dry ) */
  veg_var_struct    **veg_var[2]; /* Stores vegetation variables 
				     (wet and dry) */
} dist_prcp_struct;

/*******************************************************
  This structure stores all variables needed for output.
  *******************************************************/
typedef struct {
  double  evap;		   /* grid cell evaporation */
  double  runoff;	   /* runoff from the surface */
  double  baseflow;	   /* baseflow out of the bottom layer */
  double  moist[MAXlayer]; /* current moisture in each layer */
  double  ice[MAXlayer];   /* frozen layer ice content */
  double  inflow;          /* moisture that reaches the top of the soil
			     column */
  double  fdepth[2];	   /* depth of freezing [0], and thawing[1] fronts */
  double  shortwave;	   /* grid cell incoming shortwave flux */
  double  longwave;	   /* grid cell net longwave flux */
  double  latent;	   /* grid cell net latent heat flux */
  double  sensible;	   /* grid cell net sensible heat flux */
  double  grnd_flux;	   /* grid cell ground flux */
  double  deltaH;	   /* grid cell change in heat storage (snow only) */
  double  energy_error;    /* energy balance error */
  double  surf_temp;	   /* grid cell average daily surface temperature */
  double  r_net;           /* grid cell net radiation W/m^2  */
  double  net_short;       /* grid cell net shortwave flux */
  double  net_long;        /* grid cell net longwave flux */
  double  latent_canop;    /* grid cell net latent heat flux from canopy */
  double  latent_trans;    /* grid cell net latent heat flux from
			      transpiration */
  double  latent_bare;     /* grid cell net latent heat flux from bare soil*/
  double  latent_pet;      /* grid cell potential latent heat flux*/
  double  latent_pet_mm;   /* grid cell potential latent heat flux [mm]*/
  double  rad_temp;        /* grid cell average radiative surface 
			      temperature */
  double  aero_resist;     /* grid cell mean aerodynamic resistence  [s/m] */
  double  surf_cond;       /* grid cell mean surface conductance  [m/s] */
  double  albedo;          /* grid cell mean albedo */ 
  double  prec;            /* incoming precipitation */
  double  Wdew;            /* canopy interception of moisture */
  double *swq;		   /* snow water equivalent (mm) */
  double *snow_canopy;	   /* snow captured by canopy (mm) */
  double *snow_depth;      /* snow depth (cm) */
  double *advection;	   /* grid cell advection (snow only) (Wm-2) */
  double *deltaCC;         /* change of cold content in the snow pack [Wm-2] */
  double *snow_flux;       /* energy flux through the snow pack [Wm-2] */
  double *refreeze_energy; /* energy used to refreeze snowpack [Wm-2] */
  double *coverage;        /* fractional coverage of grid cell with snow */
} out_data_struct;

/********************************************************
  This structure holds all variables needed for the error
  handling routines.
  ********************************************************/
typedef struct {
  int rec;
  double dt;
  out_data_struct *out_data;
  outfiles_struct outfp;
  infiles_struct infp;
  soil_con_struct soil_con;
  veg_con_struct *veg_con;
  veg_var_struct *veg_var;
  energy_bal_struct *energy;
  atmos_data_struct *atmos;
  snow_data_struct *snow;
} Error_struct;

