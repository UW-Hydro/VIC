// $Id$
/**********************************************************************
  Modifications:
  2005-Mar-24 Added data structures to accomodate ALMA variables.	TJB
  2005-Apr-23 Changed ARNO_PARAMS to NIJSSEN2001_BASEFLOW.		TJB
  2005-Apr-23 Added out_data.aero_cond.					TJB
  2005-May-01 Added the ALMA vars CRainf, CSnowf, LSRainf, and LSSnowf.	TJB
  2005-May-02 Added the ALMA vars Wind_E, Wind_N.			TJB
  2005-Dec-21 Removed Trad.                                             GCT
  2006-Sep-23 Implemented flexible output configuration.
              Added output variable types.  Added binary output format
              types.  Removed all output files except the state file from
              the outfiles_struct and the filenames_struct.  Added
              Noutfiles to the option_struct.  Created new out_data_struct
              and out_data_files_struct.  Added new save_data structure.
              Organized the physical constants into one section; got rid
              of redundant Stefan-Boltzmann constant.  Implemented
              aggregation of output variables; added AGG_TYPE definitions.  TJB
  2006-Oct-10 Shortened the names of output variables whose names were
	      too long; fixed typos in others; created new OUT_IN_LONG
	      variable.							TJB
  2006-Oct-16 Merged infiles and outfiles structs into filep_struct;
	      This included merging global->statename to filenames->statefile. TJB
  2006-Nov-07 Added OUT_SOIL_TNODE.					TJB
  2006-Nov-07 Removed LAKE_MODEL option.				TJB
  2006-Nov-07 Organized model constants a bit more.			TJB
  2006-Dec-20 All atmos_data arrays are always dynamically allocated now.	TJB
  2006-Dec-29 Added REL_HUMID to list of supported met input variables.	TJB
  2007-Jan-02 Added ALMA_INPUT option; removed TAIR and PSURF from list
	      of supported met input variables.				TJB
  2007-Jan-15 Added PRT_HEADER option.					TJB
  2007-Apr-03 Added CONTINUEONERROR option.				GCT
  2007-Apr-03 Added ERROR value						KAC
  2007-Apr-24 Added IMPLICIT option.					JCA
  2007-Apr-24 Added EXP_TRANS option.					JCA
  2007-Apr-24 Added Zsum_node to soil_con structure.			JCA
  2007-Aug-08 Added features for EXCESS_ICE option.			JCA
  2007-Aug-22 Added OUTPUT_WATER_ERROR as output variable.		JCA
  2007-Sep-19 Added MAX_SUBSIDENCE parameter to EXCESS_ICE option.	JCA
  2007-Oct-24 Added surf_water to lake_var structure.			KAC via TJB
  2007-Nov-06 Updated lake_var structure with new variables.		LCB via TJB
  2008-Apr-21 Added snow surf_temp, pack_temp, and coldcontent to lake_var
	      structure.						LCB via TJB
  2008-Apr-21 Added SNOW_ALBEDO option.					KAC via TJB
  2008-Apr-21 Added SNOW_DENSITY option.				TJB
  2008-Sep-09 Added SOIL_TNODE_WL as an output variable, the soil
	      temperature in the wetland fraction of the grid cell.	LCB via TJB
*********************************************************************/

#include <user_def.h>
#include <snow.h>

/***** Model Constants *****/
#define MAXSTRING    2048
#define MINSTRING    20
#define HUGE_RESIST  1.e20	/* largest allowable double number */
#define SPVAL        1.e20	/* largest allowable double number - used to signify missing data */
#define SMALL        1.e-12	/* smallest allowable double number */
#define MISSING      -99999.	/* missing value for multipliers in BINARY format */
#define LITTLE 1		/* little-endian flag */
#define BIG 2			/* big-endian flag */
#define ERROR -999              /* Error Flag returned by subroutines */

/***** Met file formats *****/
#define ASCII 1
#define BINARY 2

/***** Snow Albedo parametrizations *****/
#define USACE   0
#define SUN1999 1

/***** Snow Density parametrizations *****/
#define DENS_BRAS   0
#define DENS_SNTHRM 1

/***** Baseflow parametrizations *****/
#define ARNO        0
#define NIJSSEN2001 1

/***** Time Constants *****/
#define DAYS_PER_YEAR 365.
#define HOURSPERDAY   24        /* number of hours per day */
#define HOURSPERYEAR  24*365    /* number of hours per year */
#define SECPHOUR     3600	/* seconds per hour */
#define SEC_PER_DAY 86400.	/* seconds per day */

/***** Physical Constants *****/
#define BARE_SOIL_ALBEDO 0.2	    /* albedo for bare soil */
#define RESID_MOIST      0.0        /* define residual moisture content 
				       of soil column */
#define MAX_ICE_INIT      0.95        /* define maximum volumetric ice fraction
				       of soil column, for EXCESS_ICE option */
#define ICE_AT_SUBSIDENCE 0.8        /* minimum ice/porosity fraction before
					subsidence occurs, for EXCESS_ICE option */
#define MAX_SUBSIDENCE    1.0        /* maximum depth of subsidence per layer per
					time-step (mm) */
#define ice_density      917.	    /* density of ice (kg/m^3) */
#define T_lapse          6.5        /* temperature lapse rate of US Std 
				       Atmos in C/km */
#define von_K        0.40	/* Von Karman constant for evapotranspiration */
#define KELVIN       273.15	/* conversion factor C to K */
#define STEFAN_B     5.6696e-8	/* stefan-boltzmann const in unit W/m^2/K^4 */
#define Lf           3.337e5	/* Latent heat of freezing (J/kg) at 0C */
#define RHO_W        999.842594	/* Density of water (kg/m^3) at 0C */
#define Cp           1010.0	/* Specific heat at constant pressure of air 
				   (J/deg/K) (H.B.H. p.4.7)*/
#define CH_ICE       2100.0e3	/* Volumetric heat capacity (J/(m3*C)) of ice */
#define CH_WATER     4186.8e3   /* volumetric heat capacity of water */
#define K_SNOW       2.9302e-6  /* conductivity of snow (W/mK) */
#define SOLAR_CONSTANT 1400.0	/* Solar constant in W/m^2 */
#define EPS          0.62196351 /* Ratio of molecular weights: M_water_vapor/M_dry_air */
#define G            9.81       /* gravity */
#define JOULESPCAL   4.1868     /* Joules per calorie */
#define GRAMSPKG     1000.      /* convert grams to kilograms */
#define kPa2Pa 1000.            /* converts kPa to Pa */
#define DtoR 0.017453293	/* degrees to radians */
#ifndef PI
#define PI 3.1415927
#endif

/* define constants for saturated vapor pressure curve (kPa) */
#define A_SVP 0.61078
#define B_SVP 17.269
#define C_SVP 237.3

/* define constants for penman evaporation */
#define CP_PM 1013		/* specific heat of moist air at constant pressure (J/kg/C)
				   (Handbook of Hydrology) */
#define PS_PM 101300		/* sea level air pressure in Pa */
#define LAPSE_PM -0.006		/* environmental lapse rate in C/m */

/***** Physical Constraints *****/
#define MINSOILDEPTH 0.001	/* minimum layer depth with which model can
					work (m) */
#define STORM_THRES  0.001      /* thresehold at which a new storm is 
				   decalred */
#define SNOW_DT       5.0	/* Used to bracket snow surface temperatures
				   while computing the snow surface energy 
				   balance (C) */
#define SURF_DT       1.0	/* Used to bracket soil surface temperatures 
                                   while computing energy balance (C) */
#define SOIL_DT       0.25      /* Used to bracket soil temperatures while
                                   solving the soil thermal flux (C) */
#define CANOPY_DT    1.0	/* Used to bracket canopy air temperatures 
                                   while computing energy balance (C) */
#define CANOPY_VP    25.0	/* Used to bracket canopy vapor pressures 
                                   while computing moisture balance (Pa) */

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


/***** Forcing Variable Types *****/
#define N_FORCING_TYPES 23
#define AIR_TEMP   0 /* air temperature per time step [C] (ALMA_INPUT: [K]) */
#define ALBEDO     1 /* surface albedo [fraction] */
#define CRAINF     2 /* convective rainfall [mm] (ALMA_INPUT: [mm/s]) */
#define CSNOWF     3 /* convective snowfall [mm] (ALMA_INPUT: [mm/s]) */
#define DENSITY    4 /* atmospheric density [kg/m3] */
#define LONGWAVE   5 /* incoming longwave radiation [W/m2] */
#define LSRAINF    6 /* large-scale rainfall [mm] (ALMA_INPUT: [mm/s]) */
#define LSSNOWF    7 /* large-scale snowfall [mm] (ALMA_INPUT: [mm/s]) */
#define PREC       8 /* total precipitation (rain and snow) [mm] (ALMA_INPUT: [mm/s]) */
#define PRESSURE   9 /* atmospheric pressure [kPa] (ALMA_INPUT: [Pa]) */
#define QAIR      10 /* specific humidity [kg/kg] */
#define RAINF     11 /* rainfall (convective and large-scale) [mm] (ALMA_INPUT: [mm/s]) */
#define REL_HUMID 12 /* relative humidity [fraction] */
#define SHORTWAVE 13 /* incoming shortwave [W/m2] */
#define SNOWF     14 /* snowfall (convective and large-scale) [mm] (ALMA_INPUT: [mm/s]) */
#define TMAX      15 /* maximum daily temperature [C] (ALMA_INPUT: [K]) */
#define TMIN      16 /* minimum daily temperature [C] (ALMA_INPUT: [K]) */
#define TSKC      17 /* cloud cover [fraction] */
#define VP        18 /* vapor pressure [kPa] (ALMA_INPUT: [Pa]) */
#define WIND      19 /* wind speed [m/s] */
#define WIND_E    20 /* zonal component of wind speed [m/s] */
#define WIND_N    21 /* meridional component of wind speed [m/s] */
#define SKIP      22 /* place holder for unused data columns */

/***** Output Variable Types *****/
#define N_OUTVAR_TYPES 120
// Water Balance Terms - state variables
#define OUT_LAKE_DEPTH       0  /* lake depth (distance between surface and deepest point) [m] */
#define OUT_LAKE_ICE         1  /* moisture stored as lake ice [mm over lake ice area] */
#define OUT_LAKE_ICE_FRACT   2  /* fractional coverage of lake ice [fraction] */
#define OUT_LAKE_ICE_HEIGHT  3  /* thickness of lake ice [cm] */
#define OUT_LAKE_MOIST       4  /* liquid water and ice stored in lake [mm over grid cell] */
#define OUT_LAKE_SURF_AREA   5  /* lake surface area [m2] */
#define OUT_LAKE_VOLUME      6  /* lake volume [m3] */
#define OUT_ROOTMOIST        7  /* root zone soil moisture  [mm] */
#define OUT_SMFROZFRAC       8  /* fraction of soil moisture (by mass) that is ice, for each soil layer */
#define OUT_SMLIQFRAC        9  /* fraction of soil moisture (by mass) that is liquid, for each soil layer */
#define OUT_SNOW_CANOPY     10  /* snow interception storage in canopy  [mm] */
#define OUT_SNOW_COVER      11  /* fractional area of snow cover [fraction] */
#define OUT_SNOW_DEPTH      12  /* depth of snow pack [cm] */
#define OUT_SOIL_ICE        13  /* soil ice content  [mm] for each soil layer */
#define OUT_SOIL_LIQ        14  /* soil liquid content  [mm] for each soil layer */
#define OUT_SOIL_MOIST      15  /* soil total moisture content  [mm] for each soil layer */
#define OUT_SOIL_WET        16  /* vertical average of (soil moisture - wilting point)/(maximum soil moisture - wilting point) [mm/mm] */
#define OUT_SURFSTOR        17  /* storage of liquid water and ice (not snow) on surface (ponding) [mm] */
#define OUT_SURF_FROST_FRAC 18  /* fraction of soil surface that is frozen [fraction] */
#define OUT_SWE             19  /* snow water equivalent in snow pack (including vegetation-intercepted snow)  [mm] */
#define OUT_WDEW            20  /* total moisture interception storage in canopy [mm] */
// Water Balance Terms - fluxes
#define OUT_BASEFLOW        21  /* baseflow out of the bottom layer  [mm] (ALMA_OUTPUT: [mm/s]) */
#define OUT_DELINTERCEPT    22  /* change in canopy interception storage  [mm] */
#define OUT_DELSOILMOIST    23  /* change in soil water content  [mm] */
#define OUT_DELSURFSTOR     24  /* change in surface liquid water storage  [mm] */
#define OUT_DELSWE          25  /* change in snow water equivalent  [mm] */
#define OUT_EVAP            26  /* total net evaporation [mm] (ALMA_OUTPUT: [mm/s]) */
#define OUT_EVAP_BARE       27  /* net evaporation from bare soil [mm] (ALMA_OUTPUT: [mm/s]) */
#define OUT_EVAP_CANOP      28  /* net evaporation from canopy interception [mm] (ALMA_OUTPUT: [mm/s]) */
#define OUT_EVAP_LAKE       29  /* net evaporation from lake surface [mm] (ALMA_OUTPUT: [mm/s]) */
#define OUT_INFLOW          30  /* moisture that reaches top of soil column [mm] (ALMA_OUTPUT: [mm/s]) */
#define OUT_PREC            31  /* incoming precipitation [mm] (ALMA_OUTPUT: [mm/s]) */
#define OUT_RAINF           32  /* rainfall  [mm] (ALMA_OUTPUT: [mm/s]) */
#define OUT_REFREEZE        33  /* refreezing of water in the snow  [mm] (ALMA_OUTPUT: [mm/s]) */
#define OUT_RUNOFF          34  /* surface runoff [mm] (ALMA_OUTPUT: [mm/s]) */
#define OUT_SNOW_MELT       35  /* snow melt  [mm] (ALMA_OUTPUT: [mm/s]) */
#define OUT_SNOWF           36  /* snowfall  [mm] (ALMA_OUTPUT: [mm/s]) */
#define OUT_SUB_BLOWING     37  /* net sublimation of blowing snow [mm] (ALMA_OUTPUT: [mm/s]) */
#define OUT_SUB_CANOP       38  /* net sublimation from snow stored in canopy [mm] (ALMA_OUTPUT: [mm/s]) */
#define OUT_SUB_SNOW        39  /* total net sublimation from snow pack (surface and blowing) [mm] (ALMA_OUTPUT: [mm/s]) */
#define OUT_SUB_SURFACE     40  /* net sublimation from snow pack surface [mm] (ALMA_OUTPUT: [mm/s]) */
#define OUT_TRANSP_VEG      41  /* net transpiration from vegetation [mm] (ALMA_OUTPUT: [mm/s]) */
#define OUT_WATER_ERROR     42  /* water budget error [mm] */
// Energy Balance Terms - state variables
#define OUT_ALBEDO          43  /* albedo [fraction] */
#define OUT_BARESOILT       44  /* bare soil surface temperature [C] (ALMA_OUTPUT: [K]) */
#define OUT_FDEPTH          45  /* depth of freezing fronts [cm] (ALMA_OUTPUT: [m]) for each freezing front */
#define OUT_LAKE_ICE_TEMP   46  /* temperature of lake ice [K] */
#define OUT_LAKE_SURF_TEMP  47  /* lake surface temperature [K] */
#define OUT_RAD_TEMP        48  /* average radiative surface temperature [K] */
#define OUT_SALBEDO         49  /* snow albedo [fraction] */
#define OUT_SNOW_PACK_TEMP  50  /* snow pack temperature [C] (ALMA_OUTPUT: [K]) */
#define OUT_SNOW_SURF_TEMP  51  /* snow surface temperature [C] (ALMA_OUTPUT: [K]) */
#define OUT_SOIL_TEMP       52  /* soil temperature [C] (ALMA_OUTPUT: [K]) for each soil layer */
#define OUT_SOIL_TNODE      53  /* soil temperature [C] (ALMA_OUTPUT: [K]) for each soil thermal node */
#define OUT_SOIL_TNODE_WL   54  /* soil temperature [C] (ALMA_OUTPUT: [K]) for each soil thermal node in the wetland */
#define OUT_SURF_TEMP       55  /* average surface temperature [C] (ALMA_OUTPUT: [K]) */
#define OUT_TDEPTH          56  /* depth of thawing fronts [cm] (ALMA_OUTPUT: [m]) for each thawing front */
#define OUT_VEGT            57  /* average vegetation canopy temperature [C] (ALMA_OUTPUT: [K]) */
// Energy Balance Terms - fluxes
#define OUT_ADV_SENS        58 /* net sensible flux advected to snow pack [W/m2] */
#define OUT_ADVECTION       59  /* advected energy [W/m2] */
#define OUT_DELTACC         60  /* rate of change in cold content in snow pack [W/m2] (ALMA_OUTPUT: [J/m2]) */
#define OUT_DELTAH          61  /* rate of change in heat storage [W/m2] (ALMA_OUTPUT: [J/m2]) */
#define OUT_ENERGY_ERROR    62  /* energy budget error [W/m2] */
#define OUT_FUSION          63  /* net energy used to melt/freeze soil moisture [W/m2] */
#define OUT_GRND_FLUX       64  /* net heat flux into ground [W/m2] */
#define OUT_IN_LONG         65  /* incoming longwave at ground surface (under veg) [W/m2] */
#define OUT_LATENT          66  /* net upward latent heat flux [W/m2] */
#define OUT_LATENT_SUB      67  /* net upward latent heat flux from sublimation [W/m2] */
#define OUT_MELT_ENERGY     68  /* energy of fusion (melting) in snowpack [W/m2] */
#define OUT_NET_LONG        69  /* net downward longwave flux [W/m2] */
#define OUT_NET_SHORT       70  /* net downward shortwave flux [W/m2] */
#define OUT_R_NET           71  /* net downward radiation flux [W/m2] */
#define OUT_RFRZ_ENERGY     72  /* net energy used to refreeze liquid water in snowpack [W/m2] */
#define OUT_SENSIBLE        73  /* net upward sensible heat flux [W/m2] */
#define OUT_SNOW_FLUX       74  /* energy flux through snow pack [W/m2] */
// Miscellaneous Terms
#define OUT_AERO_RESIST     75  /* canopy aerodynamic resistance [s/m] */
#define OUT_AERO_COND       76  /* canopy aerodynamic conductance [m/s] */
#define OUT_AIR_TEMP        77  /* air temperature [C] (ALMA_OUTPUT: [K])*/
#define OUT_DENSITY         78  /* near-surface atmospheric density [kg/m3]*/
#define OUT_LONGWAVE        79  /* incoming longwave [W/m2] */
#define OUT_PRESSURE        80  /* near surface atmospheric pressure [kPa] (ALMA_OUTPUT: [Pa])*/
#define OUT_QAIR            81  /* specific humidity [kg/kg] */
#define OUT_REL_HUMID       82  /* relative humidity [fraction]*/
#define OUT_SHORTWAVE       83  /* incoming shortwave [W/m2] */
#define OUT_SURF_COND       84  /* surface conductance [m/s] */
#define OUT_VP              85  /* near surface vapor pressure [kPa] (ALMA_OUTPUT: [Pa]) */
#define OUT_WIND            86  /* near surface wind speed [m/s] */
// Band-specific quantities
#define OUT_ADV_SENS_BAND        87  /* net sensible heat flux advected to snow pack [W/m2] */
#define OUT_ADVECTION_BAND       88  /* advected energy [W/m2] */
#define OUT_ALBEDO_BAND          89  /* albedo [fraction] */
#define OUT_DELTACC_BAND         90  /* change in cold content in snow pack [W/m2] */
#define OUT_GRND_FLUX_BAND       91  /* net heat flux into ground [W/m2] */
#define OUT_IN_LONG_BAND         92  /* incoming longwave at ground surface (under veg) [W/m2] */
#define OUT_LATENT_BAND          93  /* net upward latent heat flux [W/m2] */
#define OUT_LATENT_SUB_BAND      94  /* net upward latent heat flux due to sublimation [W/m2] */
#define OUT_MELT_ENERGY_BAND     95  /* energy of fusion (melting) in snowpack [W/m2] */
#define OUT_NET_LONG_BAND        96  /* net downward longwave flux [W/m2] */
#define OUT_NET_SHORT_BAND       97  /* net downward shortwave flux [W/m2] */
#define OUT_RFRZ_ENERGY_BAND     98  /* net energy used to refreeze liquid water in snowpack [W/m2] */
#define OUT_SENSIBLE_BAND        99  /* net upward sensible heat flux [W/m2] */
#define OUT_SNOW_CANOPY_BAND    100  /* snow interception storage in canopy [mm] */
#define OUT_SNOW_COVER_BAND     101  /* fractional area of snow cover [fraction] */
#define OUT_SNOW_DEPTH_BAND     102  /* depth of snow pack [cm] */
#define OUT_SNOW_FLUX_BAND      103  /* energy flux through snow pack [W/m2] */
#define OUT_SNOW_MELT_BAND      104  /* snow melt [mm] (ALMA_OUTPUT: [mm/s]) */
#define OUT_SNOW_PACKT_BAND     105  /* snow pack temperature [C] (ALMA_OUTPUT: [K]) */
#define OUT_SNOW_SURFT_BAND     106  /* snow surface temperature [C] (ALMA_OUTPUT: [K]) */
#define OUT_SWE_BAND            107  /* snow water equivalent in snow pack [mm] */
// Dynamic Soil Property Terms - EXCESS_ICE option
#if EXCESS_ICE
#define OUT_SOIL_DEPTH          108  /* soil moisture layer depths [m] */
#define OUT_SUBSIDENCE          109  /* subsidence of soil layer [mm] */
#define OUT_POROSITY            110  /* porosity [mm/mm] */
#define OUT_ZSUM_NODE           111  /* depths of thermal nodes [m] */
#endif // EXCESS_ICE

/***** Output BINARY format types *****/
#define OUT_TYPE_DEFAULT 0 /* Default data type */
#define OUT_TYPE_CHAR    1 /* char */
#define OUT_TYPE_SINT    2 /* short int */
#define OUT_TYPE_USINT   3 /* unsigned short int */
#define OUT_TYPE_INT     4 /* int */
#define OUT_TYPE_FLOAT   5 /* single-precision floating point */
#define OUT_TYPE_DOUBLE  6 /* double-precision floating point */

/***** Output aggregation method types *****/
#define AGG_TYPE_AVG     0 /* average over agg interval */
#define AGG_TYPE_BEG     1 /* value at beginning of agg interval */
#define AGG_TYPE_END     2 /* value at end of agg interval */
#define AGG_TYPE_MAX     3 /* maximum value over agg interval */
#define AGG_TYPE_MIN     4 /* minimum value over agg interval */
#define AGG_TYPE_SUM     5 /* sum over agg interval */

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
  FILE *init_state;     /* initial model state file */
  FILE *lakeparam;      /* lake parameter file */
  FILE *snowband;       /* snow elevation band data file */
  FILE *soilparam;      /* soil parameters for all grid cells */
  FILE *statefile;      /* output model state file */
  FILE *veglib;         /* vegetation parameters for all vege types */
  FILE *vegparam;       /* fractional coverage info for grid cell */
} filep_struct;

typedef struct {
  char  forcing[2][MAXSTRING];  /* atmospheric forcing data file names */
  char  f_path_pfx[2][MAXSTRING];  /* path and prefix for atmospheric forcing data file names */
  char  global[MAXSTRING];      /* global control file name */
  char  init_state[MAXSTRING];  /* initial model state file name */
  char  lakeparam[MAXSTRING];   /* lake model constants file */
  char  result_dir[MAXSTRING];  /* directory where results will be written */
  char  snowband[MAXSTRING];    /* snow band parameter file name */
  char  soil[MAXSTRING];        /* soil parameter file name, or name of 
				   file that has a list of all aoil 
				   ARC/INFO files */
  char  soil_dir[MAXSTRING];    /* directory from which to read ARC/INFO 
				   soil files */
  char  statefile[MAXSTRING];   /* name of file in which to store model state */
  char  veg[MAXSTRING];         /* vegetation grid coverage file */
  char  veglib[MAXSTRING];      /* vegetation parameter library file */
} filenames_struct;

typedef struct {

  // simulation modes
  char   BLOWING;        /* TRUE = calculate sublimation from blowing snow */
  char   CORRPREC;       /* TRUE = correct precipitation for gage undercatch */
  char   DIST_PRCP;      /* TRUE = Use distributed precipitation model */
  char   EQUAL_AREA;     /* TRUE = RESOLUTION stores grid cell area in km^2;
			    FALSE = RESOLUTION stores grid cell side length in degrees */
  char   FROZEN_SOIL;    /* TRUE = Use frozen soils code */
  char   FULL_ENERGY;    /* TRUE = Use full energy code */
  char   GRND_FLUX;      /* TRUE = compute ground heat flux and energy 
			    balance */
  char   LAKES;          /* TRUE = use lake energy code */
  float  MIN_WIND_SPEED; /* Minimum wind speed in m/s that can be used by 
			    the model. **/
  char   MOISTFRACT;     /* TRUE = output soil moisture as moisture content */
  int    Nlakenode;      /* Number of lake thermal nodes in the model. */
  int    Nlayer;         /* Number of layers in model */
  int    Nnode;          /* Number of soil thermal nodes in the model */
  char   NOFLUX;         /* TRUE = Use no flux lower bondary when computing 
			    soil thermal fluxes */
  char   IMPLICIT;         /* TRUE = Use implicit solution when computing 
			    soil thermal fluxes */
  char   EXP_TRANS;         /* TRUE = Uses grid transform for exponential node 
			       distribution for soil heat flux calculations*/
  float  PREC_EXPT;      /* Exponential that controls the fraction of a
			    grid cell that receives rain during a storm
			    of given intensity */
  int    ROOT_ZONES;     /* Number of root zones used in simulation */
  char   QUICK_FLUX;     /* TRUE = Use Liang et al., 1999 formulation for
			    ground heat flux, if FALSE use explicit finite
			    difference method */
  char   QUICK_SOLVE;    /* TRUE = Use Liang et al., 1999 formulation for 
			    iteration, but explicit finite difference
			    method for final step. */
  char   SNOW_ALBEDO;    /* USACE: Use algorithm of US Army Corps of Engineers, 1956; SUN1999: Use algorithm of Sun et al., JGR, 1999 */
  char   SNOW_DENSITY;   /* DENS_BRAS: Use algorithm of Bras, 1990; DENS_SNTHRM: Use algorithm of SNTHRM89 adapted for 1-layer pack */
  int    SNOW_BAND;      /* Number of elevation bands over which to solve the 
			    snow model */
  int    SNOW_STEP;      /* Time step in hours to use when solving the 
			    snow model */

  // input options
  char   ALMA_INPUT;     /* TRUE = input variables are in ALMA-compliant units; FALSE = standard VIC units */
  char   ARC_SOIL;       /* TRUE = use ARC/INFO gridded ASCII files for soil 
			    parameters*/
  char   BASEFLOW;       /* ARNO: read Ds, Dm, Ws, c; NIJSSEN2001: read d1, d2, d3, d4 */
  int    GRID_DECIMAL;   /* Number of decimal places in grid file extensions */
  char   GLOBAL_LAI;     /* TRUE = read LAI values for each vegetation type
			    from the veg param file */
  char   LAKE_PROFILE;   /* TRUE = user-specified lake/area profile */

  // state options
  char   BINARY_STATE_FILE; /* TRUE = model state file is binary (default) */
  char   INIT_STATE;     /* TRUE = initialize model state from file */
  char   SAVE_STATE;     /* TRUE = save state file */       

  // output options
  char   ALMA_OUTPUT;    /* TRUE = output variables are in ALMA-compliant units; FALSE = standard VIC units */
  char   BINARY_OUTPUT;  /* TRUE = output files are in binary, not ASCII */
  char   COMPRESS;       /* TRUE = Compress all output files */
  int    CONTINUEONERROR;/* TRUE = VIC will continue to run after a cell has an error */
  int    Noutfiles;      /* Number of output files (not including state files) */
  char   PRT_HEADER;     /* TRUE = insert header at beginning of output file; FALSE = no header */
  char   PRT_SNOW_BAND;  /* TRUE = print snow parameters for each snow band */

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
  double MAX_SNOW_TEMP; /* maximum temperature at which snow can fall (C) */
  double MIN_RAIN_TEMP; /* minimum temperature at which rain can fall (C) */
  double measure_h;  /* height of measurements (m) */
  double wind_h;     /* height of wind measurements (m) */ 
  float  resolution; /* Model resolution (degrees) */
  int    dt;         /* Time step in hours (24/dt must be an integer) */
  int    out_dt;     /* Output time step in hours (24/out_dt must be an integer) */
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

/******************************************************************
  This structure stores the lake/wetland parameters for a grid cell
  ******************************************************************/
typedef struct {
  // General information
  int    wetland_veg_class;       /* Vegetation class of the wetland */
  // Lake basin dimensions
  int    numnod;                  /* Maximum number of lake nodes for this grid cell */
  double z[MAX_LAKE_NODES+1];     /* Elevation of each lake node (when lake storage is at maximum), relative to lake's deepest point (m) */  
  double basin[MAX_LAKE_NODES+1]; /* Area of lake basin at each lake node (when lake storage is at maximum) (m^2) */
  double Cl[MAX_LAKE_NODES+1];    /* Fractional coverage of lake basin at each node (when lake storage is at maximum) (fraction of grid cell area) */
  double b;                       /* Exponent in default lake depth-area profile (y=Ax^b) */
  double maxdepth;                /* Maximum allowable lake depth (m) */
  double mindepth;                /* Minimum allowable lake depth (m) */
  double maxvolume;               /* Lake volume when lake depth is at maximum (m^3) */
  double minvolume;               /* Lake volume when lake depth is at minimum (m^3) */
  // Hydrological properties
  float  bpercent;                /* Fraction of wetland baseflow (subsurface runoff) that flows into lake */
  float  rpercent;                /* Fraction of wetland surface runoff that flows into lake */
  double eta_a;                   /* Decline of solar radiation w/ depth (m^-1) */ /* not currently used */
  double wfrac;                   /* Width of lake outlet, expressed as fraction of lake perimeter */
  // Initial conditions
  double depth_in;                /* Initial lake depth (distance from surface to deepest point) (m) */
} lake_con_struct;

/*****************************************************************
  This structure stores the lake/wetland variables for a grid cell
  *****************************************************************/
typedef struct {
  // Current lake dimensions and liquid water state variables
  int    activenod;               /* Number of nodes whose corresponding layers currently contain water */
  double dz;                      /* Vertical thickness of all horizontal water layers below the surface layer (m) */
  double surfdz;                  /* Vertical thickness of surface (top) water layer (m) */
  double ldepth;                  /* Current depth of liquid water in lake (distance from surface to deepest point) (m) */
  double surface[MAX_LAKE_NODES+1];/* Area of horizontal cross-section of lake at each node (at end of time step) (m^2) */
  double sarea;                   /* Current surface area of liquid water in lake (at beginning of time step) (m^2) */
  double volume;                  /* Current lake water volume, including liquid water equivalent of lake ice and snow (m^3) */
  double temp[MAX_LAKE_NODES];    /* Lake water temperature at each node (C) */
  double tempavg;                 /* Average liquid water temperature of entire lake (C) */
  // Current properties (state variables) specific to lake ice/snow
  double areai;                   /* Area of ice coverage (at beginning of time step) (m^2) */
  double new_ice_area;            /* Area of ice coverage (at end of time step) (m^2) */
  double ice_water_eq;            /* Liquid water equivalent volume of lake ice (m^3) */
  double hice;                    /* Height of lake ice at thickest point (m) */ 
  double tempi;                   /* Lake ice temperature (C) */
  double swe;                     /* Water equivalence of lake snow cover (m over lake ice area) */
  double surf_temp;               /* Temperature of surface snow layer (C) */
  double pack_temp;               /* Temperature of pack snow layer (C) */
  double coldcontent;             /* cold content of snow pack */
  double surf_water;              /* Water content of surface snow layer (m over lake ice area) */
  double pack_water;              /* Water content of pack snow layer (m over lake ice area) */
  double SAlbedo;                 /* Albedo of lake snow (fraction) */
  double sdepth;                  /* Depth of snow on top of ice (m over lake ice area) */
  // Other current lake properties (derived from state variables and forcings)
  double aero_resist;		  /* Aerodynamic resistance (s/m) */
  double aero_resist_used;	  /* Aerodynamic resistance (s/m) after stability correction */
  double density[MAX_LAKE_NODES]; /* Lake water density profile (kg/m^3) */
  // Moisture fluxes
  double evapw;                   /* Evaporative flux from lake (and ice/snow) surface (mm over lake area) */
  double baseflow_in;             /* Baseflow into lake from the rest of the grid cell (mm over grid cell area) */
  double baseflow_out;            /* Baseflow out of lake to channel network (mm over grid cell area) */
  double runoff_in;               /* Surface runoff into lake from the rest of the grid cell (mm over grid cell area) */
  double runoff_out;              /* Surface runoff out of lake to channel network (mm over grid cell area) */
  double snowmlt;                 /* Moisture released by melting of lake snow (mm over lake ice area) */
} lake_var_struct;

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
#if EXCESS_ICE
  double   Ds_orig;                   /* fraction of maximum subsurface flow 
					 rate */
  double   Dsmax_orig;                /* maximum subsurface flow rate 
					 (mm/day) */
  double   Ws_orig;                   /* fraction of maximum soil moisture */
#endif  
  double   alpha[MAX_NODES];          /* thermal solution constant */
  double   annual_prec;               /* annual average precipitation (mm) */
  double   avg_temp;                  /* average soil temperature (C) */
  double   b_infilt;                  /* infiltration parameter */
  double   beta[MAX_NODES];           /* thermal solution constant */
  double   bubble[MAX_LAYERS];        /* bubbling pressure, HBH 5.15 (cm) */
  double   bubble_node[MAX_NODES];    /* bubbling pressure (cm) */
  double   bulk_density[MAX_LAYERS];  /* soil bulk density (kg/m^3) */
  double   c;                         /* exponent in ARNO baseflow scheme */
  double   depth[MAX_LAYERS];         /* thickness of each soil moisture 
					 layer (m).  In the case of EXCESS_ICE,
				         this is the effective (dynamic) depth. */
#if SPATIAL_SNOW
  double   depth_full_snow_cover;     // minimum depth for full snow cover
#endif // SPATIAL_SNOW
  double   dp;                        /* soil thermal damping depth (m) */
  double   dz_node[MAX_NODES];        /* thermal node thickness (m) */
  double   Zsum_node[MAX_NODES];      /* thermal node depth (m) */
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
  char    *AboveTreeLine;             // Flag to indicate if band is above 
                                      // the treeline
#if QUICK_FS
  double **ufwc_table_layer[MAX_LAYERS];
  double **ufwc_table_node[MAX_NODES]; 
#endif
  float    elevation;                 /* grid cell elevation (m) */
  float    lat;                       /* grid cell central latitude */
  float    lng;                       /* grid cell central longitude */
  double   cell_area;                 /* Area of grid cell (m^2) */
  float    time_zone_lng;             /* central meridian of the time zone */
  float  **layer_node_fract;          /* fraction of all nodes within each 
					 layer */
  int      gridcel;                   /* grid cell number */
#if EXCESS_ICE
  double   min_depth[MAX_LAYERS];     /* soil layer depth as given in the soil file (m). 
					 The effective depth will always be >= this value. */
  double   porosity_node[MAX_NODES];  /* porosity for each thermal node */
  double   effective_porosity[MAX_LAYERS]; /* effective soil porosity (fraction)
					   when soil pores are expanded due to
					   excess ground ice */
  double   effective_porosity_node[MAX_NODES]; /* effective soil porosity (fraction)
					   when soil pores are expanded due to
					   excess ground ice */
  double   Wcr_FRACT[MAX_LAYERS];
  double   Wpwp_FRACT[MAX_LAYERS];
  double   subsidence[MAX_LAYERS];      /* subsidence of soil layer, mm*/
#endif // EXCESS_ICE
} soil_con_struct;

/*****************************************************************
  This structure stores the dynamic soil properties for a grid cell
  *****************************************************************/
#if EXCESS_ICE
typedef struct {
  double soil_depth[MAX_LAYERS];             /* soil moisture layer depths [m] */
  double subsidence[MAX_LAYERS];             /* subsidence of soil layer [mm] */
  double porosity[MAX_LAYERS];               /* porosity [mm/mm] */
  double zsum_node[MAX_NODES];               /* depths of thermal nodes [m] */
} dynamic_soil_struct;
#endif // EXCESS_ICE

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
  char   *snowflag;  /* TRUE if there is snowfall in any of the snow
                        bands during the timestep, FALSE otherwise*/
  double *air_temp;  /* air temperature (C) */
  double *density;   /* atmospheric density (kg/m^3) */
  double *longwave;  /* incoming longwave radiation (W/m^2) (net incoming
                        longwave for water balance model) */
  double out_prec;   /* Total precipitation for time step - accounts
                        for corrected precipitation totals */
  double out_rain;   /* Rainfall for time step (mm) */
  double out_snow;   /* Snowfall for time step (mm) */
  double *prec;      /* average precipitation in grid cell (mm) */
  double *pressure;  /* atmospheric pressure (kPa) */
  double *shortwave; /* incoming shortwave radiation (W/m^2) */
  double *vp;        /* atmospheric vapor pressure (kPa) */
  double *vpd;       /* atmospheric vapor pressure deficit (kPa) */
  double *wind;      /* wind speed (m/s) */
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
  double rootmoist;                    /* total of layer.moist over all layers
                                          in the root zone (mm) */
  double wetness;                      /* average of
                                          (layer.moist - Wpwp)/(porosity*depth - Wpwp)
                                          over all layers (fraction) */
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
  int    snow;              /* TRUE = snow, FALSE = no snow */
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
  double transport;	    /* flux of snow (potentially) transported from veg type */
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
  lake_var_struct     lake_var;   /* Stores lake/wetland variables */
  snow_data_struct  **snow;       /* Stores snow variables */
  veg_var_struct    **veg_var[2]; /* Stores vegetation variables (wet and 
				     dry) */
} dist_prcp_struct;

/*******************************************************
  This structure stores moisture state information for
  differencing with next time step.
  *******************************************************/
typedef struct {
  double	total_soil_moist; /* total column soil moisture [mm] */
  double	surfstor;         /* surface water storage [mm] */
  double	swe;              /* snow water equivalent [mm] */
  double	wdew;             /* canopy interception [mm] */
} save_data_struct;

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
  int		aggtype;     /* type of aggregation to use;
				AGG_TYPE_AVG    = take average value over agg interval
				AGG_TYPE_BEG    = take value at beginning of agg interval
				AGG_TYPE_END    = take value at end of agg interval
				AGG_TYPE_MAX    = take maximum value over agg interval
				AGG_TYPE_MIN    = take minimum value over agg interval
				AGG_TYPE_SUM    = take sum over agg interval */
  int		nelem;       /* number of data values */
  double	*data;       /* array of data values */
  double	*aggdata;    /* array of aggregated data values */
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
  filep_struct       filep;
  int                rec;
  out_data_struct   *out_data;
  out_data_file_struct    *out_data_files;
  snow_data_struct  *snow;
  soil_con_struct    soil_con;
  veg_con_struct    *veg_con;
  veg_var_struct    *veg_var;
} Error_struct;

