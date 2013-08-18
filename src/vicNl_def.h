// $Id$
/**********************************************************************
  This file contains "#define" statements and "typedef" statements.
  It also contains "extern" declarations for global variables.  For such
  variables, a single declaration/definition of the global variable, not
  containing the word "extern", must exist in global.h.  This is because
  global.h is only included by one file (vicNl.c), while vicNl_def.h is
  included multiple times (all *.c files) via vicNl.h.

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
  2009-Jan-12 Added COMPUTE_TREELINE and JULY_TAVG_SUPPLIED options.	TJB
  2009-Jan-16 Modified aero_resist_used and Ra_used to become arrays of
	      two elements (surface and overstory); added
	      options.AERO_RESIST_CANSNOW.				TJB
  2009-Jan-16 Added AERO_COND1&2 and AERO_RESIST1&2 to track
	      surface and overstory values; changed AERO_COND
	      and AERO_RESIST to track "scene" values.			TJB
  2009-Feb-09 Updated description of PRT_SNOW_BAND option.		TJB
  2009-Feb-22 Added OUT_VPD.						TJB
  2009-Mar-16 Added min_liq to the layer_data_struct.			TJB
  2009-May-17 Added OUT_ASAT.						TJB
  2009-May-17 Added AR_406_LS to options.AERO_RESIST_CANSNOW.		TJB
  2009-May-17 Added options.MIN_LIQ.					TJB
  2009-May-18 Added options.PLAPSE and Rd, the gas constant for dry
	      air.							TJB
  2009-May-20 Added options.GRND_FLUX_TYPE.				TJB
  2009-May-22 Added TFALLBACK value to options.CONTINUEONERROR.		TJB
  2009-Jun-09 Modified to use extension of veg_lib structure to contain
	      bare soil information.					TJB
  2009-Jun-09 Added OUT_PET_*, potential evap for various reference
	      land cover types.						TJB
  2009-Jun-09 Cell_data structure now only stores final aero_resist
	      values (called "aero_resist").  Preliminary uncorrected
	      aerodynamic resistances for current vegetation and various
	      reference land cover types for use in potential evap
	      calculations is stored in temporary array aero_resist.	TJB
  2009-Jun-19 Added T flag to indicate whether TFALLBACK occurred.	TJB
  2009-Jul-07 Added BandElev[] to soil_con_struct.			TJB
  2009-Jul-31 Added lake_idx to lake_con struct and LAKE to veg_con
	      struct.							TJB
  2009-Aug-28 OUT_LAKE_ICE_TEMP and OUT_LAKE_SURF_TEMP are [C].		TJB
  2009-Sep-19 Added T fbcount to count TFALLBACK occurrences.		TJB
  2009-Sep-19 Changed Cp to be 1013, the value for moist air.		TJB
  2009-Sep-19 Made TFALLBACK a separate option from CONTINUEONERROR.	TJB
  2009-Sep-28 Added snow and energy structures to lake_var_struct.	TJB
  2009-Sep-30 Miscellaneous fixes for lake model.			TJB
  2009-Oct-08 Extended T fallback scheme to snow and ice T.		TJB
  2009-Nov-09 Changed definition of sarea to include ice extent.	LCB via TJB
  2009-Dec-11 Removed min_liq and options.MIN_LIQ.			TJB
  2010-Feb-14 Added OUT_LAKE_AREA_FRAC.					TJB
  2010-Mar-31 Added RUNOFF_IN and OUT_RUNOFF_IN.			TJB
  2010-Apr-28 Replaced GLOBAL_LAI with VEGPARAM_LAI and LAI_SRC.	TJB
  2010-Sep-24 Renamed RUNOFF_IN and OUT_RUNOFF_IN to CHANNEL_IN and
	      OUT_LAKE_CHAN_IN, respectively.  Renamed OUT_EVAP_LAKE
	      to OUT_LAKE_EVAP.  Added other lake water balance terms
	      to set of output variables.  Added volumetric versions 
	      of these too.						TJB
  2010-Nov-02 Added OUT_LAKE_RO_IN, OUT_LAKE_RO_IN_V, OUT_LAKE_RCHRG,
	      OUT_LAKE_RCHRG_V, OUT_LAKE_VAPFLX, and OUT_LAKE_VAPFLX_V.	TJB
  2010-Nov-02 Changed units of lake_var moisture fluxes to volume (m3).	TJB
  2010-Nov-21 Added OUT_LAKE_DSTOR, OUT_LAKE_DSTOR_V, OUT_LAKE_DSWE,
	      OUT_LAKE_DSWE_V, OUT_LAKE_SWE, and OUT_LAKE_SWE_V.	TJB
  2010-Nov-26 Added lake->ice_throughfall, lake->swe_save, and
	      lake->volume_save.					TJB
  2010-Dec-01 Added cell->zwt and OUT_ZWT.				TJB
  2011-Jan-04 Added zwtvmoist_zwt and zwtvmoist_moist to soil_con_struct
	      for storing relationship between soil moisture and water
	      water table position based on soil water retention curve.	TJB
  2011-Mar-01 Added cell->zwt2, cell->zwt3, OUT_ZWT2, OUT_ZWT3, and
	      OUT_ZWTL.							TJB
  2011-May-31 Removed options.GRND_FLUX.				TJB
  2011-May-31 Increased length of zwtvmoist_* arrays.			TJB
  2011-Jun-03 Added options.ORGANIC_FRACT.  Soil properties now take
	      organic fraction into account.				TJB
  2011-Nov-04 Added options to handle new forcing estimation features.	TJB
  2012-Jan-02 Removed wetland_veg_class from lake_con_struct.		TJB
  2012-Jan-16 Removed LINK_DEBUG code					BN
  2012-Jan-28 Removed AR_COMBO and GF_FULL.				TJB
  2012-Feb-07 Removed OUT_ZWT2 and OUT_ZWTL; renamed OUT_ZWT3 to
	      OUT_ZWT_LUMPED.						TJB
  2012-Feb-08 Renamed depth_full_snow_cover to max_snow_distrib_slope
	      and clarified the descriptions of the SPATIAL_SNOW
	      option.							TJB
  2012-Mar-30 Created constant DEFAULT_WIND_SPEED.			TJB
  2013-Jul-25 Added CATM, COSZEN, FDIR, PAR, OUT_CATM, OUT_COSZEN,
	      OUT_FDIR, and OUT_PAR.					TJB
  2013-Jul-25 Added photosynthesis terms.				TJB
  2013-Jul-25 Added soil carbon terms.					TJB
  2013-Jul-25 Added SLAB_MOIST.						TJB
  2013-Jul-25 Implemented heat flux between lake and soil.		TJB
  2013-Jul-25 Added DIST_ZWT terms.					TJB
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

/***** Aerodynamic Resistance options *****/
#define AR_406      0
#define AR_406_LS   1
#define AR_406_FULL 2
#define AR_410      3

/***** Ground Flux options *****/
#define GF_406  0
#define GF_410  1

/***** VP algorithm options *****/
#define VP_ITER_NONE     0
#define VP_ITER_ALWAYS   1
#define VP_ITER_ANNUAL   2
#define VP_ITER_CONVERGE 3

/***** Longwave Clear-Sky Algorithm options *****/
#define LW_TVA        0
#define LW_ANDERSON   1
#define LW_BRUTSAERT  2
#define LW_SATTERLUND 3
#define LW_IDSO       4
#define LW_PRATA      5

/***** Longwave Cloud Algorithm options *****/
#define LW_CLOUD_BRAS       0
#define LW_CLOUD_DEARDORFF  1

/***** Potential Evap types *****/
#define N_PET_TYPES 6
#define N_PET_TYPES_NON_NAT 4
#define PET_SATSOIL 0
#define PET_H2OSURF 1
#define PET_SHORT   2
#define PET_TALL    3
#define N_PET_TYPES_NAT 2
#define PET_NATVEG  4
#define PET_VEGNOCR 5

/***** LAI source types *****/
#define LAI_FROM_VEGLIB     0
#define LAI_FROM_VEGPARAM   1

/***** Canopy resistance parametrizations *****/
#define RC_JARVIS 0
#define RC_PHOTO  1

/***** Photosynthesis parametrizations *****/
#define PS_FARQUHAR 1
#define PS_MONTEITH 2

/***** Photosynthetic pathways  *****/
#define PHOTO_C3 0
#define PHOTO_C4 1

/***** Hard-coded veg class parameters (mainly for pot_evap) *****/
#define BARE_SOIL_ALBEDO 0.2	    /* albedo for bare soil */
#define H2O_SURF_ALBEDO 0.08	    /* albedo for deep water surface */
extern char   ref_veg_over[];
extern double ref_veg_rarc[];
extern double ref_veg_rmin[];
extern double ref_veg_lai[];
extern double ref_veg_albedo[];
extern double ref_veg_rough[];
extern double ref_veg_displ[];
extern double ref_veg_wind_h[];
extern double ref_veg_RGL[];
extern double ref_veg_rad_atten[];
extern double ref_veg_wind_atten[];
extern double ref_veg_trunk_ratio[];
extern char ref_veg_ref_crop[];

/***** Time Constants *****/
#define DAYS_PER_YEAR 365.
#define HOURSPERDAY   24        /* number of hours per day */
#define HOURSPERYEAR  24*365    /* number of hours per year */
#define SECPHOUR     3600	/* seconds per hour */
#define SEC_PER_DAY 86400.	/* seconds per day */

/***** Physical Constants *****/
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
#define Cp           1013.0	/* Specific heat at constant pressure of moist air 
				   (J/deg/K) (H.B.H. p.4.13)*/
#define CH_ICE       2100.0e3	/* Volumetric heat capacity (J/(m3*C)) of ice */
#define CH_WATER     4186.8e3   /* volumetric heat capacity of water */
#define K_SNOW       2.9302e-6  /* conductivity of snow (W/mK) */
#define K_AIR        2.32e-2    /* conductivity of air (W/mK) */
#define K_ICE        2.29       /* conductivity of ice (W/mK) */
#define SOLAR_CONSTANT 1400.0	/* Solar constant in W/m^2 */
#define EPS          0.62196351 /* Ratio of molecular weights: M_water_vapor/M_dry_air */
#define G            9.81       /* gravity */
#define Rd           287        /* Gas constant of dry air (J/degC*kg) */
#define Rgas         8.3143     /* [m3 Pa mol-1 K-1] universal gas law constant */
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

/***** Carbon Cycling constants *****/

#define CatmCurrent  383       /* Current global atmospheric CO2 mixing ratio (ppm) */
#define SW2PAR       0.45      /* Empirical ratio of PAR [W/m2] to SHORTWAVE
                                 [W/m2] from Lopez et al., 2001 */

/* Photosynthesis Parameters */
#define OMEGA        0.12      /* single leaf scattering albedo */
#define LaiMax       8         /* Maximum LAI in nitrogen scaling */
#define LaiLimit     3         /* Minimum LAI in nitrogen scaling and
                                  maximum LAI in PAR computation */
#define LaiMin       1e-9      /* Minimum LAI in PAR computation */
#define Epar         2.2e5     /* Energy content of PAR [J/mol photons]
                                  = (4.6 mol/MJ PAR)^-1 */
#define FcMax        0.9       /* Maximum fractional veg cover;
                                  (1-FcMax) = min amount of ground visible */
#define FcMin        1e-3      /* Minimum fractional veg cover;
                                  (1-FcMin) = max amount of ground visible */
#define ZenithMin    0.0174524 /* Check for solar zenith angle > 89 deg */
#define ZenithMinPar 1e-3      /* Cosine of the minimum solar zenith
                                  angle for photosynthesis to take place */
#define AlbSoiParMin 0.0       /* Minimum soil reflectivity in PAR range */
#define minMaxETrans 1e-12     /* Minimum of maximum electron transport
                                  rate [10e-12 mol/(m^2 s)] */
#define minStomCond  0.0       /* Minimum stomatal conductance [mol H2O/m2s] */
#define MCO2kg       44.011e-3 /* Molar mass of CO2 (kg/mol) */
#define MAIRkg       28.97e-3  /* Molar mass of air (kg/mol) */
#define MCg          12.01     /* Molar mass of C (g/mol) */
/* Factors that relate leaf internal CO2 concentration to ambient CO2 concentration */
#define FCI1C3       0.87      /* C3 Plants */
#define FCI1C4       0.67      /* C4 Plants */

/* C3 PLANTS: FARQUHAR, G.D., S. VON CAEMMERER AND J.A. BERRY, 1980.
              A BIOCHEMICAL MODEL OF PHOTOYNTHESIS IN LEAVES OF C3 SPECIES.
              PLANTA 149, 78-90. */

#define OX     0.21    /* OXYGEN CONCENTRATION [MOL(O2) / MOL(AIR)] */
#define KC0    460.e-6 /* MICHAELIS-MENTEN CONSTANT FOR CO2 AT 25C [MOL(CO2) / MOL(AIR)] */
#define KO0    330.e-3 /* MICHAELIS-MENTEN CONSTANT FOR O2 AT 25C [MOL(O2) / MOL(AIR)] */
#define EC     59356.  /* ACTIVATION ENERGY FOR KC [J / MOL] */
#define EO     35948.  /* ACTIVATION ENERGY FOR KO [J / MOL] */
#define EV     58520.  /* ACTIVATION ENERGY FOR VCMAX [J / MOL] */
#define ER     45000.  /* ACTIVATION ENERGY FOR DARK RESPIRATION [J / MOL] */
#define ALC3   0.28    /* EFFICIENCY OF OF PHOTON CAPTURE */
#define FRDC3  0.011   /* RATIO OF DARK RESPIRATION TO "PVM" AT 25C for C3 */

/* C4 PLANTS: COLLATZ, G.J., M. RIBAS-CARBO AND J.A. BERRY, 1992.
              COUPLED PHOTOSYNTHESIS-STOMATAL CONDUCTANCE MODEL FOR LEAVES
              OF C4 PLANTS. AUST. J. PLANT PHYSIOL. 19, 519-538. */

#define EK     50967.  /*  = Q10=2 (Collatz et al. 1992) */
#define ALC4   0.04    /* EFFECTIVE QUANTUM EFFICIENCY */
#define FRDC4  0.042   /* RATIO OF DARK RESPIRATION TO "PVM" AT 25C for C4 */
#define THETA  0.83    /* CURVATURE PARAMETER */

/* Plant Maintenance and Growth Respiration Parameters */
#define FRLeaf   0.4   /* Ratio of canopy leaf respiration to whole plant maintenance respiration */
#define FRGrowth 0.25  /* Ratio of plant growth respiration to NPP */

/* Soil Respiration Parameters */
#define E0_LT	308.56 /* Lloyd-Taylor E0 parameter [K] */
#define T0_LT	227.13 /* Lloyd-Taylor T0 parameter [K] */
#define wminFM  0.0    /* minimum soil moisture (fraction) at which soil respiration can occur */
#define wmaxFM  1.0    /* maximum soil moisture (fraction) at which soil respiration can occur */
#define woptFM  0.5    /* soil moisture (fraction) at which maximum soil respiration occurs */
#define Rhsat   0.15   /* ratio of soil respiration rate under saturated conditions (w=wmaxFM) to that under optimal conditions (w=woptFM) */
#define Rfactor 0.5    /* scaling factor to account for other (non-moisture) sources of inhibition of respiration */
#define tauLitter 2.86 /* Litter pool turnover time [y] */
#define tauInter  33.3 /* Intermediate pool turnover time [y] */
#define tauSlow   1000 /* Slow pool turnover time [y] */
#define fAir    0.7    /* Fraction of respired carbon from litter pool that is lost to atmosphere */
#define fInter  0.985  /* Fraction of [respired carbon from litter pool that goes to soil] that goes to intermediate pool */

/***** Wetland Microtopography *****/
#define RIDGE_FRACT	0.5	/* area fraction of ridges on the land surface */
#define RIDGE_HEIGHT	50	/* center height (cm) of ridges on the land surface */
#define HOLLOW_DEPTH	20	/* center depth (cm) of hollows on the land surface */

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
#define DEFAULT_WIND_SPEED 3.0  /* Default wind speed [m/s] used when wind is not supplied as a forcing */
#define SLAB_MOIST_FRACT 1.0    /* Ratio of the moisture in the soil/rock below the bottom soil layer to bottom soil layer moisture */

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
#define N_FORCING_TYPES 27
#define AIR_TEMP   0 /* air temperature per time step [C] (ALMA_INPUT: [K]) */
#define ALBEDO     1 /* surface albedo [fraction] */
#define CATM       2 /* atmospheric CO2 concentration [ppm] */
#define CHANNEL_IN 3 /* incoming channel flow [m3] (ALMA_INPUT: [m3/s]) */
#define CRAINF     4 /* convective rainfall [mm] (ALMA_INPUT: [mm/s]) */
#define CSNOWF     5 /* convective snowfall [mm] (ALMA_INPUT: [mm/s]) */
#define DENSITY    6 /* atmospheric density [kg/m3] */
#define FDIR       7 /* fraction of incoming shortwave that is direct [fraction] */
#define LONGWAVE   8 /* incoming longwave radiation [W/m2] */
#define LSRAINF    9 /* large-scale rainfall [mm] (ALMA_INPUT: [mm/s]) */
#define LSSNOWF   10 /* large-scale snowfall [mm] (ALMA_INPUT: [mm/s]) */
#define PAR       11 /* incoming photosynthetically active radiation [W/m2] */
#define PREC      12 /* total precipitation (rain and snow) [mm] (ALMA_INPUT: [mm/s]) */
#define PRESSURE  13 /* atmospheric pressure [kPa] (ALMA_INPUT: [Pa]) */
#define QAIR      14 /* specific humidity [kg/kg] */
#define RAINF     15 /* rainfall (convective and large-scale) [mm] (ALMA_INPUT: [mm/s]) */
#define REL_HUMID 16 /* relative humidity [fraction] */
#define SHORTWAVE 17 /* incoming shortwave [W/m2] */
#define SNOWF     18 /* snowfall (convective and large-scale) [mm] (ALMA_INPUT: [mm/s]) */
#define TMAX      19 /* maximum daily temperature [C] (ALMA_INPUT: [K]) */
#define TMIN      20 /* minimum daily temperature [C] (ALMA_INPUT: [K]) */
#define TSKC      21 /* cloud cover fraction [fraction] */
#define VP        22 /* vapor pressure [kPa] (ALMA_INPUT: [Pa]) */
#define WIND      23 /* wind speed [m/s] */
#define WIND_E    24 /* zonal component of wind speed [m/s] */
#define WIND_N    25 /* meridional component of wind speed [m/s] */
#define SKIP      26 /* place holder for unused data columns */

/***** Output Variable Types *****/
#define N_OUTVAR_TYPES 190
// Water Balance Terms - state variables
#define OUT_ASAT             0  /* Saturated Area Fraction */
#define OUT_DISTZWT_ASAT     1  /* Array of saturated area fraction values when DIST_ZWT = TRUE; order of values = [veg][band][wt] */
#define OUT_DISTZWT_SMOIST   2  /* Array of soil moisture values [mm] when DIST_ZWT = TRUE; order of values = [veg][band][wt][layer] */
#define OUT_DISTZWT_ZWT      3  /* Array of water table positions [cm] (zwt within lowest unsaturated layer) when DIST_ZWT = TRUE; order of values = [veg][band][wt] */
#define OUT_DISTZWT_ZWT_LUMP 4  /* Array of water table positions [cm] (zwt of total moisture across all layers, lumped together) when DIST_ZWT = TRUE; order of values = [veg][band][wt] - method 2 */
#define OUT_LAKE_AREA_FRAC   5  /* lake surface area as fraction of the grid cell area [fraction] */
#define OUT_LAKE_DEPTH       6  /* lake depth (distance between surface and deepest point) [m] */
#define OUT_LAKE_ICE         7  /* moisture stored as lake ice [mm over lake ice area] */
#define OUT_LAKE_ICE_FRACT   8  /* fractional coverage of lake ice [fraction] */
#define OUT_LAKE_ICE_HEIGHT  9  /* thickness of lake ice [cm] */
#define OUT_LAKE_MOIST      10  /* liquid water and ice stored in lake [mm over grid cell] */
#define OUT_LAKE_SURF_AREA  11  /* lake surface area [m2] */
#define OUT_LAKE_SWE        12  /* liquid water equivalent of snow on top of lake ice [m over lake ice area] */
#define OUT_LAKE_SWE_V      13  /* volumetric liquid water equivalent of snow on top of lake ice [m3] */
#define OUT_LAKE_VOLUME     14  /* lake volume [m3] */
#define OUT_ROOTMOIST       15  /* root zone soil moisture  [mm] */
#define OUT_SMFROZFRAC      16  /* fraction of soil moisture (by mass) that is ice, for each soil layer */
#define OUT_SMLIQFRAC       17  /* fraction of soil moisture (by mass) that is liquid, for each soil layer */
#define OUT_SNOW_CANOPY     18  /* snow interception storage in canopy  [mm] */
#define OUT_SNOW_COVER      19  /* fractional area of snow cover [fraction] */
#define OUT_SNOW_DEPTH      20  /* depth of snow pack [cm] */
#define OUT_SOIL_ICE        21  /* soil ice content  [mm] for each soil layer */
#define OUT_SOIL_LIQ        22  /* soil liquid content  [mm] for each soil layer */
#define OUT_SOIL_MOIST      23  /* soil total moisture content  [mm] for each soil layer */
#define OUT_SOIL_WET        24  /* vertical average of (soil moisture - wilting point)/(maximum soil moisture - wilting point) [mm/mm] */
#define OUT_SURFSTOR        25  /* storage of liquid water and ice (not snow) on surface (ponding) [mm] */
#define OUT_SURF_FROST_FRAC 26  /* fraction of soil surface that is frozen [fraction] */
#define OUT_SWE             27  /* snow water equivalent in snow pack (including vegetation-intercepted snow)  [mm] */
#define OUT_WDEW            28  /* total moisture interception storage in canopy [mm] */
#define OUT_ZWT             29  /* water table position [cm] (zwt within lowest unsaturated layer) */
#define OUT_ZWT_LUMPED      30  /* lumped water table position [cm] (zwt of total moisture across all layers, lumped together) */
// Water Balance Terms - fluxes
#define OUT_BASEFLOW        31  /* baseflow out of the bottom layer  [mm] (ALMA_OUTPUT: [mm/s]) */
#define OUT_DELINTERCEPT    32  /* change in canopy interception storage  [mm] */
#define OUT_DELSOILMOIST    33  /* change in soil water content  [mm] */
#define OUT_DELSURFSTOR     34  /* change in surface liquid water storage  [mm] */
#define OUT_DELSWE          35  /* change in snow water equivalent  [mm] */
#define OUT_DISTZWT_BFLOW   36  /* Array of baseflow values [mm] when DIST_ZWT = TRUE; order of values = [veg][band][wt] */
#define OUT_DISTZWT_RUNOFF  37  /* Array of runoff values [mm] when DIST_ZWT = TRUE; order of values = [veg][band][wt] */
#define OUT_EVAP            38  /* total net evaporation [mm] (ALMA_OUTPUT: [mm/s]) */
#define OUT_EVAP_BARE       39  /* net evaporation from bare soil [mm] (ALMA_OUTPUT: [mm/s]) */
#define OUT_EVAP_CANOP      40  /* net evaporation from canopy interception [mm] (ALMA_OUTPUT: [mm/s]) */
#define OUT_INFLOW          41  /* moisture that reaches top of soil column [mm] (ALMA_OUTPUT: [mm/s]) */
#define OUT_LAKE_BF_IN      42  /* incoming baseflow from lake catchment [mm] (ALMA_OUTPUT: [mm/s]) */
#define OUT_LAKE_BF_IN_V    43  /* incoming volumetric baseflow from lake catchment [m3] (ALMA_OUTPUT: [m3/s]) */
#define OUT_LAKE_BF_OUT     44  /* outgoing baseflow from lake [mm] (ALMA_OUTPUT: [mm/s]) */
#define OUT_LAKE_BF_OUT_V   45  /* outgoing volumetric baseflow from lake [m3] (ALMA_OUTPUT: [m3/s]) */
#define OUT_LAKE_CHAN_IN    46  /* channel inflow into lake [mm] (ALMA_OUTPUT: [mm/s]) */
#define OUT_LAKE_CHAN_IN_V  47  /* volumetric channel inflow into lake [m3] (ALMA_OUTPUT: [m3/s]) */
#define OUT_LAKE_CHAN_OUT   48  /* channel outflow from lake [mm] (ALMA_OUTPUT: [mm/s]) */
#define OUT_LAKE_CHAN_OUT_V 49  /* volumetric channel outflow from lake [m3] (ALMA_OUTPUT: [m3/s]) */
#define OUT_LAKE_DSTOR      50  /* change in lake moisture storage (liquid plus ice cover) [mm] (ALMA_OUTPUT: [mm/s]) */
#define OUT_LAKE_DSTOR_V    51  /* volumetric change in lake moisture storage (liquid plus ice cover) [m3] (ALMA_OUTPUT: [m3/s]) */
#define OUT_LAKE_DSWE       52  /* change in swe on top of lake ice [mm] (ALMA_OUTPUT: [mm/s]) */
#define OUT_LAKE_DSWE_V     53  /* volumetric change in swe on top of lake ice [m3] (ALMA_OUTPUT: [m3/s]) */
#define OUT_LAKE_EVAP       54  /* net evaporation from lake surface [mm] (ALMA_OUTPUT: [mm/s]) */
#define OUT_LAKE_EVAP_V     55  /* net volumetric evaporation from lake surface [m3] (ALMA_OUTPUT: [m3/s]) */
#define OUT_LAKE_PREC_V     56  /* volumetric precipitation over lake surface [m3] (ALMA_OUTPUT: [m3/s]) */
#define OUT_LAKE_RCHRG      57  /* recharge from lake to surrounding wetland [mm] (ALMA_OUTPUT: [mm/s]) */
#define OUT_LAKE_RCHRG_V    58  /* volumetric recharge from lake to surrounding wetland [m3] (ALMA_OUTPUT: [m3/s]) */
#define OUT_LAKE_RO_IN      59  /* incoming runoff from lake catchment [mm] (ALMA_OUTPUT: [mm/s]) */
#define OUT_LAKE_RO_IN_V    60  /* incoming volumetric runoff from lake catchment [m3] (ALMA_OUTPUT: [m3/s]) */
#define OUT_LAKE_VAPFLX     61  /* outgoing sublimation from snow on top of lake ice [mm] (ALMA_OUTPUT: [mm/s]) */
#define OUT_LAKE_VAPFLX_V   62  /* outgoing volumetric sublimation from snow on top of lake ice [m3] (ALMA_OUTPUT: [m3/s]) */
#define OUT_PET_SATSOIL     63  /* potential evap from saturated bare soil [mm] (ALMA_OUTPUT: [mm/s]) */
#define OUT_PET_H2OSURF     64  /* potential evap from open water [mm] (ALMA_OUTPUT: [mm/s]) */
#define OUT_PET_SHORT       65  /* potential evap (transpiration only) from short reference crop (grass) [mm] (ALMA_OUTPUT: [mm/s]) */
#define OUT_PET_TALL        66  /* potential evap (transpiration only) from tall reference crop (alfalfa) [mm] (ALMA_OUTPUT: [mm/s]) */
#define OUT_PET_NATVEG      67  /* potential evap (transpiration only) from current vegetation and current canopy resistance [mm] (ALMA_OUTPUT: [mm/s]) */
#define OUT_PET_VEGNOCR     68  /* potential evap (transpiration only) from current vegetation and 0 canopy resistance [mm] (ALMA_OUTPUT: [mm/s]) */
#define OUT_PREC            69  /* incoming precipitation [mm] (ALMA_OUTPUT: [mm/s]) */
#define OUT_RAINF           70  /* rainfall  [mm] (ALMA_OUTPUT: [mm/s]) */
#define OUT_REFREEZE        71  /* refreezing of water in the snow  [mm] (ALMA_OUTPUT: [mm/s]) */
#define OUT_RUNOFF          72  /* surface runoff [mm] (ALMA_OUTPUT: [mm/s]) */
#define OUT_SNOW_MELT       73  /* snow melt  [mm] (ALMA_OUTPUT: [mm/s]) */
#define OUT_SNOWF           74  /* snowfall  [mm] (ALMA_OUTPUT: [mm/s]) */
#define OUT_SUB_BLOWING     75  /* net sublimation of blowing snow [mm] (ALMA_OUTPUT: [mm/s]) */
#define OUT_SUB_CANOP       76  /* net sublimation from snow stored in canopy [mm] (ALMA_OUTPUT: [mm/s]) */
#define OUT_SUB_SNOW        77  /* total net sublimation from snow pack (surface and blowing) [mm] (ALMA_OUTPUT: [mm/s]) */
#define OUT_SUB_SURFACE     78  /* net sublimation from snow pack surface [mm] (ALMA_OUTPUT: [mm/s]) */
#define OUT_TRANSP_VEG      79  /* net transpiration from vegetation [mm] (ALMA_OUTPUT: [mm/s]) */
#define OUT_WATER_ERROR     80  /* water budget error [mm] */
// Energy Balance Terms - state variables
#define OUT_ALBEDO          81  /* average surface albedo [fraction] */
#define OUT_BARESOILT       82  /* bare soil surface temperature [C] (ALMA_OUTPUT: [K]) */
#define OUT_FDEPTH          83  /* depth of freezing fronts [cm] (ALMA_OUTPUT: [m]) for each freezing front */
#define OUT_LAKE_ICE_TEMP   84  /* temperature of lake ice [C] (ALMA_OUTPUT: [K]) */
#define OUT_LAKE_SURF_TEMP  85  /* lake surface temperature [C] (ALMA_OUTPUT: [K]) */
#define OUT_LAKE_LAYER_TEMP 86  /* lake water layer temperature [C] (ALMA_OUTPUT: [K]) */
#define OUT_RAD_TEMP        87  /* average radiative surface temperature [K] */
#define OUT_SALBEDO         88  /* snow pack albedo [fraction] */
#define OUT_SNOW_PACK_TEMP  89  /* snow pack temperature [C] (ALMA_OUTPUT: [K]) */
#define OUT_SNOW_SURF_TEMP  90  /* snow surface temperature [C] (ALMA_OUTPUT: [K]) */
#define OUT_SNOWT_FBFLAG    91  /* snow surface temperature fallback flag */
#define OUT_SOIL_TEMP       92  /* soil temperature [C] (ALMA_OUTPUT: [K]) for each soil layer */
#define OUT_SOIL_TEMP_LAKE  93  /* soil temperature [C] (ALMA_OUTPUT: [K]) for each soil layer under the lake */
#define OUT_SOIL_TNODE      94  /* soil temperature [C] (ALMA_OUTPUT: [K]) for each soil thermal node */
#define OUT_SOIL_TNODE_LAKE 95  /* soil temperature [C] (ALMA_OUTPUT: [K]) for each soil thermal node under the lake */
#define OUT_SOIL_TNODE_WL   96  /* soil temperature [C] (ALMA_OUTPUT: [K]) for each soil thermal node in the wetland */
#define OUT_SOILT_FBFLAG    97  /* soil temperature flag for each soil thermal node */
#define OUT_SURF_TEMP       98  /* average surface temperature [C] (ALMA_OUTPUT: [K]) */
#define OUT_SURFT_FBFLAG    99  /* surface temperature flag */
#define OUT_TCAN_FBFLAG    100  /* Tcanopy flag */
#define OUT_TDEPTH         101  /* depth of thawing fronts [cm] (ALMA_OUTPUT: [m]) for each thawing front */
#define OUT_TFOL_FBFLAG    102  /* Tfoliage flag */
#define OUT_VEGT           103  /* average vegetation canopy temperature [C] (ALMA_OUTPUT: [K]) */
// Energy Balance Terms - fluxes
#define OUT_ADV_SENS       104  /* net sensible flux advected to snow pack [W/m2] */
#define OUT_ADVECTION      105  /* advected energy [W/m2] */
#define OUT_DELTACC        106  /* rate of change in cold content in snow pack [W/m2] (ALMA_OUTPUT: [J/m2]) */
#define OUT_DELTAH         107  /* rate of change in heat storage [W/m2] (ALMA_OUTPUT: [J/m2]) */
#define OUT_ENERGY_ERROR   108  /* energy budget error [W/m2] */
#define OUT_FUSION         109  /* net energy used to melt/freeze soil moisture [W/m2] */
#define OUT_GRND_FLUX      110  /* net heat flux into ground [W/m2] */
#define OUT_IN_LONG        111  /* incoming longwave at ground surface (under veg) [W/m2] */
#define OUT_LATENT         112  /* net upward latent heat flux [W/m2] */
#define OUT_LATENT_SUB     113  /* net upward latent heat flux from sublimation [W/m2] */
#define OUT_MELT_ENERGY    114  /* energy of fusion (melting) in snowpack [W/m2] */
#define OUT_NET_LONG       115  /* net downward longwave flux [W/m2] */
#define OUT_NET_SHORT      116  /* net downward shortwave flux [W/m2] */
#define OUT_R_NET          117  /* net downward radiation flux [W/m2] */
#define OUT_RFRZ_ENERGY    118  /* net energy used to refreeze liquid water in snowpack [W/m2] */
#define OUT_SENSIBLE       119  /* net upward sensible heat flux [W/m2] */
#define OUT_SNOW_FLUX      120  /* energy flux through snow pack [W/m2] */
// Miscellaneous Terms
#define OUT_AERO_COND      121  /* "scene" aerodynamic conductance [m/s] (tiles with overstory contribute overstory conductance; others contribute surface conductance) */
#define OUT_AERO_COND1     122  /* surface aerodynamic conductance [m/s] */
#define OUT_AERO_COND2     123  /* overstory aerodynamic conductance [m/s] */
#define OUT_AERO_RESIST    124  /* "scene"canopy aerodynamic resistance [s/m]  (tiles with overstory contribute overstory resistance; others contribute surface resistance)*/
#define OUT_AERO_RESIST1   125  /* surface aerodynamic resistance [s/m] */
#define OUT_AERO_RESIST2   126  /* overstory aerodynamic resistance [s/m] */
#define OUT_AIR_TEMP       127  /* air temperature [C] (ALMA_OUTPUT: [K])*/
#define OUT_CATM           128  /* atmospheric CO2 concentrtaion [ppm]*/
#define OUT_COSZEN         129  /* cosine of solar zenith angle [fraction]*/
#define OUT_DENSITY        130  /* near-surface atmospheric density [kg/m3]*/
#define OUT_FDIR           131  /* fraction of incoming shortwave that is direct [fraction]*/
#define OUT_LONGWAVE       132  /* incoming longwave [W/m2] */
#define OUT_PAR            133  /* incoming photosynthetically active radiation [W/m2] */
#define OUT_PRESSURE       134  /* near surface atmospheric pressure [kPa] (ALMA_OUTPUT: [Pa])*/
#define OUT_QAIR           135  /* specific humidity [kg/kg] */
#define OUT_REL_HUMID      136  /* relative humidity [fraction]*/
#define OUT_SHORTWAVE      137  /* incoming shortwave [W/m2] */
#define OUT_SURF_COND      138  /* surface conductance [m/s] */
#define OUT_TSKC           139  /* cloud cover fraction [fraction] */
#define OUT_VP             140  /* near surface vapor pressure [kPa] (ALMA_OUTPUT: [Pa]) */
#define OUT_VPD            141  /* near surface vapor pressure deficit [kPa] (ALMA_OUTPUT: [Pa]) */
#define OUT_WIND           142  /* near surface wind speed [m/s] */
// Band-specific quantities
#define OUT_ADV_SENS_BAND       143  /* net sensible heat flux advected to snow pack [W/m2] */
#define OUT_ADVECTION_BAND      144  /* advected energy [W/m2] */
#define OUT_ALBEDO_BAND         145  /* average surface albedo [fraction] */
#define OUT_DELTACC_BAND        146  /* change in cold content in snow pack [W/m2] */
#define OUT_GRND_FLUX_BAND      147  /* net heat flux into ground [W/m2] */
#define OUT_IN_LONG_BAND        148  /* incoming longwave at ground surface (under veg) [W/m2] */
#define OUT_LATENT_BAND         149  /* net upward latent heat flux [W/m2] */
#define OUT_LATENT_SUB_BAND     150  /* net upward latent heat flux due to sublimation [W/m2] */
#define OUT_MELT_ENERGY_BAND    151  /* energy of fusion (melting) in snowpack [W/m2] */
#define OUT_NET_LONG_BAND       152  /* net downward longwave flux [W/m2] */
#define OUT_NET_SHORT_BAND      153  /* net downward shortwave flux [W/m2] */
#define OUT_RFRZ_ENERGY_BAND    154  /* net energy used to refreeze liquid water in snowpack [W/m2] */
#define OUT_SENSIBLE_BAND       155  /* net upward sensible heat flux [W/m2] */
#define OUT_SNOW_CANOPY_BAND    156  /* snow interception storage in canopy [mm] */
#define OUT_SNOW_COVER_BAND     157  /* fractional area of snow cover [fraction] */
#define OUT_SNOW_DEPTH_BAND     158  /* depth of snow pack [cm] */
#define OUT_SNOW_FLUX_BAND      159  /* energy flux through snow pack [W/m2] */
#define OUT_SNOW_MELT_BAND      160  /* snow melt [mm] (ALMA_OUTPUT: [mm/s]) */
#define OUT_SNOW_PACKT_BAND     161  /* snow pack temperature [C] (ALMA_OUTPUT: [K]) */
#define OUT_SNOW_SURFT_BAND     162  /* snow surface temperature [C] (ALMA_OUTPUT: [K]) */
#define OUT_SWE_BAND            163  /* snow water equivalent in snow pack [mm] */
// Dynamic Soil Property Terms - EXCESS_ICE option
#if EXCESS_ICE
#define OUT_SOIL_DEPTH          164  /* soil moisture layer depths [m] */
#define OUT_SUBSIDENCE          165  /* subsidence of soil layer [mm] */
#define OUT_POROSITY            166  /* porosity [mm/mm] */
#define OUT_ZSUM_NODE           167  /* depths of thermal nodes [m] */
#endif // EXCESS_ICE
// Carbon-Cycling Terms
#define OUT_APAR           168  /* absorbed PAR [W/m2] */
#define OUT_GPP            169  /* gross primary productivity [g C/m2d] */
#define OUT_RAUT           170  /* autotrophic respiration [g C/m2d] */
#define OUT_NPP            171  /* net primary productivity [g C/m2d] */
#define OUT_LITTERFALL     172  /* flux of carbon from living biomass into soil [g C/m2d] */
#define OUT_RHET           173  /* soil respiration (heterotrophic respiration) [g C/m2d] */
#define OUT_NEE            174  /* net ecosystem exchange (=NPP-RHET) [g C/m2d] */
#define OUT_CLITTER        175  /* Carbon density in litter pool [g C/m2d] */
#define OUT_CINTER         176  /* Carbon density in intermediate pool [g C/m2d] */
#define OUT_CSLOW          177  /* Carbon density in slow pool [g C/m2d] */

// Miscellaneous
#define OUT_DISTZWT_AREA   178  /* area fraction of each point in zwt distribution */

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
  int    AboveTreelineVeg; /* Default veg type to use above treeline;
			      Negative number indicates bare soil. */
  char   AERO_RESIST_CANSNOW; /* "AR_406" = multiply aerodynamic resistance
					    by 10 for latent heat but not
					    for sensible heat (as in
					    VIC 4.0.6); do NOT apply stability
					    correction; use surface aero_resist
					    for ET when no snow in canopy.
				 "AR_406_LS" = multiply aerodynamic resistance
					    by 10 for BOTH latent heat AND
					    sensible heat; do NOT apply
					    stability correction;
					    use surface aero_resist
					    for ET when no snow in canopy.
				 "AR_406_FULL" = multiply aerodynamic resistance
					    by 10 for BOTH latent heat AND
					    sensible heat; do NOT apply
					    stability correction;
					    always use canopy aero_resist
					    for ET.
				 "AR_410" = do not multiply aerodynamic
					    resistance by 10 in snow-filled
					    canopy (as in VIC 4.1.0);
					    DO apply stability correction;
					    always use canopy aero_resist
					    for ET. */
  char   BLOWING;        /* TRUE = calculate sublimation from blowing snow */
  char   CARBON;         /* TRUE = simulate carbon cycling processes;
			    FALSE = no carbon cycling (default) */
  char   COMPUTE_TREELINE; /* TRUE = Determine treeline and exclude overstory
			      vegetation from higher elevations */
  char   CONTINUEONERROR;/* TRUE = VIC will continue to run after a cell has an error */
  char   CORRPREC;       /* TRUE = correct precipitation for gage undercatch */
  char   DIST_PRCP;      /* TRUE = Use distributed precipitation model */
  char   DIST_ZWT;       /* TRUE = consider a spatially-distributed water table depth */
  char   EQUAL_AREA;     /* TRUE = RESOLUTION stores grid cell area in km^2;
			    FALSE = RESOLUTION stores grid cell side length in degrees */
  char   EXP_TRANS;      /* TRUE = Uses grid transform for exponential node 
			    distribution for soil heat flux calculations*/
  char   FROZEN_SOIL;    /* TRUE = Use frozen soils code */
  char   FULL_ENERGY;    /* TRUE = Use full energy code */
  char   GRND_FLUX_TYPE; /* "GF_406"  = use (flawed) formulas for ground flux, deltaH, and fusion
                                        from VIC 4.0.6 and earlier
                            "GF_410"  = use formulas from VIC 4.1.0 */
  char   IMPLICIT;       /* TRUE = Use implicit solution when computing 
			    soil thermal fluxes */
  char   JULY_TAVG_SUPPLIED; /* If TRUE and COMPUTE_TREELINE is also true,
			        then average July air temperature will be read
			        from soil file and used in calculating treeline */
  char   LAKES;          /* TRUE = use lake energy code */
  char   LW_CLOUD;       /* Longwave cloud formulation; "LW_CLOUD_x" = code for LW cloud formulation - see LW_CLOUD codes above */
  char   LW_TYPE;        /* Longwave clear sky algorithm; "LW_x" = code for LW algorithm - see LW codes above */
  float  MIN_WIND_SPEED; /* Minimum wind speed in m/s that can be used by the model. **/
  char   MTCLIM_SWE_CORR;/* TRUE = correct MTCLIM's downward shortwave radiation estimate in presence of snow */
  int    Ncanopy;        /* Number of canopy layers in the model. */
  int    Nlayer;         /* Number of layers in model */
  int    Nnode;          /* Number of soil thermal nodes in the model */
  int    Nzwt;           /* Number of bins in the water table distribution */
  char   NOFLUX;         /* TRUE = Use no flux lower bondary when computing 
			    soil thermal fluxes */
  char   PLAPSE;         /* TRUE = If air pressure not supplied as an
			    input forcing, compute it by lapsing sea-level
			    pressure by grid cell average elevation;
			    FALSE = air pressure set to constant 95.5 kPa */
  float  PREC_EXPT;      /* Exponential that controls the fraction of a
			    grid cell that receives rain during a storm
			    of given intensity */
  char   RC_MODE;        /* RC_JARVIS = compute canopy resistance via Jarvis formulation (default)
                            RC_PHOTO = compute canopy resistance based on photosynthetic activity */
  int    ROOT_ZONES;     /* Number of root zones used in simulation */
  char   QUICK_FLUX;     /* TRUE = Use Liang et al., 1999 formulation for
			    ground heat flux, if FALSE use explicit finite
			    difference method */
  char   QUICK_SOLVE;    /* TRUE = Use Liang et al., 1999 formulation for 
			    iteration, but explicit finite difference
			    method for final step. */
  char   SHARE_LAYER_MOIST; /* TRUE = transpiration in moisture-limited layers can draw from other layers (default) */
  char   SNOW_ALBEDO;    /* USACE: Use algorithm of US Army Corps of Engineers, 1956; SUN1999: Use algorithm of Sun et al., JGR, 1999 */
  char   SNOW_DENSITY;   /* DENS_BRAS: Use algorithm of Bras, 1990; DENS_SNTHRM: Use algorithm of SNTHRM89 adapted for 1-layer pack */
  int    SNOW_BAND;      /* Number of elevation bands over which to solve the 
			    snow model */
  int    SNOW_STEP;      /* Time step in hours to use when solving the 
			    snow model */
  float  SW_PREC_THRESH; /* Minimum daily precipitation [mm] that can cause "dimming" of incoming shortwave radiation */
  char   TFALLBACK;      /* TRUE = when any temperature iterations fail to converge,
                                   use temperature from previous time step; the number
                                   of instances when this occurs will be logged and
                                   reported at the end of the cell's simulation
                            FALSE = when iterations fail to converge, report an error
                                    and abort simulation for current grid cell
                            Default = TRUE */
  char   VP_INTERP;      /* How to disaggregate VP from daily to sub-daily;
                            TRUE = linearly interpolate between daily VP values, assuming they occur at the times of Tmin;
                            FALSE = hold VP constant at the daily value */
  char   VP_ITER;        /* VP_ITER_NONE = never iterate with SW
                            VP_ITER_ALWAYS = always iterate with SW
                            VP_ITER_ANNUAL = use annual Epot/PRCP criterion
                            VP_ITER_CONVERGE = always iterate until convergence */
  int    Nlakenode;      /* Number of lake thermal nodes in the model. */

  // input options
  char   ALMA_INPUT;     /* TRUE = input variables are in ALMA-compliant units; FALSE = standard VIC units */
  char   ARC_SOIL;       /* TRUE = use ARC/INFO gridded ASCII files for soil 
			    parameters*/
  char   BASEFLOW;       /* ARNO: read Ds, Dm, Ws, c; NIJSSEN2001: read d1, d2, d3, d4 */
  int    GRID_DECIMAL;   /* Number of decimal places in grid file extensions */
  char   VEGLIB_PHOTO;   /* TRUE = veg library contains photosynthesis parameters */
  char   VEGPARAM_LAI;   /* TRUE = veg param file contains monthly LAI values */
  char   LAI_SRC;        /* LAI_FROM_VEGLIB = read LAI values from veg library file
                            LAI_FROM_VEGPARAM = read LAI values from the veg param file */
  char   LAKE_PROFILE;   /* TRUE = user-specified lake/area profile */
  char   ORGANIC_FRACT;  /* TRUE = organic matter fraction of each layer is read from the soil parameter file; otherwise set to 0.0. */

  // state options
  char   BINARY_STATE_FILE; /* TRUE = model state file is binary (default) */
  char   INIT_STATE;     /* TRUE = initialize model state from file */
  char   SAVE_STATE;     /* TRUE = save state file */       

  // output options
  char   ALMA_OUTPUT;    /* TRUE = output variables are in ALMA-compliant units; FALSE = standard VIC units */
  char   BINARY_OUTPUT;  /* TRUE = output files are in binary, not ASCII */
  char   COMPRESS;       /* TRUE = Compress all output files */
  char   MOISTFRACT;     /* TRUE = output soil moisture as fractional moisture content */
  int    Noutfiles;      /* Number of output files (not including state files) */
  char   PRT_HEADER;     /* TRUE = insert header at beginning of output file; FALSE = no header */
  char   PRT_SNOW_BAND;  /* TRUE = print snow parameters for each snow band. This is only used when default
				   output files are used (for backwards-compatibility); if outfiles and
				   variables are explicitly mentioned in global parameter file, this option
				   is ignored. */

} option_struct;

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
  float    AlbedoPar;                 /* soil albedo in PAR range (400-700nm) */
  double   alpha[MAX_NODES];          /* thermal solution constant */
  double   annual_prec;               /* annual average precipitation (mm) */
  double   avg_temp;                  /* average soil temperature (C) */
  double   avgJulyAirTemp;            /* Average July air temperature (C) */
  double   b_infilt;                  /* infiltration parameter */
  double   beta[MAX_NODES];           /* thermal solution constant */
  double   bubble[MAX_LAYERS];        /* bubbling pressure, HBH 5.15 (cm) */
  double   bubble_node[MAX_NODES];    /* bubbling pressure (cm) */
  double   bulk_density[MAX_LAYERS];  /* soil bulk density (kg/m^3) */
  double   bulk_dens_min[MAX_LAYERS]; /* bulk density of mineral soil (kg/m^3) */
  double   bulk_dens_org[MAX_LAYERS]; /* bulk density of organic soil (kg/m^3) */
  double   c;                         /* exponent in ARNO baseflow scheme */
  double   depth[MAX_LAYERS];         /* thickness of each soil moisture layer (m).  In the case of EXCESS_ICE, this is the effective (dynamic) depth. */
  double   dp;                        /* soil thermal damping depth (m) */
  double   dz_node[MAX_NODES];        /* thermal node thickness (m) */
  double   Zsum_node[MAX_NODES];      /* thermal node depth (m) */
  double   expt[MAX_LAYERS];          /* layer-specific exponent n (=3+2/lambda) in Campbell's eqn for hydraulic conductivity, HBH 5.6 */
  double   expt_node[MAX_NODES];      /* node-specific exponent n (=3+2/lambda) in Campbell's eqn for hydraulic conductivity, HBH 5.6 */
#if SPATIAL_FROST
  double   frost_fract[FROST_SUBAREAS]; /* spatially distributed frost coverage fractions */
  double   frost_slope;               // slope of frost distribution
#endif // SPATIAL_FROST
  double   gamma[MAX_NODES];          /* thermal solution constant */
  double   init_moist[MAX_LAYERS];    /* initial layer moisture level (mm) */
  double   max_infil;                 /* maximum infiltration rate */
  double   max_moist[MAX_LAYERS];     /* maximum moisture content (mm) per layer */
  double   max_moist_node[MAX_NODES]; /* maximum moisture content (mm/mm) per node */
#if SPATIAL_SNOW
  double   max_snow_distrib_slope;    /* Maximum slope of snow depth distribution [m].  This should equal 2*depth_min, where depth_min = minimum snow pack depth below which coverage < 1. */
#endif // SPATIAL_SNOW
  double   phi_s[MAX_LAYERS];         /* soil moisture diffusion parameter (mm/mm) */
  double   porosity[MAX_LAYERS];      /* porosity (fraction) */
  double   quartz[MAX_LAYERS];        /* quartz content of soil (fraction of mineral soil volume) */
  double   organic[MAX_LAYERS];       /* organic content of soil (fraction of total soil volume) */
  double   resid_moist[MAX_LAYERS];   /* residual moisture content of soil layer */
  double   rough;                     /* soil surface roughness (m) */
  double   snow_rough;                /* snow surface roughness (m) */
  double   soil_density[MAX_LAYERS];  /* soil particle density (kg/m^3) */
  double   soil_dens_min[MAX_LAYERS]; /* particle density of mineral soil (kg/m^3) */
  double   soil_dens_org[MAX_LAYERS]; /* particle density of organic soil (kg/m^3) */
  float   *BandElev;                  /* Elevation of each snow elevation band */
  double  *AreaFract;                 /* Fraction of grid cell included in each snow elevation band */
  double  *Pfactor;                   /* Change in Precipitation due to elevation (fract) in each snow elevation band */
  double  *Tfactor;                   /* Change in temperature due to elevation (C) in each snow elevation band */
  char    *AboveTreeLine;             /* Flag to indicate if band is above the treeline */
#if QUICK_FS
  double **ufwc_table_layer[MAX_LAYERS];
  double **ufwc_table_node[MAX_NODES]; 
#endif
  float    elevation;                 /* grid cell elevation (m) */
  float    lat;                       /* grid cell central latitude */
  float    lng;                       /* grid cell central longitude */
  double   cell_area;                 /* Area of grid cell (m^2) */
  float    time_zone_lng;             /* central meridian of the time zone */
  float  **layer_node_fract;          /* fraction of all nodes within each layer */
  int      gridcel;                   /* grid cell number */
  double   zwtvmoist_zwt[MAX_LAYERS+2][MAX_ZWTVMOIST]; /* zwt values in the zwt-v-moist curve for each layer */
  double   zwtvmoist_moist[MAX_LAYERS+2][MAX_ZWTVMOIST]; /* moist values in the zwt-v-moist curve for each layer */
  double   ridge_fract;               /* fractional area of ridges (of average height RIDGE_HEIGHT) on the land surface */
  double   ZwtAreaFract[MAX_NZWT];    /* array of cumulative fractional areas of bins of zwt distribution */
  double   ZwtDeltaMoist[MAX_NZWT];   /* array of moisture offsets of bins of zwt distribution */
  double   slope;
  double   aspect;
  double   ehoriz;
  double   whoriz;
#if EXCESS_ICE
  double   min_depth[MAX_LAYERS];     /* soil layer depth as given in the soil file (m).  The effective depth will always be >= this value. */
  double   porosity_node[MAX_NODES];  /* porosity for each thermal node */
  double   effective_porosity[MAX_LAYERS]; /* effective soil porosity (fraction) when soil pores are expanded due to excess ground ice */
  double   effective_porosity_node[MAX_NODES]; /* effective soil porosity (fraction) when soil pores are expanded due to excess ground ice */
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
  int     LAKE;             /* TRUE = this tile is a lake/wetland tile */
  double *CanopLayerBnd;    /* Upper boundary of each canopy layer, expressed as fraction of total LAI */
} veg_con_struct;

/******************************************************************
  This structure stores parameters for individual vegetation types.
  ******************************************************************/
typedef struct {
  char   overstory;        /* TRUE = overstory present, important for snow 
			      accumulation in canopy */
  double LAI[12];          /* leaf area index */
  double Wdmax[12];        /* maximum dew holding capacity (mm) */
  double albedo[12];       /* vegetation albedo (added for full energy) 
			      (fraction) */
  double displacement[12]; /* vegetation displacement height (m) */
  double emissivity[12];   /* vegetation emissivity (fraction) */
  int    NVegLibTypes;     /* number of vegetation classes defined in library */
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
  char   Ctype;            /* Photosynthetic pathway; can be C3 or C4 */
  double MaxCarboxRate;    /* maximum carboxlyation rate at 25 deg C (mol(CO2)/m2s) */
  double MaxETransport;    /* maximum electron transport rate at 25 deg C (mol(CO2)/m2s) (C3 plants) */
  double CO2Specificity;   /* CO2 specificity at 25 deg C (mol(CO2)/m2s) (C4 plants) */
  double LightUseEff;      /* Light-use efficiency (mol(CO2)/mol(photons)) */
  char   NscaleFlag;       /* TRUE = nitrogen-scaling factors are applicable
                              to this veg class */
  double Wnpp_inhib;       /* moisture level (fraction of maximum moisture) above which photosynthesis
                              experiencing saturation inhibition, i.e. too wet for optimal photosynthesis;
                              only applies to top soil layer */
  double NPPfactor_sat;    /* photosynthesis multiplier (fraction of maximum) when top soil
                              layer is saturated */
} veg_lib_struct;

/***************************************************************************
   This structure stores the atmospheric forcing data for each model time 
   step for a single grid cell.  Each array stores the values for the 
   SNOW_STEPs during the current model step and the value for the entire model
   step.  The latter is referred to by array[NR].  Looping over the SNOW_STEPs
   is done by for (i = 0; i < NF; i++) 
***************************************************************************/
typedef struct {
  double *air_temp;  /* air temperature (C) */
  double *Catm;      /* atmospheric CO2 mixing ratio (mol CO2/ mol air) */
  double *channel_in;/* incoming channel inflow for time step (mm) */
  double *coszen;    /* cosine of solar zenith angle (fraction) */
  double *density;   /* atmospheric density (kg/m^3) */
  double *fdir;      /* fraction of incoming shortwave that is direct (fraction) */
  double *longwave;  /* incoming longwave radiation (W/m^2) (net incoming
                        longwave for water balance model) */
  double out_prec;   /* Total precipitation for time step - accounts
                        for corrected precipitation totals */
  double out_rain;   /* Rainfall for time step (mm) */
  double out_snow;   /* Snowfall for time step (mm) */
  double *par;       /* incoming photosynthetically active radiation () */
  double *prec;      /* average precipitation in grid cell (mm) */
  double *pressure;  /* atmospheric pressure (kPa) */
  double *shortwave; /* incoming shortwave radiation (W/m^2) */
  char   *snowflag;  /* TRUE if there is snowfall in any of the snow
                        bands during the timestep, FALSE otherwise*/
  double *tskc;      /* cloud cover fraction (fraction) */
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
  double evap_dist_zwt[MAX_NZWT]; /* array of layer evapotranspiration values when DIST_ZWT = TRUE (mm) */
#if SPATIAL_FROST
  double ice[FROST_SUBAREAS]; /* ice content of the frozen sublayer (mm) */
  double ice_dist_zwt[MAX_NZWT][FROST_SUBAREAS]; /* array of layer ice contents when DIST_ZWT = TRUE (mm) */
#else
  double ice;               /* ice content of the frozen sublayer (mm) */
  double ice_dist_zwt[MAX_NZWT]; /* array of layer ice contents when DIST_ZWT = TRUE (mm) */
#endif
  double kappa;             /* average thermal conductivity of the current layer (W/m/K) */
  double moist;             /* moisture content of the unfrozen sublayer (mm) */
  double moist_dist_zwt[MAX_NZWT];   /* array of layer moisture contents when DIST_ZWT = TRUE (mm) */
  double phi;               /* moisture diffusion parameter */
  double zwt;               /* water table position relative to soil surface within the layer (cm) */
  double zwt_dist_zwt[MAX_NZWT];   /* array of layer water table positions when DIST_ZWT = TRUE (cm) */
} layer_data_struct;

/******************************************************************
  This structure stores soil variables for the complete soil column 
  for each grid cell.
  ******************************************************************/
typedef struct {
  double aero_resist[2];               /* The (stability-corrected) aerodynamic resistance (s/m) that was actually used in flux calculations.  [0] = surface (bare soil, non-overstory veg, or snow pack) [1] = overstory */
  double asat;                         /* saturated area fraction */
  double asat_dist_zwt[MAX_NZWT];      /* array of saturated area fractions when DIST_ZWT = TRUE */
  double baseflow;                     /* baseflow from current cell (mm/TS) */
  double baseflow_dist_zwt[MAX_NZWT]; /* array of baseflow values when DIST_ZWT = TRUE (mm) */
  double CLitter;                      /* carbon storage in litter pool [gC/m2] */
  double CInter;                       /* carbon storage in intermediate pool [gC/m2] */
  double CSlow;                        /* carbon storage in slow pool [gC/m2] */
  double inflow;                       /* moisture that reaches the top of 
					  the soil column (mm) */
  double pot_evap[N_PET_TYPES];        /* array of different types of potential evaporation (mm) */
  double runoff;                       /* runoff from current cell (mm/TS) */
  double runoff_dist_zwt[MAX_NZWT];   /* array of runoff values when DIST_ZWT = TRUE (mm) */
  layer_data_struct layer[MAX_LAYERS]; /* structure containing soil variables 
					  for each layer (see above) */
  double RhLitter;                     /* soil respiration from litter pool [gC/m2] */
  double RhLitter2Atm;                 /* soil respiration from litter pool [gC/m2] that goes to atmosphere */
  double RhInter;                      /* soil respiration from intermediate pool [gC/m2] */
  double RhSlow;                       /* soil respiration from slow pool [gC/m2] */
  double RhTot;                        /* total soil respiration over all pools [gC/m2] (=RhLitter2Atm+RhInter+RhSlow) */
  double rootmoist;                    /* total of layer.moist over all layers
                                          in the root zone (mm) */
  double wetness;                      /* average of
                                          (layer.moist - Wpwp)/(porosity*depth - Wpwp)
                                          over all layers (fraction) */
  double zwt;                          /* average water table position [cm] - using lowest unsaturated layer */
  double zwt_dist_zwt[MAX_NZWT];       /* array of water table positions (using lowest unsaturated layer) when DIST_ZWT = TRUE (cm) */
  double zwt_lumped;                   /* average water table position [cm] - lumping all layers' moisture together */
  double zwt_lumped_dist_zwt[MAX_NZWT];     /* array of water table positions (lumping all layers' moitures together) when DIST_ZWT = TRUE (cm) */
} cell_data_struct;

/***********************************************************************
  This structure stores energy balance components, and variables used to
  solve the thermal fluxes through the soil column.
  ***********************************************************************/
typedef struct {
  // State variables
  double  AlbedoLake;            /* albedo of lake surface (fract) */
  double  AlbedoOver;            /* albedo of intercepted snow (fract) */
  double  AlbedoUnder;           /* surface albedo (fraction) */
  double  Cs[2];                 /* heat capacity for top two layers (J/m^3/K) */
  double  Cs_node[MAX_NODES];    /* heat capacity of the soil thermal nodes (J/m^3/K) */
  double  fdepth[MAX_FRONTS];    /* all simulated freezing front depths */
  char    frozen;                /* TRUE = frozen soil present */
  double  ice[MAX_NODES];        /* thermal node ice content */
  double  kappa[2];              /* soil thermal conductivity for top two layers (W/m/K) */
  double  kappa_node[MAX_NODES]; /* thermal conductivity of the soil thermal nodes (W/m/K) */
  double  moist[MAX_NODES];      /* thermal node moisture content */
  int     Nfrost;                /* number of simulated freezing fronts */
  int     Nthaw;                 /* number of simulated thawing fronts */
  double  T[MAX_NODES];          /* thermal node temperatures (C) */
  char    T_fbflag[MAX_NODES];   /* flag indicating if previous step's temperature was used */
  int     T_fbcount[MAX_NODES];  /* running total number of times that previous step's temperature was used */
  int     T1_index;              /* soil node at the bottom of the top layer */
  double  Tcanopy;               /* temperature of the canopy air */
  char    Tcanopy_fbflag;        /* flag indicating if previous step's temperature was used */
  int     Tcanopy_fbcount;       /* running total number of times that previous step's temperature was used */
  double  tdepth[MAX_FRONTS];    /* all simulated thawing front depths */
  double  Tfoliage;              /* temperature of the overstory vegetation */
  char    Tfoliage_fbflag;       /* flag indicating if previous step's temperature was used */
  int     Tfoliage_fbcount;      /* running total number of times that previous step's temperature was used */
  double  Tlakebot;              /* temperature of the soil surface under the lake */
  char    Tlakebot_fbflag;       /* flag indicating if previous step's temperature was used */
  int     Tlakebot_fbcount;      /* running total number of times that previous step's temperature was used */
  double  Tsurf;                 /* temperature of the understory */
  char    Tsurf_fbflag;          /* flag indicating if previous step's temperature was used */
  int     Tsurf_fbcount;         /* running total number of times that previous step's temperature was used */
  double  unfrozen;              /* frozen layer water content that is unfrozen */
  // Fluxes
  double  advected_sensible;     /* net sensible heat flux advected to snowpack (Wm-2) */
  double  advection;             /* advective flux (Wm-2) */
  double  AtmosError;
  double  AtmosLatent;           /* latent heat exchange with atmosphere */
  double  AtmosLatentSub;        /* latent sub heat exchange with atmosphere */
  double  AtmosSensible;         /* sensible heat exchange with atmosphere */
  double  canopy_advection;      /* advection heat flux from the canopy (W/m^2) */
  double  canopy_latent;         /* latent heat flux from the canopy (W/m^2) */
  double  canopy_latent_sub;     /* latent heat flux of sublimation from the canopy (W/m^2) */
  double  canopy_refreeze;       /* energy used to refreeze/melt canopy intercepted snow (W/m^2) */
  double  canopy_sensible;       /* sensible heat flux from canopy interception (W/m^2) */
  double  deltaCC;               /* change in snow heat storage (Wm-2) */
  double  deltaH;                /* change in soil heat storage (Wm-2) */
  double  error;                 /* energy balance error (W/m^2) */
  double  fusion;                /* energy used to freeze/thaw soil water */
  double  grnd_flux;             /* ground heat flux (Wm-2) */
  double  latent;                /* net latent heat flux (Wm-2) */
  double  latent_sub;            /* net latent heat flux from snow (Wm-2) */
  double  longwave;              /* net longwave flux (Wm-2) */
  double  LongOverIn;            /* incoming longwave to overstory */
  double  LongUnderIn;           /* incoming longwave to understory */
  double  LongUnderOut;          /* outgoing longwave from understory */
  double  melt_energy;           /* energy used to reduce snow cover fraction (Wm-2) */
  double  NetLongAtmos;          /* net longwave radiation to the atmosphere (W/m^2) */
  double  NetLongOver;           /* net longwave radiation from the overstory (W/m^2) */
  double  NetLongUnder;          /* net longwave radiation from the understory (W/m^2) */
  double  NetShortAtmos;         /* net shortwave to the atmosphere */
  double  NetShortGrnd;          /* net shortwave penetrating snowpack */
  double  NetShortOver;          /* net shortwave radiation from the overstory (W/m^2) */
  double  NetShortUnder;         /* net shortwave radiation from the understory (W/m^2) */
  double  out_long_canopy;       /* outgoing longwave to canopy */
  double  out_long_surface;      /* outgoing longwave to surface */
  double  refreeze_energy;       /* energy used to refreeze the snowpack (Wm-2) */
  double  sensible;              /* net sensible heat flux (Wm-2) */
  double  shortwave;             /* net shortwave radiation (Wm-2) */
  double  ShortOverIn;           /* incoming shortwave to overstory */
  double  ShortUnderIn;          /* incoming shortwave to understory */
  double  snow_flux;             /* thermal flux through the snow pack (Wm-2) */
  double  lake_soil_heat_flux;   /* Heat flux from lake water into underlying soil (W/m^2) */
  double  lake_soil_net_short;   /* Net shortwave radiation from lake water into underlying soil (W/m^2) */
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
  double *NscaleFactor;         /* array of per-layer nitrogen scaling factors */
  double *aPARLayer;            /* array of per-layer absorbed PAR (mol(photons)/m2 leaf area s) */
  double *CiLayer;              /* array of per-layer leaf-internal CO2 mixing ratio (mol CO2/mol air) */
  double *rsLayer;              /* array of per-layer stomatal resistance (s/m) */
  double **rsLayer_dist_zwt;    /* 2-dimensional array of per-layer stomatal resistance (s/m) when DIST_ZWT = TRUE */
  double aPAR;                  /* whole-canopy absorbed PAR (mol(photons)/m2 leaf area s) */
  double Ci;                    /* whole-canopy leaf-internal CO2 mixing ratio (mol CO2/mol air) */
  double rc;                    /* whole-canopy stomatal resistance (s/m) */
  double *rc_dist_zwt;          /* array of whole-canopy stomatal resistances (s/m) when DIST_ZWT = TRUE */
  double NPPfactor;             /* whole-canopy photosynthesis multiplier to account for inhibition separate from stomatal resistance */
  double *NPPfactor_dist_zwt;   /* array of whole-canopy photosynthesis multipliers when DIST_ZWT = TRUE */
  double GPP;                   /* whole-canopy gross assimilation (photosynthesis) (umol(CO2)/m2s) */
  double Rphoto;                /* whole-canopy photorespiration (umol(CO2)/m2s) */
  double Rdark;                 /* whole-canopy 'dark' respiration (umol(CO2)/m2s) */
  double Rmaint;                /* plant maintenance respiration (= Rdark/FRLeaf) (umol(CO2)/m2s) */
  double Rgrowth;               /* growth respiration ( = (GPP-Rmaint)*FRGrowth/(1+FRGrowth) ) (umol(CO2)/m2s) */
  double Raut;                  /* total plant respiration (= Rmaint + Rgrowth) (umol(CO2)/m2s) */
  double NPP;                   /* net primary productivity (= GPP - Raut) (umol(CO2)/m2s) */
  double Litterfall;            /* flux of carbon from living biomass to litter pool [gC/m2] */
  double AnnualNPP;             /* running total annual NPP [gC/m2] */
  double AnnualNPPPrev;         /* total annual NPP from previous year [gC/m2] */
} veg_var_struct;

/************************************************************************
  This structure stores snow pack variables needed to run the snow model.
  ************************************************************************/
typedef struct {
  // State variables
  double albedo;            /* snow surface albedo (fraction) */
  double canopy_albedo;     /* albedo of the canopy (fract) */
  double coldcontent;       /* cold content of snow pack */
  double coverage;          /* fraction of snow band that is covered with snow */
  double density;           /* snow density (kg/m^3) */
  double depth;             /* snow depth (m) */
  int    last_snow;         /* time steps since last snowfall */
  double max_snow_depth;    /* last maximum snow depth - used to determine coverage
			       fraction during current melt period (m) */
  char   MELTING;           /* flag indicating that snowpack melted 
			       previously */
  double pack_temp;         /* depth averaged temperature of the snowpack (C) */
  double pack_water;        /* liquid water content of the snow pack (m) */
  int    snow;              /* TRUE = snow, FALSE = no snow */
  double snow_canopy;       /* amount of snow on canopy (m) */
  double store_coverage;    /* stores coverage fraction covered by new snow (m) */
  int    store_snow;        /* flag indicating whether or not new accumulation
			       is stored on top of an existing distribution */
  double store_swq;         /* stores newly accumulated snow over an 
			       established snowpack melt distribution (m) */
  double surf_temp;         /* depth averaged temperature of the snow pack surface layer (C) */
  double surf_temp_fbcount; /* running total number of times that previous step's temperature was used */
  double surf_temp_fbflag;  /* flag indicating if previous step's temperature was used */
  double surf_water;        /* liquid water content of the surface layer (m) */
  double swq;               /* snow water equivalent of the entire pack (m) */
  double snow_distrib_slope;/* current slope of uniform snow distribution (m/fract) */
  double tmp_int_storage;   /* temporary canopy storage, used in snow_canopy */
  // Fluxes
  double blowing_flux;      /* depth of sublimation from blowing snow (m) */
  double canopy_vapor_flux; /* depth of water evaporation, sublimation, or 
			       condensation from intercepted snow (m) */
  double mass_error;        /* snow mass balance error */
  double melt;              /* snowpack melt (mm) */
  double Qnet;              /* Residual of energy balance at snowpack surface */
  double surface_flux;      /* depth of sublimation from blowing snow (m) */
  double transport;	    /* flux of snow (potentially) transported from veg type */
  double vapor_flux;        /* depth of water evaporation, sublimation, or 
			       condensation from snow pack (m) */
} snow_data_struct;

/******************************************************************
  This structure stores the lake/wetland parameters for a grid cell
  ******************************************************************/
typedef struct {
  // Lake basin dimensions
  int    numnod;                  /* Maximum number of lake nodes for this grid cell */
  double z[MAX_LAKE_NODES+1];     /* Elevation of each lake node (when lake storage is at maximum), relative to lake's deepest point (m) */  
  double basin[MAX_LAKE_NODES+1]; /* Area of lake basin at each lake node (when lake storage is at maximum) (m^2) */
  double Cl[MAX_LAKE_NODES+1];    /* Fractional coverage of lake basin at each node (when lake storage is at maximum) (fraction of grid cell area) */
  double b;                       /* Exponent in default lake depth-area profile (y=Ax^b) */
  double maxdepth;                /* Maximum allowable depth of liquid portion of lake (m) */
  double mindepth;                /* Minimum allowable depth of liquid portion of lake (m) */
  double maxvolume;               /* Lake volume when lake depth is at maximum (m^3) */
  double minvolume;               /* Lake volume when lake depth is at minimum (m^3) */
  // Hydrological properties
  float  bpercent;                /* Fraction of wetland baseflow (subsurface runoff) that flows into lake */
  float  rpercent;                /* Fraction of wetland surface runoff that flows into lake */
  double eta_a;                   /* Decline of solar radiation w/ depth (m^-1) */ /* not currently used */
  double wfrac;                   /* Width of lake outlet, expressed as fraction of lake perimeter */
  // Initial conditions
  double depth_in;                /* Initial lake depth (distance from surface to deepest point) (m) */
  int    lake_idx;                /* index number of the lake/wetland veg tile */
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
  double surface[MAX_LAKE_NODES+1];/* Area of horizontal cross-section of liquid water in lake at each node (m^2) */
  double sarea;                   /* Current surface area of ice+liquid water on lake surface (m^2) */
  double sarea_save;              /* Surface area of ice+liquid water on lake surface (m^2) at beginning of time step */
  double volume;                  /* Current lake water volume, including liquid water equivalent of lake ice (m^3) */
  double volume_save;             /* Lake water volume, including liquid water equivalent of lake ice (m^3) at beginning of time step */
  double temp[MAX_LAKE_NODES];    /* Lake water temperature at each node (C) */
  double tempavg;                 /* Average liquid water temperature of entire lake (C) */
  // Current properties (state variables) specific to lake ice/snow
  double areai;                   /* Area of ice coverage (at beginning of time step) (m^2) */
  double new_ice_area;            /* Area of ice coverage (at end of time step) (m^2) */
  double ice_water_eq;            /* Liquid water equivalent volume of lake ice (m^3) */
  double hice;                    /* Height of lake ice at thickest point (m) */ 
  double tempi;                   /* Lake ice temperature (C) */
  double swe;                     /* Water equivalence of lake snow cover - end of step (m^3) */
  double swe_save;                /* Water equivalence of lake snow cover - beginning of step (m^3) */
  double surf_temp;               /* Temperature of surface snow layer (C) */
  double pack_temp;               /* Temperature of pack snow layer (C) */
  double coldcontent;             /* cold content of snow pack */
  double surf_water;              /* Water content of surface snow layer (m^3) */
  double pack_water;              /* Water content of pack snow layer (m^3) */
  double SAlbedo;                 /* Albedo of lake snow (fraction) */
  double sdepth;                  /* Depth of snow on top of ice (m^3) */
  // Other current lake properties (derived from state variables and forcings)
  double aero_resist;	          /* Aerodynamic resistance (s/m) after stability correction */
  double density[MAX_LAKE_NODES]; /* Lake water density profile (kg/m^3) */
  // Moisture fluxes
  double baseflow_in;             /* Volume of baseflow into lake from the rest of the grid cell (m3) */
  double baseflow_out;            /* Volume of baseflow out of lake to channel network (m3) */
  double channel_in;              /* Volume of channel inflow into lake (m3) */
  double evapw;                   /* Volume of evaporative flux from lake (and ice/snow) surface (m3) */
  double ice_throughfall;         /* Volume of precipitation reaching lake water surface, i.e. total precip minus accumulation of snow on ice (m3) */
  double prec;                    /* Volume of precipitation falling on lake (and ice/snow) surface (m3) */
  double recharge;                /* Volume of recharge from lake to wetland (m3) */
  double runoff_in;               /* Volume of surface runoff into lake from the rest of the grid cell (m3) */
  double runoff_out;              /* Volume of surface runoff out of lake to channel network (m3) */
  double snowmlt;                 /* Volume of moisture released by melting of lake snow (m3) */
  double vapor_flux;              /* Volume of moisture sublimated from lake snow (m3) */
  // Structures compatible with other land cover types
  // Some of this information is currently redundant with other variables in the lake_var structure
  snow_data_struct  snow;         /* Snow pack on top of lake ice */
  energy_bal_struct energy;       /* Energy fluxes and soil temperatures */
  cell_data_struct  soil;         /* Soil column below lake */
} lake_var_struct;

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

