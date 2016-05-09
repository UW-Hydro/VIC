/******************************************************************************
 * @section DESCRIPTION
 *
 * Header file for vic_driver_shared_all routines
 *
 * @section LICENSE
 *
 * The Variable Infiltration Capacity (VIC) macroscale hydrological model
 * Copyright (C) 2014 The Land Surface Hydrology Group, Department of Civil
 * and Environmental Engineering, University of Washington.
 *
 * The VIC model is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 *****************************************************************************/

#ifndef VIC_DRIVER_SHARED_H
#define VIC_DRIVER_SHARED_H

#include <vic_run.h>

#define VERSION "5.0 beta 2016 April 25"
#define SHORT_VERSION "5.0.beta"

/***** Define maximum array sizes for driver level objects *****/
#define MAX_FORCE_FILES 2

#define DEFAULT_COMPRESSION_LVL -1

/******************************************************************************
 * @brief   File formats
 *****************************************************************************/
enum
{
    UNSET_FILE_FORMAT,
    ASCII,
    BINARY,
    NETCDF3_CLASSIC,
    NETCDF3_64BIT_OFFSET,
    NETCDF4_CLASSIC,
    NETCDF4
};

/******************************************************************************
 * @brief   endian flags
 *****************************************************************************/
enum
{
    LITTLE,  /**< little-endian flag */
    BIG      /**< big-endian flag */
};

/******************************************************************************
 * @brief   Veg param sources
 *****************************************************************************/
enum
{
    FROM_DEFAULT,
    FROM_VEGLIB,
    FROM_VEGPARAM,
    FROM_VEGHIST
};

/******************************************************************************
 * @brief   Forcing Variable Types
 *****************************************************************************/
enum
{
    AIR_TEMP,    /**< air temperature per time step [C] */
    ALBEDO,      /**< surface albedo [fraction] */
    CATM,        /**< atmospheric CO2 concentration [ppm] */
    CHANNEL_IN,  /**< incoming channel flow [m3] */
    FCANOPY,     /**< fractional area covered by plant canopy [fraction] */
    FDIR,        /**< fraction of incoming shortwave that is direct [fraction] */
    LAI_IN,      /**< leaf area index [m2/m2] */
    LWDOWN,      /**< incoming longwave radiation [W/m2] */
    PAR,         /**< incoming photosynthetically active radiation [W/m2] */
    PREC,        /**< total precipitation (rain and snow) [mm] */
    PRESSURE,    /**< atmospheric pressure [kPa] */
    VP,          /**< vapor pressure [kPa] */
    SWDOWN,      /**< incoming shortwave [W/m2] */
    WIND,        /**< wind speed [m/s] */
    SKIP,        /**< place holder for unused data columns */
    // Last value of enum - DO NOT ADD ANYTHING BELOW THIS LINE!!
    // used as a loop counter and must be >= the largest value in this enum
    N_FORCING_TYPES  /**< Number of forcing types */
};


/******************************************************************************
 * @brief   Output Variable Types
 *****************************************************************************/
enum
{
    // Water Balance Terms - state variables
    OUT_ASAT,             /**< Saturated Area Fraction */
    OUT_LAKE_AREA_FRAC,   /**< lake surface area as fraction of the grid cell area [fraction] */
    OUT_LAKE_DEPTH,       /**< lake depth (distance between surface and deepest point) [m] */
    OUT_LAKE_ICE,         /**< moisture stored as lake ice [mm over lake ice area] */
    OUT_LAKE_ICE_FRACT,   /**< fractional coverage of lake ice [fraction] */
    OUT_LAKE_ICE_HEIGHT,  /**< thickness of lake ice [cm] */
    OUT_LAKE_MOIST,       /**< liquid water and ice stored in lake [mm over grid cell] */
    OUT_LAKE_SURF_AREA,   /**< lake surface area [m2] */
    OUT_LAKE_SWE,         /**< liquid water equivalent of snow on top of lake ice [m over lake ice area] */
    OUT_LAKE_SWE_V,       /**< volumetric liquid water equivalent of snow on top of lake ice [m3] */
    OUT_LAKE_VOLUME,      /**< lake volume [m3] */
    OUT_ROOTMOIST,        /**< root zone soil moisture  [mm] */
    OUT_SMFROZFRAC,       /**< fraction of soil moisture (by mass) that is ice, for each soil layer */
    OUT_SMLIQFRAC,        /**< fraction of soil moisture (by mass) that is liquid, for each soil layer */
    OUT_SNOW_CANOPY,      /**< snow interception storage in canopy  [mm] */
    OUT_SNOW_COVER,       /**< fractional area of snow cover [fraction] */
    OUT_SNOW_DEPTH,       /**< depth of snow pack [cm] */
    OUT_SOIL_ICE,         /**< soil ice content  [mm] for each soil layer */
    OUT_SOIL_LIQ,         /**< soil liquid content  [mm] for each soil layer */
    OUT_SOIL_ICE_FRAC,    /**< soil ice content fraction of column volume [1] for each soil layer */
    OUT_SOIL_LIQ_FRAC,    /**< soil liquid content fraction of column volume [1] for each soil layer */
    OUT_SOIL_MOIST,       /**< soil total moisture content  [mm] for each soil layer */
    OUT_SOIL_WET,         /**< vertical average of (soil moisture - wilting point)/(maximum soil moisture - wilting point) [mm/mm] */
    OUT_SURFSTOR,         /**< storage of liquid water and ice (not snow) on surface (ponding) [mm] */
    OUT_SURF_FROST_FRAC,  /**< fraction of soil surface that is frozen [fraction] */
    OUT_SWE,              /**< snow water equivalent in snow pack (including vegetation-intercepted snow)  [mm] */
    OUT_WDEW,             /**< total moisture interception storage in canopy [mm] */
    OUT_ZWT,              /**< water table position [cm] (zwt within lowest unsaturated layer) */
    OUT_ZWT_LUMPED,       /**< lumped water table position [cm] (zwt of total moisture across all layers, lumped together) */
    // Water Balance Terms - fluxes
    OUT_BASEFLOW,         /**< baseflow out of the bottom layer  [mm] */
    OUT_DELINTERCEPT,     /**< change in canopy interception storage  [mm] */
    OUT_DELSOILMOIST,     /**< change in soil water content  [mm] */
    OUT_DELSURFSTOR,      /**< change in surface liquid water storage  [mm] */
    OUT_DELSWE,           /**< change in snow water equivalent  [mm] */
    OUT_EVAP,             /**< total net evaporation [mm] */
    OUT_EVAP_BARE,        /**< net evaporation from bare soil [mm] */
    OUT_EVAP_CANOP,       /**< net evaporation from canopy interception [mm] */
    OUT_INFLOW,           /**< moisture that reaches top of soil column [mm] */
    OUT_LAKE_BF_IN,       /**< incoming baseflow from lake catchment [mm] */
    OUT_LAKE_BF_IN_V,     /**< incoming volumetric baseflow from lake catchment [m3] */
    OUT_LAKE_BF_OUT,      /**< outgoing baseflow from lake [mm] */
    OUT_LAKE_BF_OUT_V,    /**< outgoing volumetric baseflow from lake [m3] */
    OUT_LAKE_CHAN_IN,     /**< channel inflow into lake [mm] */
    OUT_LAKE_CHAN_IN_V,   /**< volumetric channel inflow into lake [m3] */
    OUT_LAKE_CHAN_OUT,    /**< channel outflow from lake [mm] */
    OUT_LAKE_CHAN_OUT_V,  /**< volumetric channel outflow from lake [m3] */
    OUT_LAKE_DSTOR,       /**< change in lake moisture storage (liquid plus ice cover) [mm] */
    OUT_LAKE_DSTOR_V,     /**< volumetric change in lake moisture storage (liquid plus ice cover) [m3] */
    OUT_LAKE_DSWE,        /**< change in swe on top of lake ice [mm] */
    OUT_LAKE_DSWE_V,      /**< volumetric change in swe on top of lake ice [m3] */
    OUT_LAKE_EVAP,        /**< net evaporation from lake surface [mm] */
    OUT_LAKE_EVAP_V,      /**< net volumetric evaporation from lake surface [m3] */
    OUT_LAKE_PREC_V,      /**< volumetric precipitation over lake surface [m3] */
    OUT_LAKE_RCHRG,       /**< recharge from lake to surrounding wetland [mm] */
    OUT_LAKE_RCHRG_V,     /**< volumetric recharge from lake to surrounding wetland [m3] */
    OUT_LAKE_RO_IN,       /**< incoming runoff from lake catchment [mm] */
    OUT_LAKE_RO_IN_V,     /**< incoming volumetric runoff from lake catchment [m3] */
    OUT_LAKE_VAPFLX,      /**< outgoing sublimation from snow on top of lake ice [mm] */
    OUT_LAKE_VAPFLX_V,    /**< outgoing volumetric sublimation from snow on top of lake ice [m3] */
    OUT_PET,              /**< Potential evapotranspiration (= area-weighted sum of potential transpiration and potential soil evaporation).  Potential transpiration is computed using the Penman-Monteith eqn with architectural resistance and LAI of the current veg cover. [mm] */
    OUT_PREC,             /**< incoming precipitation [mm] */
    OUT_RAINF,            /**< rainfall  [mm] */
    OUT_REFREEZE,         /**< refreezing of water in the snow  [mm] */
    OUT_RUNOFF,           /**< surface runoff [mm] */
    OUT_SNOW_MELT,        /**< snow melt  [mm] */
    OUT_SNOWF,            /**< snowfall  [mm] */
    OUT_SUB_BLOWING,      /**< net sublimation of blowing snow [mm] */
    OUT_SUB_CANOP,        /**< net sublimation from snow stored in canopy [mm] */
    OUT_SUB_SNOW,         /**< total net sublimation from snow pack (surface and blowing) [mm] */
    OUT_SUB_SURFACE,      /**< net sublimation from snow pack surface [mm] */
    OUT_TRANSP_VEG,       /**< net transpiration from vegetation [mm] */
    OUT_WATER_ERROR,      /**< water budget error [mm] */
    // Energy Balance Terms - state variables
    OUT_ALBEDO,           /**< average surface albedo [fraction] */
    OUT_BARESOILT,        /**< bare soil surface temperature [C] */
    OUT_FDEPTH,           /**< depth of freezing fronts [cm] */
    OUT_LAKE_ICE_TEMP,    /**< temperature of lake ice [C] */
    OUT_LAKE_SURF_TEMP,   /**< lake surface temperature [C] */
    OUT_RAD_TEMP,         /**< average radiative surface temperature [K] */
    OUT_SALBEDO,          /**< snow pack albedo [fraction] */
    OUT_SNOW_PACK_TEMP,   /**< snow pack temperature [C] */
    OUT_SNOW_SURF_TEMP,   /**< snow surface temperature [C] */
    OUT_SNOWT_FBFLAG,     /**< snow surface temperature fallback flag */
    OUT_SOIL_TEMP,        /**< soil temperature [C] */
    OUT_SOIL_TNODE,       /**< soil temperature [C] */
    OUT_SOIL_TNODE_WL,    /**< soil temperature [C] */
    OUT_SOILT_FBFLAG,     /**< soil temperature flag for each soil thermal node */
    OUT_SURF_TEMP,        /**< average surface temperature [C] */
    OUT_SURFT_FBFLAG,     /**< surface temperature flag */
    OUT_TCAN_FBFLAG,      /**< Tcanopy flag */
    OUT_TDEPTH,           /**< depth of thawing fronts [cm] */
    OUT_TFOL_FBFLAG,      /**< Tfoliage flag */
    OUT_VEGT,             /**< average vegetation canopy temperature [C] */
    // Energy Balance Terms - fluxes
    OUT_ADV_SENS,         /**< net sensible flux advected to snow pack [W/m2] */
    OUT_ADVECTION,        /**< advected energy [W/m2] */
    OUT_DELTACC,          /**< rate of change in cold content in snow pack [W/m2] */
    OUT_DELTAH,           /**< rate of change in heat storage [W/m2] */
    OUT_ENERGY_ERROR,     /**< energy budget error [W/m2] */
    OUT_FUSION,           /**< net energy used to melt/freeze soil moisture [W/m2] */
    OUT_GRND_FLUX,        /**< net heat flux into ground [W/m2] */
    OUT_IN_LONG,          /**< incoming longwave at ground surface (under veg) [W/m2] */
    OUT_LATENT,           /**< net upward latent heat flux [W/m2] */
    OUT_LATENT_SUB,       /**< net upward latent heat flux from sublimation [W/m2] */
    OUT_MELT_ENERGY,      /**< energy of fusion (melting) in snowpack [W/m2] */
    OUT_LWNET,            /**< net downward longwave flux [W/m2] */
    OUT_SWNET,            /**< net downward shortwave flux [W/m2] */
    OUT_R_NET,            /**< net downward radiation flux [W/m2] */
    OUT_RFRZ_ENERGY,      /**< net energy used to refreeze liquid water in snowpack [W/m2] */
    OUT_SENSIBLE,         /**< net upward sensible heat flux [W/m2] */
    OUT_SNOW_FLUX,        /**< energy flux through snow pack [W/m2] */
    // Miscellaneous Terms
    OUT_AERO_COND,        /**< "scene" aerodynamic conductance [m/s] (tiles with overstory contribute overstory conductance; others contribute surface conductance) */
    OUT_AERO_COND1,       /**< surface aerodynamic conductance [m/s] */
    OUT_AERO_COND2,       /**< overstory aerodynamic conductance [m/s] */
    OUT_AERO_RESIST,      /**< "scene"canopy aerodynamic resistance [s/m]  (tiles with overstory contribute overstory resistance; others contribute surface resistance)*/
    OUT_AERO_RESIST1,     /**< surface aerodynamic resistance [s/m] */
    OUT_AERO_RESIST2,     /**< overstory aerodynamic resistance [s/m] */
    OUT_AIR_TEMP,         /**< air temperature [C] */
    OUT_CATM,             /**< atmospheric CO2 concentrtaion [ppm]*/
    OUT_DENSITY,          /**< near-surface atmospheric density [kg/m3]*/
    OUT_FCANOPY,          /**< fractional area covered by plant canopy [fraction] */
    OUT_FDIR,             /**< fraction of incoming shortwave that is direct [fraction]*/
    OUT_LAI,              /**< leaf area index [m2/m2] */
    OUT_LWDOWN,           /**< incoming longwave [W/m2] */
    OUT_PAR,              /**< incoming photosynthetically active radiation [W/m2] */
    OUT_PRESSURE,         /**< near surface atmospheric pressure [kPa] */
    OUT_QAIR,             /**< specific humidity [kg/kg] */
    OUT_REL_HUMID,        /**< relative humidity [%]*/
    OUT_SWDOWN,           /**< incoming shortwave [W/m2] */
    OUT_SURF_COND,        /**< surface conductance [m/s] */
    OUT_VP,               /**< near surface vapor pressure [kPa] */
    OUT_VPD,              /**< near surface vapor pressure deficit [kPa] */
    OUT_WIND,             /**< near surface wind speed [m/s] */
    // Band-specific quantities
    OUT_ADV_SENS_BAND,    /**< net sensible heat flux advected to snow pack [W/m2] */
    OUT_ADVECTION_BAND,   /**< advected energy [W/m2] */
    OUT_ALBEDO_BAND,      /**< average surface albedo [fraction] */
    OUT_DELTACC_BAND,     /**< change in cold content in snow pack [W/m2] */
    OUT_GRND_FLUX_BAND,   /**< net heat flux into ground [W/m2] */
    OUT_IN_LONG_BAND,     /**< incoming longwave at ground surface (under veg) [W/m2] */
    OUT_LATENT_BAND,      /**< net upward latent heat flux [W/m2] */
    OUT_LATENT_SUB_BAND,  /**< net upward latent heat flux due to sublimation [W/m2] */
    OUT_MELT_ENERGY_BAND, /**< energy of fusion (melting) in snowpack [W/m2] */
    OUT_LWNET_BAND,       /**< net downward longwave flux [W/m2] */
    OUT_SWNET_BAND,       /**< net downward shortwave flux [W/m2] */
    OUT_RFRZ_ENERGY_BAND, /**< net energy used to refreeze liquid water in snowpack [W/m2] */
    OUT_SENSIBLE_BAND,    /**< net upward sensible heat flux [W/m2] */
    OUT_SNOW_CANOPY_BAND, /**< snow interception storage in canopy [mm] */
    OUT_SNOW_COVER_BAND,  /**< fractional area of snow cover [fraction] */
    OUT_SNOW_DEPTH_BAND,  /**< depth of snow pack [cm] */
    OUT_SNOW_FLUX_BAND,   /**< energy flux through snow pack [W/m2] */
    OUT_SNOW_MELT_BAND,   /**< snow melt [mm] */
    OUT_SNOW_PACKT_BAND,  /**< snow pack temperature [C] */
    OUT_SNOW_SURFT_BAND,  /**< snow surface temperature [C] */
    OUT_SWE_BAND,         /**< snow water equivalent in snow pack [mm] */
    // Carbon-Cycling Terms
    OUT_APAR,             /**< absorbed PAR [W/m2] */
    OUT_GPP,              /**< gross primary productivity [g C/m2d] */
    OUT_RAUT,             /**< autotrophic respiration [g C/m2d] */
    OUT_NPP,              /**< net primary productivity [g C/m2d] */
    OUT_LITTERFALL,       /**< flux of carbon from living biomass into soil [g C/m2d] */
    OUT_RHET,             /**< soil respiration (heterotrophic respiration) [g C/m2d] */
    OUT_NEE,              /**< net ecosystem exchange (=NPP-RHET) [g C/m2d] */
    OUT_CLITTER,          /**< Carbon density in litter pool [g C/m2] */
    OUT_CINTER,           /**< Carbon density in intermediate pool [g C/m2] */
    OUT_CSLOW,            /**< Carbon density in slow pool [g C/m2] */
    // Last value of enum - DO NOT ADD ANYTHING BELOW THIS LINE!!
    // used as a loop counter and must be >= the largest value in this enum
    N_OUTVAR_TYPES        /**< used as a loop counter*/
};


/******************************************************************************
 * @brief   Output BINARY format types
 *****************************************************************************/
enum
{
    OUT_TYPE_DEFAULT,  /**< Default data type */
    OUT_TYPE_CHAR,     /**< char */
    OUT_TYPE_SINT,     /**< short int */
    OUT_TYPE_USINT,    /**< unsigned short int */
    OUT_TYPE_INT,      /**< int */
    OUT_TYPE_FLOAT,    /**< single-precision floating point */
    OUT_TYPE_DOUBLE    /**< double-precision floating point */
};

/******************************************************************************
 * @brief   Output aggregation method types
 *****************************************************************************/
enum
{
    AGG_TYPE_AVG,     /**< average over agg interval */
    AGG_TYPE_BEG,     /**< value at beginning of agg interval */
    AGG_TYPE_END,     /**< value at end of agg interval */
    AGG_TYPE_MAX,     /**< maximum value over agg interval */
    AGG_TYPE_MIN,     /**< minimum value over agg interval */
    AGG_TYPE_SUM     /**< sum over agg interval */
};


/******************************************************************************
 * @brief   Codes for displaying version information
 *****************************************************************************/
enum
{
    DISP_VERSION,
    DISP_COMPILE_TIME,
    DISP_ALL
};

/******************************************************************************
 * @brief   Codes for calendar option.
 *****************************************************************************/
enum calendars
{
    CALENDAR_STANDARD,
    CALENDAR_GREGORIAN,
    CALENDAR_PROLEPTIC_GREGORIAN,
    CALENDAR_NOLEAP,
    CALENDAR_365_DAY,
    CALENDAR_360_DAY,
    CALENDAR_JULIAN,
    CALENDAR_ALL_LEAP,
    CALENDAR_366_DAY
};

/******************************************************************************
 * @brief   Codes for time units option.
 *****************************************************************************/
enum time_units
{
    TIME_UNITS_SECONDS,
    TIME_UNITS_MINUTES,
    TIME_UNITS_HOURS,
    TIME_UNITS_DAYS
};

/******************************************************************************
 * @brief   file structures
 *****************************************************************************/
typedef struct {
    FILE *forcing[MAX_FORCE_FILES];   /**< atmospheric forcing data files */
    FILE *globalparam;  /**< global parameters file */
    FILE *constants;    /**< model constants parameter file */
    FILE *domain;       /**< domain file */
    FILE *init_state;   /**< initial model state file */
    FILE *lakeparam;    /**< lake parameter file */
    FILE *snowband;     /**< snow elevation band data file */
    FILE *soilparam;    /**< soil parameters for all grid cells */
    FILE *statefile;    /**< output model state file */
    FILE *veglib;       /**< vegetation parameters for all vege types */
    FILE *vegparam;     /**< fractional coverage info for grid cell */
    FILE *logfile;      /**< log file */
} filep_struct;

/******************************************************************************
 * @brief   This structure stores input and output filenames.
 *****************************************************************************/
typedef struct {
    char forcing[MAX_FORCE_FILES][MAXSTRING];    /**< atmospheric forcing data file names */
    char f_path_pfx[MAX_FORCE_FILES][MAXSTRING]; /**< path and prefix for atmospheric forcing data file names */
    char global[MAXSTRING];        /**< global control file name */
    char domain[MAXSTRING];        /**< domain file name */
    char constants[MAXSTRING];     /**< model constants file name */
    char init_state[MAXSTRING];    /**< initial model state file name */
    char lakeparam[MAXSTRING];     /**< lake model constants file */
    char result_dir[MAXSTRING];    /**< directory where results will be written */
    char snowband[MAXSTRING];      /**< snow band parameter file name */
    char soil[MAXSTRING];          /**< soil parameter file name */
    char statefile[MAXSTRING];     /**< name of file in which to store model state */
    char veg[MAXSTRING];           /**< vegetation grid coverage file */
    char veglib[MAXSTRING];        /**< vegetation parameter library file */
    char log_path[MAXSTRING];      /**< Location to write log file to*/
} filenames_struct;

/******************************************************************************
 * @brief    Stores forcing file input information.
 *****************************************************************************/
typedef struct {
    size_t N_ELEM; /**< number of elements per record; for LAI and ALBEDO,
                        1 element per veg tile; for others N_ELEM = 1; */
    bool SIGNED;
    bool SUPPLIED;
    double multiplier;
    char varname[MAXSTRING];
} force_type_struct;

/******************************************************************************
 * @brief    This structure records the parameters set by the forcing file
             input routines.  Those filled, are used to estimate the paramters
             needed for the model run in initialize_atmos.c.
 *****************************************************************************/
typedef struct {
    force_type_struct TYPE[N_FORCING_TYPES];
    double FORCE_DT[2];    /**< forcing file time step */
    size_t force_steps_per_day[2];    /**< forcing file timesteps per day */
    unsigned short int FORCE_ENDIAN[2];  /**< endian-ness of input file, used for
                                            DAILY_BINARY format */
    int FORCE_FORMAT[2];            /**< ASCII or BINARY */
    int FORCE_INDEX[2][N_FORCING_TYPES];
    size_t N_TYPES[2];
} param_set_struct;

/******************************************************************************
 * @brief   This structure stores output information for one output file.
 *****************************************************************************/
typedef struct {
    char prefix[MAXSTRING];                   /**< prefix of the file name, e.g. "fluxes" */
    char filename[MAXSTRING];                 /**< complete file name */
    FILE *fh;                               /**< filehandle */
    unsigned short int file_format;    /**< output file format */
    short int compress;            /**< Compress all output files */
    size_t output_steps_per_day;                /**< Number of output timesteps per day */
    double out_dt;                /**< output timestep in seconds */
    unsigned short int skipyear;    /**< Number of years to skip before writing
                                       output data */
    unsigned short int *type;   /**< type, when written to a binary file;
                                   OUT_TYPE_USINT  = unsigned short int
                                   OUT_TYPE_SINT   = short int
                                   OUT_TYPE_FLOAT  = single precision floating point
                                   OUT_TYPE_DOUBLE = double precision floating point */
    double *mult;               /**< multiplier, when written to a binary file [shape=(nvars, )] */
    char **format;                /**< format, when written to disk [shape=(nvars, )] */
} stream_file_struct;

/******************************************************************************
 * @brief   Thiss structure stores output information for one output stream.
 *****************************************************************************/
typedef struct {
    size_t nvars;            /**< number of variables to store in the file */
    unsigned int counter;
    unsigned int nextagg;
    unsigned int *varid;         /**< id numbers of the variables to store in the file
                                    (a variable's id number is its index in the out_data array).
                                    The order of the id numbers in the varid array
                                    is the order in which the variables will be written. */
    unsigned short int *aggtype;   /**< type of aggregation to use [shape=(nvars, )] */
    double ***aggdata;             /**< array of aggregated data values [shape=(nvars, nelem, nsteps)] */
} stream_struct;

/******************************************************************************
 * @brief   This structure stores moisture state information for differencing
 *          with next time step.
 *****************************************************************************/
typedef struct {
    double total_moist_storage;   /**< total moisture storage [mm] */
    double total_soil_moist;      /**< total column soil moisture [mm] */
    double surfstor;              /**< surface water storage [mm] */
    double swe;                   /**< snow water equivalent [mm] */
    double wdew;                  /**< canopy interception [mm] */
} save_data_struct;

/******************************************************************************
 * @brief   This structure stores output information for individual history streams
 *****************************************************************************/
typedef struct {
    char varname[MAXSTRING];  /**< name of variable */
    char long_name[MAXSTRING];  /**< cf long_name of variable */
    char units[MAXSTRING];  /**< units of variable */
    char description[MAXSTRING];  /**< descripition of variable */
    size_t nelem;          /**< number of data values */
} out_metadata_struct;

/******************************************************************************
 * @brief   This structure holds all variables needed for the error handling
 *          routines.
 *****************************************************************************/
typedef struct {
    atmos_data_struct *atmos;
    double dt;
    energy_bal_struct *energy;
    filep_struct filep;
    size_t rec;
    double **out_data;
    stream_struct *output_streams;
    snow_data_struct *snow;
    soil_con_struct soil_con;
    veg_con_struct *veg_con;
    veg_var_struct *veg_var;
} Error_struct;

double air_density(double t, double p);
void agg_stream_data(stream_struct *stream, double **out_data);
double all_30_day_from_dmy(dmy_struct *dmy);
double all_leap_from_dmy(dmy_struct *dmy);
double average(double *ar, size_t n);
double calc_energy_balance_error(double, double, double, double, double);
void calc_root_fractions(veg_con_struct *veg_con, soil_con_struct *soil_con);
double calc_water_balance_error(double, double, double, double);
unsigned short int calendar_from_chars(char *cal_chars);
bool check_write_flag(int rec);
void collect_eb_terms(energy_bal_struct, snow_data_struct, cell_data_struct,
                      double, double, double, bool, bool, double, bool, int,
                      double *, double, double **);
void collect_wb_terms(cell_data_struct, veg_var_struct, snow_data_struct,
                      double, double, double, bool, double, bool, double *,
                      double *, double **);
void compute_derived_state_vars(all_vars_struct *, soil_con_struct *,
                                veg_con_struct *);
void compute_lake_params(lake_con_struct *, soil_con_struct);
void compute_treeline(atmos_data_struct *, dmy_struct *, double, double *,
                      bool *);
size_t count_force_vars(FILE *gp);
void cmd_proc(int argc, char **argv, char *globalfilename);
void compress_files(char string[], short int level);
double **create_outdata(void);
stream_struct create_outstream(stream_struct *output_streams);
void get_current_datetime(char *cdt);
double date2num(double origin, dmy_struct *date, double tzoffset,
                unsigned short int calendar, unsigned short int time_units);
void dmy_all_30_day(double julian, dmy_struct *dmy);
void dmy_all_leap(double julian, dmy_struct *dmy);
void dmy_julian_day(double julian, unsigned short int calendar,
                    dmy_struct *dmy);
void dmy_no_leap_day(double julian, dmy_struct *dmy);
void dt_seconds_to_time_units(unsigned short int time_units, double dt_seconds,
                              double *dt_time_units);
void display_current_settings(int);
double fractional_day_from_dmy(dmy_struct *dmy);
void free_all_vars(all_vars_struct *all_vars, int Nveg);
void free_dmy(dmy_struct **dmy);
void free_out_data(double **);
void free_out_data_files(stream_file_struct **out_data_files);
void free_out_data_streams(stream_struct **streams);
void free_vegcon(veg_con_struct **veg_con);
void generate_default_state(all_vars_struct *, soil_con_struct *,
                            veg_con_struct *);
void generate_default_lake_state(all_vars_struct *, soil_con_struct *,
                                 lake_con_struct);
void get_parameters(FILE *paramfile);
void init_output_list(double **out_data, int write, char *format, int type,
                      double mult);
void initialize_energy(energy_bal_struct **energy, size_t nveg);
void initialize_filenames(void);
void initialize_fileps(void);
void initialize_global(void);
void initialize_options(void);
void initialize_parameters(void);
void initialize_put_data();
void initialize_snow(snow_data_struct **snow, size_t veg_num);
void initialize_soil(cell_data_struct **cell, size_t veg_num);
void initialize_time(void);
void initialize_veg(veg_var_struct **veg_var, size_t nveg);
double julian_day_from_dmy(dmy_struct *dmy, unsigned short int calendar);
bool leap_year(unsigned short int year, unsigned short int calendar);
all_vars_struct make_all_vars(size_t nveg);
cell_data_struct **make_cell_data(size_t veg_type_num);
dmy_struct *make_dmy(global_param_struct *global);
energy_bal_struct **make_energy_bal(size_t nveg);
void make_lastday(unsigned short int calendar, unsigned short int year,
                  unsigned short int lastday[]);
snow_data_struct **make_snow_data(size_t nveg);
veg_var_struct **make_veg_var(size_t veg_type_num);
double no_leap_day_from_dmy(dmy_struct *dmy);
void num2date(double origin, double time_value, double tzoffset,
              unsigned short int calendar, unsigned short int time_units,
              dmy_struct *date);
FILE *open_file(char string[], char type[]);
void parse_nc_time_units(char *nc_unit_chars, unsigned short int *units,
                         dmy_struct *dmy);
void put_data(all_vars_struct *, atmos_data_struct *, soil_con_struct *,
              veg_con_struct *, veg_lib_struct *veg_lib, lake_con_struct *,
              double **out_data, save_data_struct *);
void print_cell_data(cell_data_struct *cell, size_t nlayers, size_t nfrost);
void print_dmy(dmy_struct *dmy);
void print_energy_bal(energy_bal_struct *eb, size_t nnodes, size_t nfronts);
void print_filenames(filenames_struct *fnames);
void print_filep(filep_struct *fp);
void print_force_type(force_type_struct *force_type);
void print_global_param(global_param_struct *gp);
void print_lake_con(lake_con_struct *lcon, size_t nlnodes);
void print_lake_var(lake_var_struct *lvar, size_t nlnodes, size_t nfronts,
                    size_t nlayers, size_t nnodes, size_t nfrost);
void print_layer_data_states(layer_data_struct *ldata, size_t nfrost);
void print_layer_data_fluxes(layer_data_struct *ldata);
void print_license(void);
void print_option(option_struct *option);
void print_out_data(double **out_data, out_metadata_struct *metadata);
void print_out_metadata(out_metadata_struct *metadata);
void print_output_streams(stream_struct *outf);
void print_param_set(param_set_struct *param_set);
void print_parameters(parameters_struct *param);
void print_save_data(save_data_struct *save);
void print_snow_data(snow_data_struct *snow);
void print_soil_con(soil_con_struct *scon, size_t nlayers, size_t nnodes,
                    size_t nfrost, size_t nbands, size_t nzwt);
void print_stream_file(stream_file_struct *sf, stream_struct *stream,
                       out_metadata_struct *metadata);
void print_veg_con(veg_con_struct *vcon, size_t nroots, char blowing, char lake,
                   char carbon, size_t ncanopy);
void print_veg_lib(veg_lib_struct *vlib, char carbon);
void print_veg_var(veg_var_struct *vvar, size_t ncanopy);
void print_version(char *);
void print_usage(char *);
double q_to_vp(double q, double p);
void reset_stream(stream_struct *stream);
void set_output_var(stream_struct *stream, stream_file_struct *out_data_file,
                    char *varname, size_t varnum, char *format,
                    unsigned short int type, double mult);
unsigned int get_default_outvar_aggtype(unsigned int varid);
void set_output_met_data_info();
void setup_stream(stream_struct *stream, stream_file_struct *stream_file,
                  size_t nvars, unsigned int nextagg);
void soil_moisture_from_water_table(soil_con_struct *soil_con, size_t nlayers);
void sprint_dmy(char *str, dmy_struct *dmy);
bool strtotrueorfalse(char str[]);
unsigned short int timeunits_from_chars(char *units_chars);
int update_step_vars(all_vars_struct *, veg_con_struct *, veg_hist_struct *);
int invalid_date(unsigned short int calendar, dmy_struct *dmy);
void validate_parameters(void);
void validate_stream_settings(stream_file_struct **out_data_file);
char will_it_snow(double *t, double t_offset, double max_snow_temp,
                  double *prcp, size_t n);
void zero_output_list(double **);

#endif
