/******************************************************************************
 * @section DESCRIPTION
 *
 * Definition header file
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

#ifndef VIC_DEF_H
#define VIC_DEF_H

#include <float.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

#include <vic_physical_constants.h>
#include <vic_log.h>

/***** Model Constants *****/
#define MAXSTRING    2048
#define MISSING      -99999.   /**< missing value */
#define NODATA_VH    -1        /**< missing value for veg_hist inputs */
#define ERROR        -999      /**< Error Flag returned by subroutines */

/***** Define maximum array sizes for model source code *****/
#define MAX_VEG         12     /**< maximum number of vegetation types per cell */
#define MAX_LAYERS      3      /**< maximum number of soil moisture layers */
#define MAX_NODES       50     /**< maximum number of soil thermal nodes */
#define MAX_BANDS       10     /**< maximum number of snow bands */
#define MAX_FRONTS      3      /**< maximum number of freezing and thawing front depths to store */
#define MAX_FROST_AREAS 10     /**< maximum number of frost sub-areas */
#define MAX_LAKE_NODES  20     /**< maximum number of lake thermal nodes */
#define MAX_ZWTVMOIST   11     /**< maximum number of points in water table vs moisture curve for each soil layer; should include points at lower and upper boundaries of the layer */

/***** Define minimum values for model parameters *****/
#define MINSOILDEPTH    0.001  /**< Minimum layer depth with which model can work (m) */
#define MIN_VEGCOVER    0.0001 /**< Minimum allowable vegcover fraction */

/***** Potential Evap types *****/
#define N_PET_TYPES 0
#define N_PET_TYPES_NON_NAT 0
#define PET_SATSOIL 0
#define PET_H2OSURF 1
#define PET_SHORT   2
#define PET_TALL    3
#define N_PET_TYPES_NAT 0
#define PET_NATVEG  4
#define PET_VEGNOCR 5

/***** Hard-coded veg class parameters (mainly for pot_evap) *****/
extern bool   ref_veg_over[];
extern double ref_veg_rarc[];
extern double ref_veg_rmin[];
extern double ref_veg_lai[];
extern double ref_veg_albedo[];
extern double ref_veg_vegcover[];
extern double ref_veg_rough[];
extern double ref_veg_displ[];
extern double ref_veg_wind_h[];
extern double ref_veg_RGL[];
extern double ref_veg_rad_atten[];
extern double ref_veg_wind_atten[];
extern double ref_veg_trunk_ratio[];
extern bool   ref_veg_ref_crop[];

#ifndef WET
#define WET 0
#define DRY 1
#endif

#ifndef SNOW
#define RAIN 0
#define SNOW 1
#endif

#define min(a, b) (a < b) ? a : b
#define max(a, b) (a > b) ? a : b

extern size_t NR;       /**< array index for atmos struct that indicates
                             the model step avarage or sum */
extern size_t NF;       /**< array index loop counter limit for atmos
                             struct that indicates the SNOW_STEP values */

/******************************************************************************
 * @brief   Met file formats
 *****************************************************************************/
enum
{
    ASCII,
    BINARY
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
 * @brief   Snow Density parametrizations
 *****************************************************************************/
enum
{
    DENS_BRAS,
    DENS_SNTHRM
};

/******************************************************************************
 * @brief   Baseflow parametrizations
 *****************************************************************************/
enum
{
    ARNO,
    NIJSSEN2001
};

/******************************************************************************
 * @brief   Aerodynamic Resistance options
 *****************************************************************************/
enum
{
    AR_406,
    AR_406_LS,
    AR_406_FULL,
    AR_410
};

/******************************************************************************
 * @brief   Ground Flux options
 *****************************************************************************/
enum
{
    GF_406,
    GF_410
};

/******************************************************************************
 * @brief   VP algorithm options
 *****************************************************************************/
enum
{
    VP_ITER_NONE,
    VP_ITER_ALWAYS,
    VP_ITER_ANNUAL,
    VP_ITER_CONVERGE
};

/******************************************************************************
 * @brief   Longwave Clear-Sky Algorithm options
 *****************************************************************************/
enum
{
    LW_TVA,
    LW_ANDERSON,
    LW_BRUTSAERT,
    LW_SATTERLUND,
    LW_IDSO,
    LW_PRATA
};

/******************************************************************************
 * @brief   Longwave Cloud Algorithm options
 *****************************************************************************/
enum
{
    LW_CLOUD_BRAS,
    LW_CLOUD_DEARDORFF
};

/******************************************************************************
 * @brief   Veg param sources
 *****************************************************************************/
enum
{
    FROM_VEGLIB,
    FROM_VEGPARAM
};

/******************************************************************************
 * @brief   LAI sources (for backwards compatibility)
 *****************************************************************************/
enum
{
    LAI_FROM_VEGLIB,
    LAI_FROM_VEGPARAM
};

/******************************************************************************
 * @brief   Canopy resistance parametrizations
 *****************************************************************************/
enum
{
    RC_JARVIS,
    RC_PHOTO
};

/******************************************************************************
 * @brief   Photosynthesis parametrizations
 *****************************************************************************/
enum
{
    PS_FARQUHAR,
    PS_MONTEITH
};

/******************************************************************************
 * @brief   Photosynthetic pathways
 *****************************************************************************/
enum
{
    PHOTO_C3,
    PHOTO_C4
};

/******************************************************************************
 * @brief   Forcing Variable Types
 *****************************************************************************/
enum
{
    AIR_TEMP,    /**< air temperature per time step [C] (ALMA_INPUT: [K]) */
    ALBEDO,      /**< surface albedo [fraction] */
    CATM,        /**< atmospheric CO2 concentration [ppm] */
    CHANNEL_IN,  /**< incoming channel flow [m3] (ALMA_INPUT: [m3/s]) */
    CRAINF,      /**< convective rainfall [mm] (ALMA_INPUT: [mm/s]) */
    CSNOWF,      /**< convective snowfall [mm] (ALMA_INPUT: [mm/s]) */
    DENSITY,     /**< atmospheric density [kg/m3] */
    FDIR,        /**< fraction of incoming shortwave that is direct [fraction] */
    LAI_IN,      /**< leaf area index [m2/m2] */
    LONGWAVE,    /**< incoming longwave radiation [W/m2] */
    LSRAINF,     /**< large-scale rainfall [mm] (ALMA_INPUT: [mm/s]) */
    LSSNOWF,     /**< large-scale snowfall [mm] (ALMA_INPUT: [mm/s]) */
    PAR,         /**< incoming photosynthetically active radiation [W/m2] */
    PREC,        /**< total precipitation (rain and snow) [mm] (ALMA_INPUT: [mm/s]) */
    PRESSURE,    /**< atmospheric pressure [kPa] (ALMA_INPUT: [Pa]) */
    QAIR,        /**< specific humidity [kg/kg] */
    RAINF,       /**< rainfall (convective and large-scale) [mm] (ALMA_INPUT: [mm/s]) */
    REL_HUMID,   /**< relative humidity [%] */
    SHORTWAVE,   /**< incoming shortwave [W/m2] */
    SNOWF,       /**< snowfall (convective and large-scale) [mm] (ALMA_INPUT: [mm/s]) */
    TMAX,        /**< maximum daily temperature [C] (ALMA_INPUT: [K]) */
    TMIN,        /**< minimum daily temperature [C] (ALMA_INPUT: [K]) */
    TSKC,        /**< cloud cover fraction [fraction] */
    VEGCOVER,    /**< fraction of each veg tile covered by plants [fraction] */
    VP,          /**< vapor pressure [kPa] (ALMA_INPUT: [Pa]) */
    WIND,        /**< wind speed [m/s] */
    WIND_E,      /**< zonal component of wind speed [m/s] */
    WIND_N,      /**< meridional component of wind speed [m/s] */
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
    OUT_SOIL_MOIST,       /**< soil total moisture content  [mm] for each soil layer */
    OUT_SOIL_WET,         /**< vertical average of (soil moisture - wilting point)/(maximum soil moisture - wilting point) [mm/mm] */
    OUT_SURFSTOR,         /**< storage of liquid water and ice (not snow) on surface (ponding) [mm] */
    OUT_SURF_FROST_FRAC,  /**< fraction of soil surface that is frozen [fraction] */
    OUT_SWE,              /**< snow water equivalent in snow pack (including vegetation-intercepted snow)  [mm] */
    OUT_WDEW,             /**< total moisture interception storage in canopy [mm] */
    OUT_ZWT,              /**< water table position [cm] (zwt within lowest unsaturated layer) */
    OUT_ZWT_LUMPED,       /**< lumped water table position [cm] (zwt of total moisture across all layers, lumped together) */
    // Water Balance Terms - fluxes
    OUT_BASEFLOW,         /**< baseflow out of the bottom layer  [mm] (ALMA_OUTPUT: [mm/s]) */
    OUT_DELINTERCEPT,     /**< change in canopy interception storage  [mm] */
    OUT_DELSOILMOIST,     /**< change in soil water content  [mm] */
    OUT_DELSURFSTOR,      /**< change in surface liquid water storage  [mm] */
    OUT_DELSWE,           /**< change in snow water equivalent  [mm] */
    OUT_EVAP,             /**< total net evaporation [mm] (ALMA_OUTPUT: [mm/s]) */
    OUT_EVAP_BARE,        /**< net evaporation from bare soil [mm] (ALMA_OUTPUT: [mm/s]) */
    OUT_EVAP_CANOP,       /**< net evaporation from canopy interception [mm] (ALMA_OUTPUT: [mm/s]) */
    OUT_INFLOW,           /**< moisture that reaches top of soil column [mm] (ALMA_OUTPUT: [mm/s]) */
    OUT_LAKE_BF_IN,       /**< incoming baseflow from lake catchment [mm] (ALMA_OUTPUT: [mm/s]) */
    OUT_LAKE_BF_IN_V,     /**< incoming volumetric baseflow from lake catchment [m3] (ALMA_OUTPUT: [m3/s]) */
    OUT_LAKE_BF_OUT,      /**< outgoing baseflow from lake [mm] (ALMA_OUTPUT: [mm/s]) */
    OUT_LAKE_BF_OUT_V,    /**< outgoing volumetric baseflow from lake [m3] (ALMA_OUTPUT: [m3/s]) */
    OUT_LAKE_CHAN_IN,     /**< channel inflow into lake [mm] (ALMA_OUTPUT: [mm/s]) */
    OUT_LAKE_CHAN_IN_V,   /**< volumetric channel inflow into lake [m3] (ALMA_OUTPUT: [m3/s]) */
    OUT_LAKE_CHAN_OUT,    /**< channel outflow from lake [mm] (ALMA_OUTPUT: [mm/s]) */
    OUT_LAKE_CHAN_OUT_V,  /**< volumetric channel outflow from lake [m3] (ALMA_OUTPUT: [m3/s]) */
    OUT_LAKE_DSTOR,       /**< change in lake moisture storage (liquid plus ice cover) [mm] (ALMA_OUTPUT: [mm/s]) */
    OUT_LAKE_DSTOR_V,     /**< volumetric change in lake moisture storage (liquid plus ice cover) [m3] (ALMA_OUTPUT: [m3/s]) */
    OUT_LAKE_DSWE,        /**< change in swe on top of lake ice [mm] (ALMA_OUTPUT: [mm/s]) */
    OUT_LAKE_DSWE_V,      /**< volumetric change in swe on top of lake ice [m3] (ALMA_OUTPUT: [m3/s]) */
    OUT_LAKE_EVAP,        /**< net evaporation from lake surface [mm] (ALMA_OUTPUT: [mm/s]) */
    OUT_LAKE_EVAP_V,      /**< net volumetric evaporation from lake surface [m3] (ALMA_OUTPUT: [m3/s]) */
    OUT_LAKE_PREC_V,      /**< volumetric precipitation over lake surface [m3] (ALMA_OUTPUT: [m3/s]) */
    OUT_LAKE_RCHRG,       /**< recharge from lake to surrounding wetland [mm] (ALMA_OUTPUT: [mm/s]) */
    OUT_LAKE_RCHRG_V,     /**< volumetric recharge from lake to surrounding wetland [m3] (ALMA_OUTPUT: [m3/s]) */
    OUT_LAKE_RO_IN,       /**< incoming runoff from lake catchment [mm] (ALMA_OUTPUT: [mm/s]) */
    OUT_LAKE_RO_IN_V,     /**< incoming volumetric runoff from lake catchment [m3] (ALMA_OUTPUT: [m3/s]) */
    OUT_LAKE_VAPFLX,      /**< outgoing sublimation from snow on top of lake ice [mm] (ALMA_OUTPUT: [mm/s]) */
    OUT_LAKE_VAPFLX_V,    /**< outgoing volumetric sublimation from snow on top of lake ice [m3] (ALMA_OUTPUT: [m3/s]) */
    OUT_PET_SATSOIL,      /**< potential evap from saturated bare soil [mm] (ALMA_OUTPUT: [mm/s]) */
    OUT_PET_H2OSURF,      /**< potential evap from open water [mm] (ALMA_OUTPUT: [mm/s]) */
    OUT_PET_SHORT,        /**< potential evap (transpiration only) from short reference crop (grass) [mm] (ALMA_OUTPUT: [mm/s]) */
    OUT_PET_TALL,         /**< potential evap (transpiration only) from tall reference crop (alfalfa) [mm] (ALMA_OUTPUT: [mm/s]) */
    OUT_PET_NATVEG,       /**< potential evap (transpiration only) from current vegetation and current canopy resistance [mm] (ALMA_OUTPUT: [mm/s]) */
    OUT_PET_VEGNOCR,      /**< potential evap (transpiration only) from current vegetation and 0 canopy resistance [mm] (ALMA_OUTPUT: [mm/s]) */
    OUT_PREC,             /**< incoming precipitation [mm] (ALMA_OUTPUT: [mm/s]) */
    OUT_RAINF,            /**< rainfall  [mm] (ALMA_OUTPUT: [mm/s]) */
    OUT_REFREEZE,         /**< refreezing of water in the snow  [mm] (ALMA_OUTPUT: [mm/s]) */
    OUT_RUNOFF,           /**< surface runoff [mm] (ALMA_OUTPUT: [mm/s]) */
    OUT_SNOW_MELT,        /**< snow melt  [mm] (ALMA_OUTPUT: [mm/s]) */
    OUT_SNOWF,            /**< snowfall  [mm] (ALMA_OUTPUT: [mm/s]) */
    OUT_SUB_BLOWING,      /**< net sublimation of blowing snow [mm] (ALMA_OUTPUT: [mm/s]) */
    OUT_SUB_CANOP,        /**< net sublimation from snow stored in canopy [mm] (ALMA_OUTPUT: [mm/s]) */
    OUT_SUB_SNOW,         /**< total net sublimation from snow pack (surface and blowing) [mm] (ALMA_OUTPUT: [mm/s]) */
    OUT_SUB_SURFACE,      /**< net sublimation from snow pack surface [mm] (ALMA_OUTPUT: [mm/s]) */
    OUT_TRANSP_VEG,       /**< net transpiration from vegetation [mm] (ALMA_OUTPUT: [mm/s]) */
    OUT_WATER_ERROR,      /**< water budget error [mm] */
    // Energy Balance Terms - state variables
    OUT_ALBEDO,           /**< average surface albedo [fraction] */
    OUT_BARESOILT,        /**< bare soil surface temperature [C] (ALMA_OUTPUT: [K]) */
    OUT_FDEPTH,           /**< depth of freezing fronts [cm] (ALMA_OUTPUT: [m]) for each freezing front */
    OUT_LAKE_ICE_TEMP,    /**< temperature of lake ice [C] (ALMA_OUTPUT: [K]) */
    OUT_LAKE_SURF_TEMP,   /**< lake surface temperature [C] (ALMA_OUTPUT: [K]) */
    OUT_RAD_TEMP,         /**< average radiative surface temperature [K] */
    OUT_SALBEDO,          /**< snow pack albedo [fraction] */
    OUT_SNOW_PACK_TEMP,   /**< snow pack temperature [C] (ALMA_OUTPUT: [K]) */
    OUT_SNOW_SURF_TEMP,   /**< snow surface temperature [C] (ALMA_OUTPUT: [K]) */
    OUT_SNOWT_FBFLAG,     /**< snow surface temperature fallback flag */
    OUT_SOIL_TEMP,        /**< soil temperature [C] (ALMA_OUTPUT: [K]) for each soil layer */
    OUT_SOIL_TNODE,       /**< soil temperature [C] (ALMA_OUTPUT: [K]) for each soil thermal node */
    OUT_SOIL_TNODE_WL,    /**< soil temperature [C] (ALMA_OUTPUT: [K]) for each soil thermal node in the wetland */
    OUT_SOILT_FBFLAG,     /**< soil temperature flag for each soil thermal node */
    OUT_SURF_TEMP,        /**< average surface temperature [C] (ALMA_OUTPUT: [K]) */
    OUT_SURFT_FBFLAG,     /**< surface temperature flag */
    OUT_TCAN_FBFLAG,      /**< Tcanopy flag */
    OUT_TDEPTH,           /**< depth of thawing fronts [cm] (ALMA_OUTPUT: [m]) for each thawing front */
    OUT_TFOL_FBFLAG,      /**< Tfoliage flag */
    OUT_VEGT,             /**< average vegetation canopy temperature [C] (ALMA_OUTPUT: [K]) */
    // Energy Balance Terms - fluxes
    OUT_ADV_SENS,         /**< net sensible flux advected to snow pack [W/m2] */
    OUT_ADVECTION,        /**< advected energy [W/m2] */
    OUT_DELTACC,          /**< rate of change in cold content in snow pack [W/m2] (ALMA_OUTPUT: [J/m2]) */
    OUT_DELTAH,           /**< rate of change in heat storage [W/m2] (ALMA_OUTPUT: [J/m2]) */
    OUT_ENERGY_ERROR,     /**< energy budget error [W/m2] */
    OUT_FUSION,           /**< net energy used to melt/freeze soil moisture [W/m2] */
    OUT_GRND_FLUX,        /**< net heat flux into ground [W/m2] */
    OUT_IN_LONG,          /**< incoming longwave at ground surface (under veg) [W/m2] */
    OUT_LATENT,           /**< net upward latent heat flux [W/m2] */
    OUT_LATENT_SUB,       /**< net upward latent heat flux from sublimation [W/m2] */
    OUT_MELT_ENERGY,      /**< energy of fusion (melting) in snowpack [W/m2] */
    OUT_NET_LONG,         /**< net downward longwave flux [W/m2] */
    OUT_NET_SHORT,        /**< net downward shortwave flux [W/m2] */
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
    OUT_AIR_TEMP,         /**< air temperature [C] (ALMA_OUTPUT: [K])*/
    OUT_CATM,             /**< atmospheric CO2 concentrtaion [ppm]*/
    OUT_COSZEN,           /**< cosine of solar zenith angle [fraction]*/
    OUT_DENSITY,          /**< near-surface atmospheric density [kg/m3]*/
    OUT_FDIR,             /**< fraction of incoming shortwave that is direct [fraction]*/
    OUT_LAI,              /**< leaf area index [m2/m2] */
    OUT_LONGWAVE,         /**< incoming longwave [W/m2] */
    OUT_PAR,              /**< incoming photosynthetically active radiation [W/m2] */
    OUT_PRESSURE,         /**< near surface atmospheric pressure [kPa] (ALMA_OUTPUT: [Pa])*/
    OUT_QAIR,             /**< specific humidity [kg/kg] */
    OUT_REL_HUMID,        /**< relative humidity [%]*/
    OUT_SHORTWAVE,        /**< incoming shortwave [W/m2] */
    OUT_SURF_COND,        /**< surface conductance [m/s] */
    OUT_TSKC,             /**< cloud cover fraction [fraction] */
    OUT_VEGCOVER,         /**< fractional area of plants [fraction] */
    OUT_VP,               /**< near surface vapor pressure [kPa] (ALMA_OUTPUT: [Pa]) */
    OUT_VPD,              /**< near surface vapor pressure deficit [kPa] (ALMA_OUTPUT: [Pa]) */
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
    OUT_NET_LONG_BAND,    /**< net downward longwave flux [W/m2] */
    OUT_NET_SHORT_BAND,   /**< net downward shortwave flux [W/m2] */
    OUT_RFRZ_ENERGY_BAND, /**< net energy used to refreeze liquid water in snowpack [W/m2] */
    OUT_SENSIBLE_BAND,    /**< net upward sensible heat flux [W/m2] */
    OUT_SNOW_CANOPY_BAND, /**< snow interception storage in canopy [mm] */
    OUT_SNOW_COVER_BAND,  /**< fractional area of snow cover [fraction] */
    OUT_SNOW_DEPTH_BAND,  /**< depth of snow pack [cm] */
    OUT_SNOW_FLUX_BAND,   /**< energy flux through snow pack [W/m2] */
    OUT_SNOW_MELT_BAND,   /**< snow melt [mm] (ALMA_OUTPUT: [mm/s]) */
    OUT_SNOW_PACKT_BAND,  /**< snow pack temperature [C] (ALMA_OUTPUT: [K]) */
    OUT_SNOW_SURFT_BAND,  /**< snow surface temperature [C] (ALMA_OUTPUT: [K]) */
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

/***** Data Structures *****/

/******************************************************************************
 * @brief   file structures
 *****************************************************************************/
typedef struct {
    FILE *forcing[2];   /**< atmospheric forcing data files */
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
    char forcing[2][MAXSTRING];    /**< atmospheric forcing data file names */
    char f_path_pfx[2][MAXSTRING]; /**< path and prefix for atmospheric forcing data file names */
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
 * @brief   This structure stores model options.
 *****************************************************************************/
typedef struct {
    // simulation modes
    short AboveTreelineVeg;  /**< Default veg type to use above treeline;
                                Negative number indicates bare soil. */
    unsigned short AERO_RESIST_CANSNOW;  /**< "AR_406" = multiply aerodynamic resistance
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
    bool BLOWING;        /**< TRUE = calculate sublimation from blowing snow */
    bool BLOWING_VAR_THRESHOLD;
    bool BLOWING_CALC_PROB;
    bool BLOWING_SIMPLE;
    bool BLOWING_FETCH;
    bool BLOWING_SPATIAL_WIND;
    bool CARBON;         /**< TRUE = simulate carbon cycling processes;
                            FALSE = no carbon cycling (default) */
    bool CLOSE_ENERGY;   /**< TRUE = all energy balance calculations are
                            iterated to minimize the total column (air,
                            canopy, snow and ground) error; FALSE = no
                            iteration is used and the model estimates the new
                            fluxes based on those from the previous time step,
                            results should be similar, however, the model will
                            report energy balance errors. */
    bool COMPUTE_TREELINE; /**< TRUE = Determine treeline and exclude overstory
                              vegetation from higher elevations */
    bool CONTINUEONERROR; /**< TRUE = VIC will continue to run after a cell has an error */
    bool CORRPREC;       /**< TRUE = correct precipitation for gage undercatch */
    bool EQUAL_AREA;     /**< TRUE = RESOLUTION stores grid cell area in km^2;
                            FALSE = RESOLUTION stores grid cell side length in degrees */
    bool EXP_TRANS;      /**< TRUE = Uses grid transform for exponential node
                            distribution for soil heat flux calculations*/
    bool FROZEN_SOIL;    /**< TRUE = Use frozen soils code */
    bool FULL_ENERGY;    /**< TRUE = Use full energy code */
    unsigned short GRND_FLUX_TYPE; /**< "GF_406"  = use (flawed) formulas for ground flux, deltaH, and fusion
                                        from VIC 4.0.6 and earlier
                                      "GF_410"  = use formulas from VIC 4.1.0 */
    bool IMPLICIT;       /**< TRUE = Use implicit solution when computing
                            soil thermal fluxes */
    bool JULY_TAVG_SUPPLIED; /**< If TRUE and COMPUTE_TREELINE is also true,
                                then average July air temperature will be read
                                from soil file and used in calculating treeline */
    bool LAKES;          /**< TRUE = use lake energy code */
    unsigned short LW_CLOUD;       /**< Longwave cloud formulation; "LW_CLOUD_x" = code for LW cloud formulation - see LW_CLOUD codes above */
    unsigned short LW_TYPE;        /**< Longwave clear sky algorithm; "LW_x" = code for LW algorithm - see LW codes above */
    bool MTCLIM_SWE_CORR; /**< TRUE = correct MTCLIM's downward shortwave radiation estimate in presence of snow */
    size_t Ncanopy;      /**< Number of canopy layers in the model. */
    size_t Nfrost;       /**< Number of frost subareas in model */
    size_t Nlakenode;    /**< Number of lake thermal nodes in the model. */
    size_t Nlayer;       /**< Number of layers in model */
    size_t Nnode;        /**< Number of soil thermal nodes in the model */
    bool NOFLUX;         /**< TRUE = Use no flux lower bondary when computing
                            soil thermal fluxes */
    size_t NVEGTYPES;    /**< number of vegetation types (used by image driver) */
    bool PLAPSE;         /**< TRUE = If air pressure not supplied as an
                            input forcing, compute it by lapsing sea-level
                            pressure by grid cell average elevation;
                            FALSE = air pressure set to constant 95.5 kPa */
    unsigned short RC_MODE;        /**< RC_JARVIS = compute canopy resistance via Jarvis formulation (default)
                                      RC_PHOTO = compute canopy resistance based on photosynthetic activity */
    size_t ROOT_ZONES;   /**< Number of root zones used in simulation */
    bool QUICK_FLUX;     /**< TRUE = Use Liang et al., 1999 formulation for
                            ground heat flux, if FALSE use explicit finite
                            difference method */
    bool QUICK_SOLVE;    /**< TRUE = Use Liang et al., 1999 formulation for
                            iteration, but explicit finite difference
                            method for final step. */
    bool SHARE_LAYER_MOIST; /**< TRUE = transpiration in moisture-limited layers can draw from other layers (default) */
    unsigned short SNOW_DENSITY;   /**< DENS_BRAS: Use algorithm of Bras, 1990; DENS_SNTHRM: Use algorithm of SNTHRM89 adapted for 1-layer pack */
    size_t SNOW_BAND;    /**< Number of elevation bands over which to solve the
                            snow model */
    unsigned SNOW_STEP;  /**< Time step in hours to use when solving the snow model */
    bool SPATIAL_FROST;   /**< TRUE = use a uniform distribution to simulate the
                             spatial distribution of soil frost; FALSE = assume
                             that the entire grid cell is frozen uniformly. */
    bool SPATIAL_SNOW;    /**< TRUE = use a uniform distribution to simulate the
                             partial coverage of the surface by a thin snowpack.
                             Coverage is assumed to be uniform after snowfall
                             until the pack begins to melt. */
    bool TFALLBACK;      /**< TRUE = when any temperature iterations fail to converge,
                                   use temperature from previous time step; the number
                                   of instances when this occurs will be logged and
                                   reported at the end of the cell's simulation
                            FALSE = when iterations fail to converge, report an error
                                    and abort simulation for current grid cell
                            Default = TRUE */
    bool VP_INTERP;      /**< How to disaggregate VP from daily to sub-daily;
                            TRUE = linearly interpolate between daily VP values, assuming they occur at the times of Tmin;
                            FALSE = hold VP constant at the daily value */
    unsigned short VP_ITER;        /**< VP_ITER_NONE = never iterate with SW
                                      VP_ITER_ALWAYS = always iterate with SW
                                      VP_ITER_ANNUAL = use annual Epot/PRCP criterion
                                      VP_ITER_CONVERGE = always iterate until convergence */

    // input options
    bool ALMA_INPUT;     /**< TRUE = input variables are in ALMA-compliant units; FALSE = standard VIC units */
    bool BASEFLOW;       /**< ARNO: read Ds, Dm, Ws, c; NIJSSEN2001: read d1, d2, d3, d4 */
    unsigned short GRID_DECIMAL; /**< Number of decimal places in grid file extensions */
    bool VEGLIB_PHOTO;   /**< TRUE = veg library contains photosynthesis parameters */
    bool VEGLIB_VEGCOVER; /**< TRUE = veg library file contains monthly vegcover values */
    bool VEGPARAM_ALB;   /**< TRUE = veg param file contains monthly albedo values */
    bool VEGPARAM_LAI;   /**< TRUE = veg param file contains monthly LAI values */
    bool VEGPARAM_VEGCOVER; /**< TRUE = veg param file contains monthly vegcover values */
    unsigned short ALB_SRC;        /**< FROM_VEGLIB = use albedo values from veg library file
                                      FROM_VEGPARAM = use albedo values from the veg param file */
    unsigned short LAI_SRC;        /**< FROM_VEGLIB = use LAI values from veg library file
                                      FROM_VEGPARAM = use LAI values from the veg param file */
    unsigned short VEGCOVER_SRC;   /**< FROM_VEGLIB = use vegcover values from veg library file
                                      FROM_VEGPARAM = use vegcover values from the veg param file */
    bool LAKE_PROFILE;   /**< TRUE = user-specified lake/area profile */
    bool ORGANIC_FRACT;  /**< TRUE = organic matter fraction of each layer is read from the soil parameter file; otherwise set to 0.0. */

    // state options
    bool BINARY_STATE_FILE; /**< TRUE = model state file is binary (default) */
    bool INIT_STATE;     /**< TRUE = initialize model state from file */
    bool SAVE_STATE;     /**< TRUE = save state file */

    // output options
    bool ALMA_OUTPUT;    /**< TRUE = output variables are in ALMA-compliant units; FALSE = standard VIC units */
    bool BINARY_OUTPUT;  /**< TRUE = output files are in binary, not ASCII */
    bool COMPRESS;       /**< TRUE = Compress all output files */
    bool MOISTFRACT;     /**< TRUE = output soil moisture as fractional moisture content */
    size_t Noutfiles;    /**< Number of output files (not including state files) */
    bool OUTPUT_FORCE;   /**< TRUE = perform disaggregation of forcings, skip
                            the simulation, and output the disaggregated
                            forcings. */
    bool PRT_HEADER;     /**< TRUE = insert header at beginning of output file; FALSE = no header */
    bool PRT_SNOW_BAND;  /**< TRUE = print snow parameters for each snow band. This is only used when default
                                   output files are used (for backwards-compatibility); if outfiles and
                                   variables are explicitly mentioned in global parameter file, this option
                                   is ignored. */
} option_struct;

/******************************************************************************
 * @brief   This structure stores all model run global parameters.
 *****************************************************************************/
typedef struct {
    double measure_h;              /**< height of measurements (m) */
    double wind_h;                 /**< height of wind measurements (m) */
    double resolution;             /**< Model resolution (degrees) */
    unsigned dt;             /**< Time step in hours (24/dt must be an
                                      integer) */
    unsigned out_dt;          /**< Output time step in hours (24/out_dt must
                                      be an integer) */
    unsigned short endday;          /**< Last day of model simulation */
    unsigned short endmonth;        /**< Last month of model simulation */
    unsigned short endyear;         /**< Last year of model simulation */
    unsigned short forceday[2];    /**< day forcing files starts */
    unsigned short forcehour[2];   /**< hour forcing files starts */
    unsigned short forcemonth[2];  /**< month forcing files starts */
    unsigned short forceoffset[2]; /**< counter to keep track of offset in reading
                                      forcing files; updated after every read */
    unsigned short forceskip[2];   /**< number of model time steps to skip at
                                      the start of the forcing file */
    unsigned short forceyear[2];   /**< year forcing files start */
    unsigned nrecs;                /**< Number of time steps simulated */
    unsigned short skipyear;       /**< Number of years to skip before writing
                                      output data */
    unsigned short startday;       /**< Starting day of the simulation */
    unsigned short starthour;      /**< Starting hour of the simulation */
    unsigned short startmonth;     /**< Starting month of the simulation */
    unsigned short startyear;      /**< Starting year of the simulation */
    unsigned short stateday;       /**< Day of the simulation at which to save
                                      model state */
    unsigned short statemonth;     /**< Month of the simulation at which to save
                                      model state */
    unsigned short stateyear;      /**< Year of the simulation at which to save
                                      model state */
} global_param_struct;

/******************************************************************************
 * @brief    This structure holds the model parameters.
 *****************************************************************************/
typedef struct {
    // Lapse Rate
    double LAPSE_RATE;  /**< temperature lapse rate (C/m) */

    // Precipitation Guage Height
    double GAUGE_HEIGHT;   /**< precipitation gauge height (m) */

    // Default Wind Speed
    double WIND_SPEED_DEFAULT;  /**< Default wind speed (m/s) used when wind is not supplied as a forcing */
    double WIND_SPEED_MIN;  /**< Minimum wind speed in m/s that can be used by the model. */

    // Huge Resistance Term
    double HUGE_RESIST;  /**< Extermely large resistance term (s/m) */

    // Surface Albedo Parameters
    double ALBEDO_BARE_SOIL;  /**< Broadband albedo of bare soil */
    double ALBEDO_H20_SURF;  /**< Broadband albedo of open water surface */

    // Surface Emissivities
    double EMISS_GRND;  /**< Emissivity of bare soil */
    double EMISS_VEG;  /**< Emissivity of vegetation */
    double EMISS_ICE;  /**< Emissivity of bare ice */
    double EMISS_SNOW;  /**< Emissivity of snow */
    double EMISS_H2O;  /**< Emissivity of open water surface */

    // Soil Constraints
    double SOIL_RESID_MOIST;  /**< Default residual moisture content of soil colum */
    double SOIL_SLAB_MOIST_FRACT;  /**< Volumetric moisture content (fraction of porosity) in the soil/rock below the bottom soil layer; this assumes that the soil below the bottom layer has the same texture as the bottom layer. */

    // Vegetation Parameters
    double VEG_LAI_SNOW_MULTIPLIER;  /**< multiplier to calculate the amount of available snow interception as a function of LAI (m) */
    double VEG_MIN_INTERCEPTION_STORAGE;  /**< the amount of snow on the canopy that can only be melted off. (m) */
    double VEG_LAI_WATER_FACTOR;  /**< Coefficient multiplied by the LAI to determine the amount of water that can be stored in the canopy */

    // Canopy Parameters
    double CANOPY_CLOSURE;  /**< Threshold vapor pressure deficit for stomatal closure (Pa) */
    double CANOPY_RSMAX;  /**< Maximum allowable resistance (s/m) */
    double CANOPY_VPDMINFACTOR;  /**< Minimum allowable vapor pressure deficit factor */

    // MTCLIM Parameters
    double MTCLIM_TDAYCOEF;  /**< (dim) daylight air temperature coefficient */
    double MTCLIM_SOLAR_CONSTANT;  /**< solar constants (W/m2) */
    double MTCLIM_SNOW_TCRIT;  /**<  (deg C) critical temperature for snowmelt */
    double MTCLIM_SNOW_TRATE;  /**< (cm/degC/day) snowmelt rate */
    double MTCLIM_TBASE;  /**< (dim) max inst. trans., 0m, nadir, dry atm */
    double MTCLIM_ABASE;  /**< (1/Pa) vapor pressure effect on transmittance */
    double MTCLIM_C;  /**< (dim) radiation parameter */
    double MTCLIM_B0;  /**< (dim) radiation parameter */
    double MTCLIM_B1;  /**< (dim) radiation parameter */
    double MTCLIM_B2; // (dim) radiation parameter
    double MTCLIM_RAIN_SCALAR;  /**< (dim) correction to trans. for rain day */
    double MTCLIM_DIF_ALB;  /**< (dim) diffuse albedo for horizon correction */
    double MTCLIM_SC_INT;  /**< (MJ/m2/day) snow correction intercept */
    double MTCLIM_SC_SLOPE;  /**<  (MJ/m2/day) snow correction intercept */
    double MTCLIM_SRADDT;  /**< timestep for radiation routine (seconds).  Note: Make sure that 3600 % SRADDT == 0 */
    double MTCLIM_SW_PREC_THRESH;  /**< Minimum daily precip (mm) that can cause dimming of incoming shortwave; default = 0. */

    // Lake Parameters
    double LAKE_TMELT;
    double LAKE_MAX_SURFACE;  /**< max. surface layer thickness for E-B (m) */
    double LAKE_BETA;  /**< Curve shape parameter for lake profile. */
    double LAKE_FRACMIN;  /**< min ice thickness in meters */
    double LAKE_FRACLIM;  /**< lower limit on fractional ice cover */
    double LAKE_DM;  /**< molecular diffusivity of water */
    double LAKE_SNOWCRIT;  /**< for albedo (m) */
    double LAKE_ZWATER;
    double LAKE_ZSNOW;
    double LAKE_RHOSNOW;  /**< density of snow (kg m-3) */
    double LAKE_CONDI;  /**< thermal conductivity of ice */
    double LAKE_CONDS;  /**< thermal conductivity of snow */
    double LAKE_LAMISW;  /**< attenuation of shortwave radiation through ice (1/m) */
    double LAKE_LAMILW;  /**< attenuation of longwave radiation through ice (1/m) */
    double LAKE_LAMSSW;  /**< attenuation of shortwave radiation through snow (1/m) */
    double LAKE_LAMSLW;  /**< attenuation of longwave radiation through snow (1/m) */
    double LAKE_LAMWSW;  /**< attenuation of shortwave radiation through water (1/m) */
    double LAKE_LAMWLW;  /**< attenuation of longwave radiation through water (1/m) */
    double LAKE_A1;  /**< Percent of radiation in visible band. */
    double LAKE_A2;  /**< Percent of radiation in infrared band. */
    double LAKE_QWTAU;  /**< D. Pollard sub-ice time constant. */
    int LAKE_MAX_ITER;

    // Saturation Vapor Pressure Parameters
    double SVP_A;  /**< constant for saturated vapor pressure curve (kPa) */
    double SVP_B;  /**< constant for saturated vapor pressure curve (kPa) */
    double SVP_C;  /**< constant for saturated vapor pressure curve (kPa) */

    // Carbon Parameters
    double CARBON_CATMCURRENT;  /**< Current global atmospheric CO2 mixing ratio (ppm) */
    double CARBON_SW2PAR;  /**< Empirical ratio of PAR (W/m2) to SHORTWAVE (W/m2) from Lopez et al., 2001 */

    // Photosynthesis Parameters
    double PHOTO_OMEGA;  /**< single leaf scattering albedo */
    double PHOTO_LAIMAX;  /**< Maximum LAI in nitrogen scaling */
    double PHOTO_LAILIMIT;  /**< Minimum LAI in nitrogen scaling and maximum LAI in PAR computation */
    double PHOTO_LAIMIN;  /**< Minimum LAI in PAR computation */
    double PHOTO_EPAR;  /**< Energy content of PAR [J/mol photons] = (4.6 mol/MJ PAR)^-1 */
    double PHOTO_FCMAX;  /**< Maximum fractional veg cover; (1-FcMax) = min amount of ground visible */
    double PHOTO_FCMIN;  /**< Minimum fractional veg cover; (1-FcMin) = max amount of ground visible */
    double PHOTO_ZENITHMIN;  /**< Check for solar zenith angle > 89 deg */
    double PHOTO_ZENITHMINPAR;  /**< Cosine of the minimum solar zenith angle for photosynthesis to take place */
    double PHOTO_ALBSOIPARMIN;  /**< Minimum soil reflectivity in PAR range */
    double PHOTO_MINMAXETRANS;  /**< Minimum of maximum electron transport rate [10e-12 mol/(m^2 s)] */
    double PHOTO_MINSTOMCOND;  /**< Minimum stomatal conductance [mol H2O/m2s] */
    double PHOTO_FCI1C3;  /**< C3 Plants factor that relate leaf internal CO2 concentration to ambient CO2 concentration */
    double PHOTO_FCI1C4;  /**< C4 Plants factor that relate leaf internal CO2 concentration to ambient CO2 concentration */
    double PHOTO_OX;  /**< OXYGEN CONCENTRATION [MOL(O2) / MOL(AIR)] */
    double PHOTO_KC;  /**< MICHAELIS-MENTEN CONSTANT FOR CO2 AT 25C [MOL(CO2) / MOL(AIR)] */
    double PHOTO_KO;  /**< MICHAELIS-MENTEN CONSTANT FOR O2 AT 25C [MOL(O2) / MOL(AIR)] */
    double PHOTO_EC;  /**< ACTIVATION ENERGY FOR KC [J / MOL] */
    double PHOTO_EO;  /**< ACTIVATION ENERGY FOR KO [J / MOL] */
    double PHOTO_EV;  /**< ACTIVATION ENERGY FOR VCMAX [J / MOL] */
    double PHOTO_ER;  /**< ACTIVATION ENERGY FOR DARK RESPIRATION [J / MOL] */
    double PHOTO_ALC3;  /**< EFFICIENCY OF OF PHOTON CAPTURE */
    double PHOTO_FRDC3;  /**< RATIO OF DARK RESPIRATION TO "PVM" AT 25C for C3 */
    double PHOTO_EK;  /**< = Q10=2 (Collatz et al. 1992) */
    double PHOTO_ALC4;  /**< EFFECTIVE QUANTUM EFFICIENCY */
    double PHOTO_FRDC4;  /**< RATIO OF DARK RESPIRATION TO "PVM" AT 25C for C4 */
    double PHOTO_THETA;  /**< CURVATURE PARAMETER */
    double PHOTO_FRLEAF;  /**< Ratio of canopy leaf respiration to whole plant maintenance respiration */
    double PHOTO_FRGROWTH;  /**< Ratio of plant growth respiration to NPP */

    // Soil Respiration Parameters
    double SRESP_E0_LT;  /**< Lloyd-Taylor E0 parameter [K] */
    double SRESP_T0_LT;  /**< Lloyd-Taylor T0 parameter [K] */
    double SRESP_WMINFM;  /**< minimum soil moisture (fraction) at which soil respiration can occur */
    double SRESP_WMAXFM;  /**< maximum soil moisture (fraction) at which soil respiration can occur */
    double SRESP_WOPTFM;  /**< soil moisture (fraction) at which maximum soil respiration occurs */
    double SRESP_RHSAT;  /**< ratio of soil respiration rate under saturated conditions (w=wmaxFM) to that under optimal conditions (w=woptFM) */
    double SRESP_RFACTOR;  /**< scaling factor to account for other (non-moisture) sources of inhibition of respiration */
    double SRESP_TAULITTER;  /**< Litter pool turnover time [y] */
    double SRESP_TAUINTER;  /**< Intermediate pool turnover time [y] */
    double SRESP_TAUSLOW;  /**< Slow pool turnover time [y] */
    double SRESP_FAIR;  /**< Fraction of respired carbon from litter pool that is lost to atmosphere */
    double SRESP_FINTER;  /**< Fraction of [respired carbon from litter pool that goes to soil] that goes to intermediate pool */

    // Snow Parameters
    double SNOW_MAX_SURFACE_SWE;  /**< maximum depth of the surface layer in water equivalent (m) */
    double SNOW_LIQUID_WATER_CAPACITY;  /**< water holding capacity of snow as a fraction of snow-water-equivalent */
    double SNOW_NEW_SNOW_DENSITY;  /**< density of new fallen snow */
    double SNOW_DENS_DMLIMIT;  /**< Density limit used in calculation of destructive metamorphism (kg/m^3) */
    double SNOW_DENS_MAX_CHANGE;  /**< maximum change in snowfall depth (fraction of swe) */
    double SNOW_DENS_ETA0;  /**< viscosity of snow at T = 0C and density = 0 used in calculation of true viscosity (Ns/m2) */
    double SNOW_DENS_C1;  /**< Constant in snow density computation */
    double SNOW_DENS_C2;  /**< Constant in snow density computation */
    double SNOW_DENS_C5;  /**< constant used in snow viscosity calculation, taken from SNTHRM.89 (/C) */
    double SNOW_DENS_C6;  /**< constant used in snow viscosity calculation, taken from SNTHRM.89 (kg/m3) */
    double SNOW_DENS_F;  /**< internal compaction rate coefficient */
    double SNOW_MIN_SWQ_EB_THRES;  /**< Minimum SWQ for which the snowpack energy balance is computed independent of the soil surface temperature */
    double SNOW_A1;  /**< Attenuation coefficient for shortwave in a snowpack. Value and equation taken from Patterson and Hamblin, 1988 */
    double SNOW_A2;  /**< Attenuation coefficient for shortwave in a snowpack. Value and equation taken from Patterson and Hamblin, 1988 */
    double SNOW_L1;  /**< Attenuation coefficient for shortwave in a snowpack. Value and equation taken from Patterson and Hamblin, 1988 (1/m) */
    double SNOW_L2;  /**< Attenuation coefficient for shortwave in a snowpack. Value and equation taken from Patterson and Hamblin, 1988 (1/m) */
    double SNOW_NEW_SNOW_ALB;  /**< Snow albedo curve parameters. */
    double SNOW_ALB_ACCUM_A;  /**< Snow albedo curve parameters. */
    double SNOW_ALB_ACCUM_B;  /**< Snow albedo curve parameters. */
    double SNOW_ALB_THAW_A;  /**< Snow albedo curve parameters. */
    double SNOW_ALB_THAW_B;  /**< Snow albedo curve parameters. */
    double SNOW_TRACESNOW;  /**< Defines the minimum amount of new snow (mm) which will reset the snowpack albedo to new snow */
    double SNOW_CONDUCT;  /**< conductivity of snow (W/mK) */
    double SNOW_MAX_SNOW_TEMP;  /**< maximum temperature (C) at which snow can fall */
    double SNOW_MIN_RAIN_TEMP;  /**< minimum temperature (C) at which rain can fall */

    // Blowing Snow Parameters
    double BLOWING_KA;  /**< thermal conductivity of air (W/mK) */
    double BLOWING_CSALT;  /**< saltation constant m/s */
    double BLOWING_UTHRESH;  /**< threshold shear velocity m/s */
    double BLOWING_KIN_VIS;  /**< Kinemativ viscosity of air (m2/s) */
    int BLOWING_MAX_ITER;     /**< Max. iterations for numerical integration */
    int BLOWING_K;
    double BLOWING_SETTLING;  /**< Particle settling velocity m/s */
    int BLOWING_NUMINCS;     /**< Number of prob intervals to solve for wind. */

    // Treeline temperature
    double TREELINE_TEMPERATURE;  /**< Number of prob intervals to solve for wind. */

    // Iteration Bracket Widths
    double SNOW_DT;  /**< Used to bracket snow surface temperatures while computing the snow surface energy balance (C) */
    double SURF_DT;  /**< Used to bracket soil surface temperatures while computing energy balance (C) */
    double SOIL_DT;  /**< Used to bracket soil temperatures while solving the soil thermal flux (C) */
    double CANOPY_DT;  /**< Used to bracket canopy air temperatures while computing energy balance (C) */
    double CANOPY_VP;  /**< Used to bracket canopy vapor pressures while computing moisture balance (Pa) */

    // Convergence Tolerances
    double TOL_GRND;
    double TOL_OVER;

    // Frozen Soil Parameters
    int FROZEN_MAXITER;

    // Newton-Raphson Solver Parameters
    int NEWT_RAPH_MAXTRIAL;
    double NEWT_RAPH_TOLX;
    double NEWT_RAPH_TOLF;
    double NEWT_RAPH_R_MAX;
    double NEWT_RAPH_R_MIN;
    double NEWT_RAPH_RELAX1;
    double NEWT_RAPH_RELAX2;
    double NEWT_RAPH_RELAX3;
    double NEWT_RAPH_EPS2;

    // Root-Brent parameters
    int ROOT_BRENT_MAXTRIES;
    int ROOT_BRENT_MAXITER;
    double ROOT_BRENT_TSTEP;
    double ROOT_BRENT_T;
} parameters_struct;

/******************************************************************************
 * @brief   This structure stores the soil parameters for a grid cell.
 *****************************************************************************/
typedef struct {
    bool FS_ACTIVE;                   /**< if TRUE frozen soil algorithm is
                                         active in current grid cell */
    double Ds;                        /**< fraction of maximum subsurface flow
                                         rate */
    double Dsmax;                     /**< maximum subsurface flow rate
                                         (mm/day) */
    double Ksat[MAX_LAYERS];          /**< saturated hydraulic  conductivity
                                         (mm/day) */
    double Wcr[MAX_LAYERS];           /**< critical moisture level for soil
                                         layer, evaporation is no longer
                                         affected moisture stress in the
                                         soil (mm) */
    double Wpwp[MAX_LAYERS];          /**< soil moisture content at permanent
                                         wilting point (mm) */
    double Ws;                        /**< fraction of maximum soil moisture */
    double AlbedoPar;                 /**< soil albedo in PAR range (400-700nm) */
    double alpha[MAX_NODES];          /**< thermal solution constant */
    double annual_prec;               /**< annual average precipitation (mm) */
    double avg_temp;                  /**< average soil temperature (C) */
    double avgJulyAirTemp;            /**< Average July air temperature (C) */
    double b_infilt;                  /**< infiltration parameter */
    double beta[MAX_NODES];           /**< thermal solution constant */
    double bubble[MAX_LAYERS];        /**< bubbling pressure, HBH 5.15 (cm) */
    double bubble_node[MAX_NODES];    /**< bubbling pressure (cm) */
    double bulk_density[MAX_LAYERS];  /**< soil bulk density (kg/m^3) */
    double bulk_dens_min[MAX_LAYERS]; /**< bulk density of mineral soil (kg/m^3) */
    double bulk_dens_org[MAX_LAYERS]; /**< bulk density of organic soil (kg/m^3) */
    double c;                         /**< exponent in ARNO baseflow scheme */
    double depth[MAX_LAYERS];         /**< thickness of each soil moisture layer (m) */
    double dp;                        /**< soil thermal damping depth (m) */
    double dz_node[MAX_NODES];        /**< thermal node thickness (m) */
    double Zsum_node[MAX_NODES];      /**< thermal node depth (m) */
    double expt[MAX_LAYERS];          /**< layer-specific exponent n (=3+2/lambda) in Campbell's eqn for hydraulic conductivity, HBH 5.6 */
    double expt_node[MAX_NODES];      /**< node-specific exponent n (=3+2/lambda) in Campbell's eqn for hydraulic conductivity, HBH 5.6 */
    double frost_fract[MAX_FROST_AREAS]; /**< spatially distributed frost coverage fractions */
    double frost_slope;               /**< slope of frost distribution */
    double gamma[MAX_NODES];          /**< thermal solution constant */
    double init_moist[MAX_LAYERS];    /**< initial layer moisture level (mm) */
    double max_infil;                 /**< maximum infiltration rate */
    double max_moist[MAX_LAYERS];     /**< maximum moisture content (mm) per layer */
    double max_moist_node[MAX_NODES]; /**< maximum moisture content (mm/mm) per node */
    double max_snow_distrib_slope;    /**< Maximum slope of snow depth distribution [m].  This should equal 2*depth_min, where depth_min = minimum snow pack depth below which coverage < 1.  Comment, ported from user_def.h, with questionable units: SiB uses 0.076; Rosemount data imply 0.155cm depth ~ 0.028mm swq. */
    double phi_s[MAX_LAYERS];         /**< soil moisture diffusion parameter (mm/mm) */
    double porosity[MAX_LAYERS];      /**< porosity (fraction) */
    double quartz[MAX_LAYERS];        /**< quartz content of soil (fraction of mineral soil volume) */
    double organic[MAX_LAYERS];       /**< organic content of soil (fraction of total soil volume) */
    double resid_moist[MAX_LAYERS];   /**< residual moisture content of soil layer */
    double rough;                     /**< soil surface roughness (m) */
    double snow_rough;                /**< snow surface roughness (m) */
    double soil_density[MAX_LAYERS];  /**< soil particle density (kg/m^3) */
    double soil_dens_min[MAX_LAYERS]; /**< particle density of mineral soil (kg/m^3) */
    double soil_dens_org[MAX_LAYERS]; /**< particle density of organic soil (kg/m^3) */
    double *BandElev;                 /**< Elevation of each snow elevation band */
    double *AreaFract;                /**< Fraction of grid cell included in each snow elevation band */
    double *Pfactor;                  /**< Change in Precipitation due to elevation (fract) in each snow elevation band */
    double *Tfactor;                  /**< Change in temperature due to elevation (C) in each snow elevation band */
    bool *AboveTreeLine;             /**< Flag to indicate if band is above the treeline */
    double elevation;                 /**< grid cell elevation (m) */
    double lat;                       /**< grid cell central latitude */
    double lng;                       /**< grid cell central longitude */
    double cell_area;                 /**< Area of grid cell (m^2) */
    double time_zone_lng;             /**< central meridian of the time zone */
    unsigned gridcel;                 /**< grid cell number */
    double zwtvmoist_zwt[MAX_LAYERS + 2][MAX_ZWTVMOIST]; /**< zwt values in the zwt-v-moist curve for each layer */
    double zwtvmoist_moist[MAX_LAYERS + 2][MAX_ZWTVMOIST]; /**< moist values in the zwt-v-moist curve for each layer */
    double slope;
    double aspect;
    double ehoriz;
    double whoriz;
} soil_con_struct;

/******************************************************************************
 * @brief   This structure stores information about the vegetation coverage of
 *          the current grid cell.
 *****************************************************************************/
typedef struct {
    double Cv;              /**< fraction of vegetation coverage */
    double Cv_sum;          /**< total fraction of vegetation coverage */
    double root[MAX_LAYERS]; /**< percent of roots in each soil layer (fraction) */
    double *zone_depth;     /**< depth of root zone */
    double *zone_fract;     /**< fraction of roots within root zone */
    int veg_class;          /**< vegetation class reference number */
    size_t vegetat_type_num; /**< number of vegetation types in the grid cell */
    double sigma_slope;     /**< Std. deviation of terrain slope for each vegetation class. */
    double lag_one;         /**< Lag one gradient autocorrelation of terrain slope */
    double fetch;           /**< Average fetch length for each vegetation class. */
    int LAKE;               /**< TRUE = this tile is a lake/wetland tile */
    double *CanopLayerBnd;  /**< Upper boundary of each canopy layer, expressed as fraction of total LAI */
} veg_con_struct;

/******************************************************************************
 * @brief   This structure stores parameters for individual vegetation types.
 *****************************************************************************/
typedef struct {
    bool overstory;        /**< TRUE = overstory present, important for snow
                              accumulation in canopy */
    double LAI[MONTHS_PER_YEAR];  /**< leaf area index */
    double vegcover[MONTHS_PER_YEAR];  /**< fractional area covered by plants within the tile (fraction) */
    double Wdmax[MONTHS_PER_YEAR];  /**< maximum dew holding capacity (mm) */
    double albedo[MONTHS_PER_YEAR];  /**< vegetation albedo (added for full energy)
                                                           (fraction) */
    double displacement[MONTHS_PER_YEAR]; /**< vegetation displacement height (m) */
    double emissivity[MONTHS_PER_YEAR]; /**< vegetation emissivity (fraction) */
    size_t NVegLibTypes;   /**< number of vegetation classes defined in library */
    double rad_atten;      /**< radiation attenuation due to canopy,
                              default = 0.5 (N/A) */
    double rarc;           /**< architectural resistance (s/m) */
    double rmin;           /**< minimum stomatal resistance (s/m) */
    double roughness[MONTHS_PER_YEAR];  /**< vegetation roughness length (m) */
    double trunk_ratio;    /**< ratio of trunk height to tree height,
                                              default = 0.2 (fraction) */
    double wind_atten;     /**< wind attenuation through canopy,
                              default = 0.5 (N/A) */
    double wind_h;         /**< height at which wind is measured (m) */
    double RGL;            /**< Value of solar radiation below which there
                              will be no transpiration (ranges from
                              ~30 W/m^2 for trees to ~100 W/m^2 for crops) */
    unsigned short veg_class; /**< vegetation class reference number */
    char Ctype;              /**< Photosynthetic pathway; can be C3 or C4 */
    double MaxCarboxRate;  /**< maximum carboxlyation rate at 25 deg C (mol(CO2)/m2s) */
    double MaxETransport;  /**< maximum electron transport rate at 25 deg C (mol(CO2)/m2s) (C3 plants) */
    double CO2Specificity; /**< CO2 specificity at 25 deg C (mol(CO2)/m2s) (C4 plants) */
    double LightUseEff;    /**< Light-use efficiency (mol(CO2)/mol(photons)) */
    bool NscaleFlag;         /**< TRUE = nitrogen-scaling factors are applicable
                                to this veg class */
    double Wnpp_inhib;     /**< moisture level (fraction of maximum moisture) above which photosynthesis
                              experiencing saturation inhibition, i.e. too wet for optimal photosynthesis;
                              only applies to top soil layer */
    double NPPfactor_sat;  /**< photosynthesis multiplier (fraction of maximum) when top soil
                              layer is saturated */
} veg_lib_struct;

/******************************************************************************
 * @brief   This structure stores historical timeseries of vegetation
 *          parameters for a given vegetation tile.  Structure is similar to
 *          atmos_data_struct.
 *****************************************************************************/
typedef struct {
    double *albedo;  /**< timeseries of vegetation albedo (fraction) */
    double *LAI;     /**< timeseries of leaf area index (m2/m2) */
    double *vegcover; /**< timeseries of fractional area of plants within veg tile (fraction) */
} veg_hist_struct;

/******************************************************************************
 * @brief   This structure stores the atmospheric forcing data for each model
 * time step for a single grid cell.  Each array stores the values for the
 * SNOW_STEPs during the current model step and the value for the entire model
 * step.  The latter is referred to by array[NR].  Looping over the SNOW_STEPs
 * is done by for (i = 0; i < NF; i++)
 *****************************************************************************/
typedef struct {
    double *air_temp; /**< air temperature (C) */
    double *Catm;    /**< atmospheric CO2 mixing ratio (mol CO2/ mol air) */
    double *channel_in; /**< incoming channel inflow for time step (mm) */
    double *coszen;  /**< cosine of solar zenith angle (fraction) */
    double *density; /**< atmospheric density (kg/m^3) */
    double *fdir;    /**< fraction of incoming shortwave that is direct (fraction) */
    double *longwave; /**< incoming longwave radiation (W/m^2) (net incoming
                         longwave for water balance model) */
    double out_prec;  /**< Total precipitation for time step - accounts
                         for corrected precipitation totals */
    double out_rain;  /**< Rainfall for time step (mm) */
    double out_snow;  /**< Snowfall for time step (mm) */
    double *par;     /**< incoming photosynthetically active radiation () */
    double *prec;    /**< average precipitation in grid cell (mm) */
    double *pressure; /**< atmospheric pressure (kPa) */
    double *shortwave; /**< incoming shortwave radiation (W/m^2) */
    bool *snowflag;    /**< TRUE if there is snowfall in any of the snow
                          bands during the timestep, FALSE otherwise*/
    double *tskc;    /**< cloud cover fraction (fraction) */
    double *vp;      /**< atmospheric vapor pressure (kPa) */
    double *vpd;     /**< atmospheric vapor pressure deficit (kPa) */
    double *wind;    /**< wind speed (m/s) */
} atmos_data_struct;

/******************************************************************************
 * @brief   This structure stores information about the time and date of the
 *          current time step.
 *****************************************************************************/
typedef struct {
    unsigned short day;         /**< current day */
    unsigned short day_in_year; /**< julian day in year */
    unsigned short hour;        /**< beginning of current hour */
    unsigned short month;       /**< current month */
    unsigned short year;        /**< current year */
} dmy_struct;                   /**< array of length nrec created */

/******************************************************************************
 * @brief   This structure stores all soil variables for each layer in the
 *          soil column.
 *****************************************************************************/
typedef struct {
    double Cs;              /**< average volumetric heat capacity of the
                               current layer (J/m^3/K) */
    double T;               /**< temperature of the unfrozen sublayer (C) */
    double bare_evap_frac;  /**< fraction of evapotranspiration coming from bare soil evap, from soil layer (mm) */
    double evap;            /**< evapotranspiration from soil layer (mm) */
    double ice[MAX_FROST_AREAS]; /**< ice content of the frozen sublayer (mm) */
    double kappa;           /**< average thermal conductivity of the current
                               layer (W/m/K) */
    double moist;           /**< moisture content of the unfrozen sublayer
                               (mm) */
    double phi;             /**< moisture diffusion parameter */
    double zwt;             /**< water table position relative to soil surface within the layer (cm) */
} layer_data_struct;

/******************************************************************************
 * @brief   This structure stores soil variables for the complete soil column
 *          for each grid cell.
 *****************************************************************************/
typedef struct {
    double aero_resist[2];             /**< The (stability-corrected) aerodynamic
                                          resistance (s/m) that was actually used
                                          in flux calculations.
                                          [0] = surface (bare soil, non-overstory veg, or snow pack)
                                          [1] = overstory */
    double asat;                       /**< saturated area fraction */
    double baseflow;                   /**< baseflow from current cell (mm/TS) */
    double CLitter;                    /**< carbon storage in litter pool [gC/m2] */
    double CInter;                     /**< carbon storage in intermediate pool [gC/m2] */
    double CSlow;                      /**< carbon storage in slow pool [gC/m2] */
    double inflow;                     /**< moisture that reaches the top of
                                          the soil column (mm) */
    double pot_evap[N_PET_TYPES];      /**< array of different types of potential evaporation (mm) */
    double runoff;                     /**< runoff from current cell (mm/TS) */
    layer_data_struct layer[MAX_LAYERS]; /**< structure containing soil variables
                                            for each layer (see above) */
    double RhLitter;                   /**< soil respiration from litter pool [gC/m2] */
    double RhLitter2Atm;               /**< soil respiration from litter pool [gC/m2] that goes to atmosphere */
    double RhInter;                    /**< soil respiration from intermediate pool [gC/m2] */
    double RhSlow;                     /**< soil respiration from slow pool [gC/m2] */
    double RhTot;                      /**< total soil respiration over all pools [gC/m2] (=RhLitter2Atm+RhInter+RhSlow) */
    double rootmoist;                  /**< total of layer.moist over all layers
                                          in the root zone (mm) */
    double wetness;                    /**< average of
                                          (layer.moist - Wpwp)/(porosity*depth - Wpwp)
                                          over all layers (fraction) */
    double zwt;                        /**< average water table position [cm] - using lowest unsaturated layer */
    double zwt_lumped;                 /**< average water table position [cm] - lumping all layers' moisture together */
} cell_data_struct;

/******************************************************************************
 * @brief   This structure stores energy balance components, and variables used
 *          to solve the thermal fluxes through the soil column.
 *****************************************************************************/
typedef struct {
    // State variables
    double AlbedoLake;           /**< albedo of lake surface (fract) */
    double AlbedoOver;           /**< albedo of intercepted snow (fract) */
    double AlbedoUnder;          /**< surface albedo (fraction) */
    double Cs[2];                /**< heat capacity for top two layers (J/m^3/K) */
    double Cs_node[MAX_NODES];   /**< heat capacity of the soil thermal nodes (J/m^3/K) */
    double fdepth[MAX_FRONTS];   /**< all simulated freezing front depths */
    bool frozen;                   /**< TRUE = frozen soil present */
    double ice[MAX_NODES];       /**< thermal node ice content */
    double kappa[2];             /**< soil thermal conductivity for top two layers (W/m/K) */
    double kappa_node[MAX_NODES]; /**< thermal conductivity of the soil thermal nodes (W/m/K) */
    double moist[MAX_NODES];     /**< thermal node moisture content */
    size_t Nfrost;               /**< number of simulated freezing fronts */
    size_t Nthaw;                /**< number of simulated thawing fronts */
    double T[MAX_NODES];         /**< thermal node temperatures (C) */
    bool T_fbflag[MAX_NODES];      /**< flag indicating if previous step's temperature was used */
    unsigned T_fbcount[MAX_NODES]; /**< running total number of times that previous step's temperature was used */
    int T1_index;                   /**< soil node at the bottom of the top layer */
    double Tcanopy;              /**< temperature of the canopy air */
    bool Tcanopy_fbflag;           /**< flag indicating if previous step's temperature was used */
    unsigned Tcanopy_fbcount;    /**< running total number of times that previous step's temperature was used */
    double tdepth[MAX_FRONTS];   /**< all simulated thawing front depths */
    double Tfoliage;             /**< temperature of the overstory vegetation */
    bool Tfoliage_fbflag;            /**< flag indicating if previous step's temperature was used */
    unsigned Tfoliage_fbcount;   /**< running total number of times that previous step's temperature was used */
    double Tsurf;                /**< temperature of the understory */
    bool Tsurf_fbflag;           /**< flag indicating if previous step's temperature was used */
    unsigned Tsurf_fbcount;      /**< running total number of times that previous step's temperature was used */
    double unfrozen;             /**< frozen layer water content that is unfrozen */
    // Fluxes
    double advected_sensible;    /**< net sensible heat flux advected to snowpack (Wm-2) */
    double advection;            /**< advective flux (Wm-2) */
    double AtmosError;
    double AtmosLatent;          /**< latent heat exchange with atmosphere */
    double AtmosLatentSub;       /**< latent sub heat exchange with atmosphere */
    double AtmosSensible;        /**< sensible heat exchange with atmosphere */
    double canopy_advection;     /**< advection heat flux from the canopy (W/m^2) */
    double canopy_latent;        /**< latent heat flux from the canopy (W/m^2) */
    double canopy_latent_sub;    /**< latent heat flux of sublimation from the canopy (W/m^2) */
    double canopy_refreeze;      /**< energy used to refreeze/melt canopy intercepted snow (W/m^2) */
    double canopy_sensible;      /**< sensible heat flux from canopy interception (W/m^2) */
    double deltaCC;              /**< change in snow heat storage (Wm-2) */
    double deltaH;               /**< change in soil heat storage (Wm-2) */
    double error;                /**< energy balance error (W/m^2) */
    double fusion;               /**< energy used to freeze/thaw soil water */
    double grnd_flux;            /**< ground heat flux (Wm-2) */
    double latent;               /**< net latent heat flux (Wm-2) */
    double latent_sub;           /**< net latent heat flux from snow (Wm-2) */
    double longwave;             /**< net longwave flux (Wm-2) */
    double LongOverIn;           /**< incoming longwave to overstory */
    double LongUnderIn;          /**< incoming longwave to understory */
    double LongUnderOut;         /**< outgoing longwave from understory */
    double melt_energy;          /**< energy used to reduce snow cover fraction (Wm-2) */
    double NetLongAtmos;         /**< net longwave radiation to the atmosphere (W/m^2) */
    double NetLongOver;          /**< net longwave radiation from the overstory (W/m^2) */
    double NetLongUnder;         /**< net longwave radiation from the understory (W/m^2) */
    double NetShortAtmos;        /**< net shortwave to the atmosphere */
    double NetShortGrnd;         /**< net shortwave penetrating snowpack */
    double NetShortOver;         /**< net shortwave radiation from the overstory (W/m^2) */
    double NetShortUnder;        /**< net shortwave radiation from the understory (W/m^2) */
    double out_long_canopy;      /**< outgoing longwave to canopy */
    double out_long_surface;     /**< outgoing longwave to surface */
    double refreeze_energy;      /**< energy used to refreeze the snowpack (Wm-2) */
    double sensible;             /**< net sensible heat flux (Wm-2) */
    double shortwave;            /**< net shortwave radiation (Wm-2) */
    double ShortOverIn;          /**< incoming shortwave to overstory */
    double ShortUnderIn;         /**< incoming shortwave to understory */
    double snow_flux;            /**< thermal flux through the snow pack (Wm-2) */
} energy_bal_struct;

/******************************************************************************
 * @brief   This structure stores vegetation variables for each vegetation type
 *          in a grid cell.
 *****************************************************************************/
typedef struct {
    double albedo;              /**< current vegetation albedo (fraction) */
    double canopyevap;          /**< evaporation from canopy (mm/TS) */
    double LAI;                 /**< current leaf area index (m2/m2) */
    double throughfall;         /**< water that reaches the ground through the canopy (mm/TS) */
    double vegcover;            /**< current fractional area of plants within veg tile (fraction) */
    double Wdew;                /**< dew trapped on vegetation (mm) */
    double Wdmax;               /**< current maximum dew holding capacity (mm) */
    double *NscaleFactor;       /**< array of per-layer nitrogen scaling factors */
    double *aPARLayer;          /**< array of per-layer absorbed PAR (mol(photons)/m2 leaf area s) */
    double *CiLayer;            /**< array of per-layer leaf-internal CO2 mixing ratio (mol CO2/mol air) */
    double *rsLayer;            /**< array of per-layer stomatal resistance (s/m) */
    double aPAR;                /**< whole-canopy absorbed PAR (mol(photons)/m2 leaf area s) */
    double Ci;                  /**< whole-canopy leaf-internal CO2 mixing ratio (mol CO2/mol air) */
    double rc;                  /**< whole-canopy stomatal resistance (s/m) */
    double NPPfactor;           /**< whole-canopy photosynthesis multiplier to account for inhibition separate from stomatal resistance */
    double GPP;                 /**< whole-canopy gross assimilation (photosynthesis) (umol(CO2)/m2s) */
    double Rphoto;              /**< whole-canopy photorespiration (umol(CO2)/m2s) */
    double Rdark;               /**< whole-canopy 'dark' respiration (umol(CO2)/m2s) */
    double Rmaint;              /**< plant maintenance respiration (= Rdark/FRLeaf) (umol(CO2)/m2s) */
    double Rgrowth;             /**< growth respiration ( = (GPP-Rmaint)*FRGrowth/(1+FRGrowth) ) (umol(CO2)/m2s) */
    double Raut;                /**< total plant respiration (= Rmaint + Rgrowth) (umol(CO2)/m2s) */
    double NPP;                 /**< net primary productivity (= GPP - Raut) (umol(CO2)/m2s) */
    double Litterfall;          /**< flux of carbon from living biomass to litter pool [gC/m2] */
    double AnnualNPP;           /**< running total annual NPP [gC/m2] */
    double AnnualNPPPrev;       /**< total annual NPP from previous year [gC/m2] */
} veg_var_struct;

/******************************************************************************
 * @brief   This structure stores snow pack variables needed to run the snow
 *          model.
 *****************************************************************************/
typedef struct {
    // State variables
    double albedo;          /**< snow surface albedo (fraction) */
    double canopy_albedo;   /**< albedo of the canopy (fract) */
    double coldcontent;     /**< cold content of snow pack */
    double coverage;        /**< fraction of snow band that is covered with snow */
    double density;         /**< snow density (kg/m^3) */
    double depth;           /**< snow depth (m) */
    unsigned last_snow;     /**< time steps since last snowfall */
    double max_snow_depth;  /**< last maximum snow depth - used to determine coverage
                               fraction during current melt period (m) */
    bool MELTING;           /**< flag indicating that snowpack melted
                               previously */
    double pack_temp;       /**< depth averaged temperature of the snowpack (C) */
    double pack_water;      /**< liquid water content of the snow pack (m) */
    bool snow;              /**< TRUE = snow, FALSE = no snow */
    double snow_canopy;     /**< amount of snow on canopy (m) */
    double store_coverage;  /**< stores coverage fraction covered by new snow (m) */
    bool store_snow;        /**< flag indicating whether or not new accumulation
                               is stored on top of an existing distribution */
    double store_swq;       /**< stores newly accumulated snow over an
                               established snowpack melt distribution (m) */
    double surf_temp;       /**< depth averaged temperature of the snow pack surface layer (C) */
    unsigned surf_temp_fbcount; /**< running total number of times that previous step's temperature was used */
    bool surf_temp_fbflag;    /**< flag indicating if previous step's temperature was used */
    double surf_water;      /**< liquid water content of the surface layer (m) */
    double swq;             /**< snow water equivalent of the entire pack (m) */
    double snow_distrib_slope; /**< current slope of uniform snow distribution (m/fract) */
    double tmp_int_storage; /**< temporary canopy storage, used in snow_canopy */
    // Fluxes
    double blowing_flux;    /**< depth of sublimation from blowing snow (m) */
    double canopy_vapor_flux; /**< depth of water evaporation, sublimation, or
                                 condensation from intercepted snow (m) */
    double mass_error;      /**< snow mass balance error */
    double melt;            /**< snowpack melt (mm) */
    double Qnet;            /**< Residual of energy balance at snowpack surface */
    double surface_flux;    /**< depth of sublimation from blowing snow (m) */
    double transport;       /**< flux of snow (potentially) transported from veg type */
    double vapor_flux;      /**< depth of water evaporation, sublimation, or
                               condensation from snow pack (m) */
} snow_data_struct;

/******************************************************************************
 * @brief   This structure stores the lake/wetland parameters for a grid cell
 *****************************************************************************/
typedef struct {
    // Lake basin dimensions
    size_t numnod;                /**< Maximum number of lake nodes for this grid cell */
    double z[MAX_LAKE_NODES + 1]; /**< Elevation of each lake node (when lake storage is at maximum), relative to lake's deepest point (m) */
    double basin[MAX_LAKE_NODES + 1]; /**< Area of lake basin at each lake node (when lake storage is at maximum) (m^2) */
    double Cl[MAX_LAKE_NODES + 1]; /**< Fractional coverage of lake basin at each node (when lake storage is at maximum) (fraction of grid cell area) */
    double b;                     /**< Exponent in default lake depth-area profile (y=Ax^b) */
    double maxdepth;              /**< Maximum allowable depth of liquid portion of lake (m) */
    double mindepth;              /**< Minimum allowable depth of liquid portion of lake (m) */
    double maxvolume;             /**< Lake volume when lake depth is at maximum (m^3) */
    double minvolume;             /**< Lake volume when lake depth is at minimum (m^3) */
    // Hydrological properties
    double bpercent;              /**< Fraction of wetland baseflow (subsurface runoff) that flows into lake */
    double rpercent;              /**< Fraction of wetland surface runoff that flows into lake */
    double wfrac;                 /**< Width of lake outlet, expressed as fraction of lake perimeter */
    // Initial conditions
    double depth_in;              /**< Initial lake depth (distance from surface to deepest point) (m) */
    int lake_idx;                 /**< index number of the lake/wetland veg tile */
} lake_con_struct;

/******************************************************************************
 * @brief   This structure stores the lake/wetland variables for a grid cell
 *****************************************************************************/
typedef struct {
    // Current lake dimensions and liquid water state variables
    unsigned short activenod;     /**< Number of nodes whose corresponding layers currently contain water */
    double dz;                    /**< Vertical thickness of all horizontal water layers below the surface layer (m) */
    double surfdz;                /**< Vertical thickness of surface (top) water layer (m) */
    double ldepth;                /**< Current depth of liquid water in lake (distance from surface to deepest point) (m) */
    double surface[MAX_LAKE_NODES + 1]; /**< Area of horizontal cross-section of liquid water in lake at each node (m^2) */
    double sarea;                 /**< Current surface area of ice+liquid water on lake surface (m^2) */
    double sarea_save;            /**< Surface area of ice+liquid water on lake surface (m^2) at beginning of time step */
    double volume;                /**< Current lake water volume, including liquid water equivalent of lake ice (m^3) */
    double volume_save;           /**< Lake water volume, including liquid water equivalent of lake ice (m^3) at beginning of time step */
    double temp[MAX_LAKE_NODES];  /**< Lake water temperature at each node (C) */
    double tempavg;               /**< Average liquid water temperature of entire lake (C) */
    // Current properties (state variables) specific to lake ice/snow
    double areai;                 /**< Area of ice coverage (at beginning of time step) (m^2) */
    double new_ice_area;          /**< Area of ice coverage (at end of time step) (m^2) */
    double ice_water_eq;          /**< Liquid water equivalent volume of lake ice (m^3) */
    double hice;                  /**< Height of lake ice at thickest point (m) */
    double tempi;                 /**< Lake ice temperature (C) */
    double swe;                   /**< Water equivalence of lake snow cover - end of step (m^3) */
    double swe_save;              /**< Water equivalence of lake snow cover - beginning of step (m^3) */
    double surf_temp;             /**< Temperature of surface snow layer (C) */
    double pack_temp;             /**< Temperature of pack snow layer (C) */
    double coldcontent;           /**< cold content of snow pack */
    double surf_water;            /**< Water content of surface snow layer (m^3) */
    double pack_water;            /**< Water content of pack snow layer (m^3) */
    double SAlbedo;               /**< Albedo of lake snow (fraction) */
    double sdepth;                /**< Depth of snow on top of ice (m^3) */
    // Other current lake properties (derived from state variables and forcings)
    double aero_resist;           /**< Aerodynamic resistance (s/m) after stability correction */
    double density[MAX_LAKE_NODES]; /**< Lake water density profile (kg/m^3) */
    // Moisture fluxes
    double baseflow_in;           /**< Volume of baseflow into lake from the rest of the grid cell (m3) */
    double baseflow_out;          /**< Volume of baseflow out of lake to channel network (m3) */
    double channel_in;            /**< Volume of channel inflow into lake (m3) */
    double evapw;                 /**< Volume of evaporative flux from lake (and ice/snow) surface (m3) */
    double ice_throughfall;       /**< Volume of precipitation reaching lake water surface, i.e. total precip minus accumulation of snow on ice (m3) */
    double prec;                  /**< Volume of precipitation falling on lake (and ice/snow) surface (m3) */
    double recharge;              /**< Volume of recharge from lake to wetland (m3) */
    double runoff_in;             /**< Volume of surface runoff into lake from the rest of the grid cell (m3) */
    double runoff_out;            /**< Volume of surface runoff out of lake to channel network (m3) */
    double snowmlt;               /**< Volume of moisture released by melting of lake snow (m3) */
    double vapor_flux;            /**< Volume of moisture sublimated from lake snow (m3) */
    // Structures compatible with other land cover types
    // Some of this information is currently redundant with other variables in the lake_var structure
    snow_data_struct snow;        /**< Snow pack on top of lake ice */
    energy_bal_struct energy;     /**< Energy fluxes and soil temperatures */
    cell_data_struct soil;        /**< Soil column below lake */
} lake_var_struct;

/******************************************************************************
 * @brief   This structure stores all variables needed to solve, or save
 *          solututions for all versions of this model.
 *****************************************************************************/
typedef struct {
    cell_data_struct **cell;      /**< Stores soil layer variables */
    energy_bal_struct **energy;   /**< Stores energy balance variables */
    lake_var_struct lake_var;     /**< Stores lake/wetland variables */
    snow_data_struct **snow;      /**< Stores snow variables */
    veg_var_struct **veg_var;     /**< Stores vegetation variables */
} all_vars_struct;

/******************************************************************************
 * @brief   This structure stores moisture state information for differencing
 *          with next time step.
 *****************************************************************************/
typedef struct {
    double total_soil_moist;      /**< total column soil moisture [mm] */
    double surfstor;              /**< surface water storage [mm] */
    double swe;                   /**< snow water equivalent [mm] */
    double wdew;                  /**< canopy interception [mm] */
} save_data_struct;

/******************************************************************************
 * @brief   This structure stores output information for one variable.
 *****************************************************************************/
typedef struct {
    char varname[20];        /**< name of variable */
    bool write;              /**< FALSE = don't write; TRUE = write */
    char format[10];         /**< format, when written to an ascii file;
                                should match the desired fprintf format specifier, e.g. %.4f */
    unsigned short type;     /**< type, when written to a binary file;
                                OUT_TYPE_USunsigned short  = unsigned short int
                                OUT_TYPE_SINT   = short int
                                OUT_TYPE_FLOAT  = single precision floating point
                                OUT_TYPE_DOUBLE = double precision floating point */
    double mult;             /**< multiplier, when written to a binary file */
    unsigned short aggtype;  /**< type of aggregation to use;
                                AGG_TYPE_AVG    = take average value over agg interval
                                AGG_TYPE_BEG    = take value at beginning of agg interval
                                AGG_TYPE_END    = take value at end of agg interval
                                AGG_TYPE_MAX    = take maximum value over agg interval
                                AGG_TYPE_MIN    = take minimum value over agg interval
                                AGG_TYPE_SUM    = take sum over agg interval */
    unsigned nelem;          /**< number of data values */
    double *data;            /**< array of data values */
    double *aggdata;         /**< array of aggregated data values */
} out_data_struct;

/******************************************************************************
 * @brief   This structure stores output information for one output file.
 *****************************************************************************/
typedef struct {
    char prefix[20];         /**< prefix of the file name, e.g. "fluxes" */
    char filename[MAXSTRING];        /**< complete file name */
    FILE *fh;                /**< filehandle */
    size_t nvars;            /**< number of variables to store in the file */
    unsigned *varid;         /**< id numbers of the variables to store in the file
                                (a variable's id number is its index in the out_data array).
                                The order of the id numbers in the varid array
                                is the order in which the variables will be written. */
} out_data_file_struct;

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
    out_data_struct *out_data;
    out_data_file_struct *out_data_files;
    snow_data_struct *snow;
    soil_con_struct soil_con;
    veg_con_struct *veg_con;
    veg_var_struct *veg_var;
} Error_struct;

#endif
