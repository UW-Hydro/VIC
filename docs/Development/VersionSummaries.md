# Summary of VIC Versions

The tables below compare features (and bug fixes) of current and previous versions of VIC.

* * *

## VIC 4.2 (hotfix 4.2.d)

* **Partial Vegetation Cover**
    *   Accounts for gaps in vegetative cover, through which bare soil evaporation can occur In terms of radiation attenuation and turbulence, the bare soil between plants is treated like understory, rather than a large open area of bare soil.
    *   The vegetation and the bare soil in the gaps all share the same soil moisture reservoirs.
    *   The new vegetation parameter, "VEGCOVER", tells VIC the fractional area of vegetation; the remainder is treated as bare soil.
    *   12 climatological monthly values of VEGCOVER can be optionally specified in the vegetation library file in between the LAI and albedo values (if so, VEGLIB_VEGCOVER must be set to TRUE; to use these, VEGCOVER_SRC must be set to FROM_VEGLIB).
    *   Grid-cell-specific climatological monthly values of VEGCOVER can be optionally specified in the vegetation parameter file (if so, VEGPARAM_VEGCOVER must be set to TRUE; to use these; VEGCOVER_SRC must be set to FROM_VEGPARAM).
    *   A daily timeseries of VEGCOVER values can be given to VIC in its input forcing files, one column per veg tile (VEGCOVER must be specified as one of the forcing variables in the global parameter file); this will override the climatological values in the veg library and veg param files (see Non-Climatological Time-Varying Vegetation Parameters).

* **Non-Climatological Time-Varying Vegetation Parameters**
    *   Daily timeseries of LAI, VEGCOVER, and ALBEDO values can be given to VIC in its input forcing files (for LAI, the name of the input variable is LAI_IN).
    *   For each of these variables supplied in the forcing files, there must be one column per veg tile (e.g., for a cell with 3 veg tiles, there will need to be 3 columns of LAI, 3 columns of VEGCOVER, and/or 3 columns of ALBEDO).
    *   These daily values will override the climatological values in the veg library and veg param files.
    *   For any of these variables that are not supplied in a forcing file, the 12 climatological monthly values will be used instead; these will come from either the vegetation library file or the vegetation parameter file, depending on the value of LAI_SRC, VEGCOVER_SRC, or ALBEDO_SRC, respectively.
    *   New output variables: OUT_LAI, OUT_VEGCOVER (OUT_ALBEDO already existed).

* **Simulation of Carbon Cycle**
    *   Simulates photosynthesis, autotrophic respiration, heterotrophic (soil) respiration, and soil carbon storage.
    *   Phenology is not dynamic; LAI is still prescribed as before.
    *   Non-leaf vegetative biomass is not simulated; it is assumed to be constant in time.
    *   3 soil carbon pools: litter, intermediate, and slow.
    *   Previous year's NPP is transferred to litter pool (as litterfall) at a constant rate over the current year.
    *   New input variables: CATM (atmospheric carbon dioxide concentration), FDIR (fraction of solar radiation that is direct), PAR (photosynthetically active radiation); all of these are optional, with VIC using default values if they are not supplied New output variables: OUT_CATM, OUT_FDIR, OUT_PAR, OUT_GPP (gross primary productivity), OUT_RAUT (autotrophic respiration), OUT_NPP (GPP-RAUT), OUT_RHET (heterotrophic respiration), OUT_NEE (NPP-RHET), OUT_APAR (absorbed PAR), OUT_LITTERFALL (flux of C from live vegetation to litter pool), OUT_CLITTER (C density in litter pool), OUT_CINTER (C density in intermediate pool), and OUT_CSLOW (C density in slow pool).

* **Soil Moisture Content for Half-Space below Bottom Soil Layer**
    *   Previously, moisture content in the soil thermal nodes below the bottom of VIC's hydrologic soil column was set to the moisture content at the bottom of VIC's hydrologic soil column; but since VIC's soil temperature profile can extend far below the hydrologic layers, this was unrealistic.
    *   Now, VIC assumes a constant soil moisture content below the bottom of VIC's hydrologic soil column; this constant amount is stored in the variable SLAB_MOIST_FRACT (units = fraction of maximum content), which is defined in vicNl_def.h

* **Better Out-of-the-Box Behavior for Soil Temperature Scheme**
    *   Set default values of IMPLICIT and EXP_TRANS to TRUE.
    *   Made "cold" (no-spinup) initial soil temperatures more consistent with air temperature and bottom boundary temperature.
    *   Added validation of option.Nnodes for EXP_TRANS=TRUE to guarantee that, for the given soil temperature bottom boundary depth "dp" (also known as the damping depth), there are at least 3 nodes within the top 50 cm of the soil column. This is to constrain errors to a reasonable size. To satisfy this condition, the following relationship must hold:

        `Nnodes >= 5*ln(dp+1)+1`

* **Removed or Depreciated Options**
    *   Removed the DIST_PRCP option
    *   Removed user_def.h
    *   Removed ARC_SOIL and SOIL_DIR options from global.param.sample
    *   Removed LOG_MATRIC option
    *   Removed EXCESS_ICE option
    *   Removed SUN1999 albedo option
    *   Require JULY_TAVG_SUPPLIED=TRUE when COMPUTE_TREELINE=TRUE

## Bug Fixes in VIC 4.2.a-4.2.d

* **Fixed uninitialized dmy struct when writing binary output when OUTPUT_FORCE is True**
    *   Previously, the `dmy_struct` was not allocated or initializaed when `OUTPUT_FORCE == TRUE` and the output format is binary.  This fix corrects this bug.

* **Fixed uninitialized vegetation albedo when VEGPARAM_LAI is False**
    *   Previously, some vegetation parameters were left uninitialized when `VEGPARAM_LAI == FALSE`.  This bug has been corrected.

* **Fixed uninitialized bare soil albedo**
    *   Previously, bare_albedo was unset for the bare soil case (`iveg!=Nveg`). This fix sets the bare_albedo to the global variable value of `BARE_SOIL_ALBEDO`.

* **Cleanup of frozen soil option constraints**
    *   Removed hardcoded, behind the scenes checks for the `EXP_TRANS` and `NO_FLUX` global parameter values for case of `QUICK_SOLVE=TRUE` in `calc_surf_energy_bal`.

* **Fixed memory error in `initialize_atmos` when OUTPUT_FORCE = TRUE ***
    *   Previously, access to uninitialized elements of the veg_con and veg_hist structure was attempted when OUTPUT_FORCE = TRUE, causing a memory error and the model to crash.  This fix sets these elements inside a `if (!options.OUTPUT_FORCE)` block allowing the OUTPUT_FORCE option to work as expected.

* **Documented how VIC 4.2 needs user to specify veg_lib and veg_param files when OUTPUT_FORCE = TRUE**
    *   Prior to release 4.2, a user could run VIC in OUTPUT_FORCE mode with only a soil parameter file and forcing files.  This functionality is now broken as of release 4.2 and will not be fixed.  Users must either supply veg_lib and veg_parameter files (which the user is likely to have anyway) or use the standalone forcing disaggregator under cevelopment for use with release 5.0.  The documentation was updated to describe this issue as of release 4.2.c.

* **Added architectural resistance of 100 s/m for soil evaporation**
    *   Testing at approx. 60 eddy covariance towers ([Bohn and Vivoni, 2016](../Documentation/References.md#other-historical-references)) has indicated that soil evaporation is too high with the prior architectural resistance of 0 s/m and too low with a value of 200 s/m.  Further refinement would be ideal but this is a good ballpark figure.

* **Compute aerodynamic conductance of each veg tile as area-weighted average of conductances of vegetated and exposed soil fractions of the tile**
    *   The prior formulation was not the final version used in [Bohn and Vivoni (2016)](../Documentation/References.md#other-historical-references), but was mistakenly added to the codebase instead of the formulation used here.  This fixes the mistake.

* **Fix overwriting of veg_lib structure with values of current cell in veg_param file**
    *   Previously, VIC overwrote the LAI, albedo, and vegcover values in the copy of the veg library stored in memory (which is supposed to be constant reference values that apply to all grid cells) with those from the veg_parameter file pertaining to the current grid cell.  Values for veg classes not present in the current grid cell therefore were those of the last grid cell that contained those veg classes.  This did not affect performance but interfered with diagnostics while debugging.

* **Lake parameter validation**
    *   Previously, there were minimal checks performed on the values of the depth-area relationship.  This allowed unphysical values to be specified, leading to all manner of unphysical behaviors.  This has been fixed.

* **Fix lake water balance errors**
    *   Previously, precipitation over the lake was scaled by the lake area fraction twice, resulting in water balance errors.  This has been fixed.


## VIC 4.1.2

* **Global Parameter File**
    *   Re-organized the sample global parameter file that comes with the code distribution (global.param.sample). The large number of options available in VIC 4.1.2 has begun to clutter the global parameter file, making it difficult for users to determine which options need to be set. Now the sample global parameter file groups the various options into several categories, making a clear distinction between those options that need to be changed by the user for each simulation, and those that rarely need to be changed. Those options that rarely need to be changed are commented out, and VIC now sets those options to the most commonly-used default values.
        *   Removed the **AR_COMBO** value for **AERO_RESIST_CANSNOW** and **GF_FULL** value for **GRND_FLUX_TYPE**
        *   Removed all output options related to **LINK_DEBUG**

* **Forcing estimation/ disaggregation**
    *   We have improved some of the algorithms VIC uses for estimating daily and sub-daily short- and long-wave radiaton and humidity (see [Bohn et al., 2013a](../Documentation/References.md) for more details)
    *   To access these improvements, we have added the following new global parameter file options:
        *   **MTCLIM_SWE_CORR**: optional correction of downward shortwave for the effect of snow, taken from MTCLIM 4.3 (Thornton et al., 2000). Default value is TRUE.
        *   **VP_ITER**: optional changes in the iteration between shortwave and VP estimates. Allowed values (**some of which may be removed before release**) are:
            *   **VP_ITER_NEVER** = never perform the final update of Tdew and VP. (our testing indicates that this option generally performs worse than VP_ITER_ALWAYS).
            *   **VP_ITER_ALWAYS** = always perform the final update of Tdew and VP. This is the traditional behavior in previous VIC releases.
            *   **VP_ITER_ANNUAL** = use an annual aridity criterion (annual PET < 2.5 * annual precip) do decide whether to perform the final update. Taken from MTCLIM 4.3 (Thornton et al., 2000). (our testing indicates that this option generally performs worse than VP_ITER_ALWAYS).
            *   **VP_ITER_CONVERGE** = continue iteratively updating both shortwave and Tdew/VP until all values stabilize within a given tolerance. (our testing indicates that this option yields negligible differences from VP_ITER_ALWAYS).
            *   Default value is **VP_ITER_ALWAYS**.
        *   **VP_INTERP**: optional linear interpolation between the daily VP estimates to get a sub-daily varying VP (as opposed to holding VP constant throughout the day, as in previous versions of VIC). Default value is TRUE.
        *   **LW_TYPE**: optional alternative clear-sky longwave radiation algorithms. Allowed values (**some of which may be removed before release**) are:
            *   **LW_TVA** = algorithm of Tennessee Valley Authority (TVA, 1972) (This is what all previous versions of VIC have used. Our tests indicate that this algorithm is still the best when observed cloud fractions are unavailable and are estimated by MTCLIM, which is the current situation for VIC.)
            *   **LW_ANDERSON** = algorithm of Anderson (1964)
            *   **LW_BRUTSAERT** = algorithm of Brutseart (1975)
            *   **LW_SATTERLUND** = algorithm of Satterlund (1979)
            *   **LW_IDSO** = algorithm of Idso (1981)
            *   **LW_PRATA** = algorithm of Prata (1996) (Our tests indicate that this algorithm is best when cloud fractions are supplied as a forcing.)
            *   Default is set to **LW_TVA**.
        *   **LW_CLOUD**: optional alternative cloud longwave radiation algorithms. Allowed values are:
            *   **LW_CLOUD_BRAS** = algorithm composed of equations 2.29 and 2.43 from the Bras textbook (Bras, R. L., "Hydrology, an introduction to hydrologic science", Addison-Wesley, Reading, Massachusetts, 1990). This was the algorithm used in all previous releases of VIC. (Our tests indicate that this algorithm introduces substantial temperature-dependent biases in longwave estimates outside of the temperate zone)
            *   **LW_CLOUD_DEARDORFF** = algorithm used in Deardorff (1978) in which cloud_fraction is assumed equal to (1 - actual_shortwave/theoretical_clear_sky_shortwave) and total sky emissivity is represented as the weighted average: [ cloud_fraction*cloud_emissivity + (1-cloud_fraction)*clear_sky_emissivity ] (Our tests indicate that this algorithm is superior)
            *   Default is set to **LW_CLOUD_DEARDORFF**.
        *   **SW_PREC_THRESH**: optional minimum daily precipitation [mm] that must be exceeded to cause a dimming of 25% in estimated daily incoming shortwave radiation. Previously, any precipitation would cause estimated daily incoming shortwave radiation to dim by 25%. Because the smoothing/resampling methods used in creating gridded forcings can artificially "bleed" trace amounts of precipitation into neighboring grid cells, we introduced this threshold to counteract any resulting low biases in shortwave. Default is set to 0.1 mm.

* **State Files**
    *   Tool to convert state files created by older versions of VIC to the format used by the latest version of VIC (VIC_4.1.2). This stand-alone perl script is now available on the VIC web site: [convert_state_file.pl](ftp://ftp.hydro.washington.edu/pub/HYDRO/models/VIC/Utility_Programs/convert_state_file.pl).

* **Features related to lake model**
    *   Added output variable OUT_LAKE_AREA_FRAC, the lake surface area as a fraction of the grid cell area.
    *   Added input/output variables CHANNEL_IN and OUT_LAKE_CHAN_IN, allowing input of channel flow from upstream grid cells into the lake model.
        *   New input forcing variable CHANNEL_IN allows VIC to read channel flow or runoff volume (m^3; or m^3/s if ALMA_INPUT=TRUE) from upstream grid cells and use this as an inflow into the lake model
        *   New output variable OUT_LAKE_CHAN_IN reports this incoming flow (as mm over grid cell; mm/s if ALMA_OUTPUT=TRUE)
    *   Added several other lake-related water budget terms to the list of output variables. Now, the complete water budget of the lake can be monitored as follows (in mm over grid cell area; or mm/s if ALMA_OUTPUT=TRUE):
        *   OUT_LAKE_BF_IN baseflow into lake from catchment
        *   OUT_LAKE_BF_OUT baseflow out of lake
        *   OUT_LAKE_CHAN_IN channel inflow into lake
        *   OUT_LAKE_CHAN_OUT channel outflow from lake
        *   OUT_LAKE_DSTOR change in lake water storage (liquid plus ice cover)
        *   OUT_LAKE_DSWE change in snowpack on top of lake ice
        *   OUT_LAKE_EVAP (was OUT_EVAP_LAKE) evaporation at lake surface
        *   OUT_PREC (already existed) precipitation over lake (== precip over grid cell)
        *   OUT_LAKE_RCHRG recharge from lake to wetland
        *   OUT_LAKE_RO_IN runoff into lake from catchment
        *   OUT_LAKE_SWE water equivalent of snowpack on top of lake ice; OUT_LAKE_DSWE is the derivative of this
        *   OUT_LAKE_VAPFLX sublimation from lake snow
    *   To monitor the lake water budget in terms of volume (m^3; or m^3/s if ALMA_OUTPUT=TRUE):
        *   OUT_LAKE_BF_IN_V
        *   OUT_LAKE_BF_OUT_V
        *   OUT_LAKE_CHAN_IN_V
        *   OUT_LAKE_CHAN_OUT_V
        *   OUT_LAKE_DSTOR_V
        *   OUT_LAKE_DSWE_V
        *   OUT_LAKE_EVAP_V
        *   OUT_LAKE_PREC_V note: this is the volume of precip falling on lake only
        *   OUT_LAKE_RCHRG_V
        *   OUT_LAKE_RO_IN_V
        *   OUT_LAKE_VAPFLX_V
        *   OUT_LAKE_VOLUME water stored in lake (liquid plus ice cover); OUT_LAKE_DSTOR_V is the derivative of this
    *   see [Gao et al., 2011](../Documentation/References.md) for an example application using CHANNEL_IN to convey runoff from upstream into a lake

* **Water table position**
    *   Added computation of water table position [cm]. Some details:
        *   We use the formulation described in Frolking et al. (2002), and elements of the CLASS formulation of Letts et al. (2000). The water table's position within a given soil layer is computed using a Van Genuchten soil water retention relationship, which depends on the soil layer's porosity, bubble pressure, Brooks-Corey exponent, etc.
        *   Water table position is 0 at the soil surface, and negative below the soil surface. It should never be a positive number.
        *   The water table position within the entire soil column is computed in 2 different ways, which can be saved as output via the new variables OUT_ZWT and OUT_ZWT_LUMPED:
            *   OUT_ZWT: Water table position falls in the highest layer having complete saturation in all layers below it (this is usually the lowest layer, since the lowest layer is rarely completely saturated).
            *   OUT_ZWT_LUMPED: Water table position is calculated by taking all soil moisture from the entire soil column and using it to fill the soil layers from the bottom layer upwards, as if all the water in the soil column were to sink downwards as far as it can go. The reason for this method is to eliminate jumps in the water table position as the soil moisture in the top layer passes the minimum threshold required to have a water table. This seems somewhat artificial at first glance, but the shape of the timeseries of soil moisture given by this method looks more realistic. It simply may be biased in certain circumstances.
    *   see [Bohn et al., 2012a](../Documentation/References.md) for more details

* **Soil temperatures**
    *   Extended the computation of soil temperatures, ice contents, and ground fluxes to all modes of model operation.
        *   Previously, soil temperatures (at the nodes in the vertical temperature profile) only explicitly computed and output for FROZEN_SOIL TRUE and QUICK_FLUX FALSE; average temperatures of the hydrologic soil layers were never output. Now all of these are computed and output for all settings of FROZEN_SOIL and QUICK_FLUX.
        *   Ice contents will always be 0 except when FROZEN_SOIL is TRUE.
        *   The default method of computing soil temperature profiles depends on the settings of FROZEN_SOIL and QUICK_FLUX:
            *   FROZEN_SOIL = TRUE: default method is the finite element method of [_Cherkauer and Lettenmaier_ (1999)](../Documentation/References.md). This corresponds to QUICK_FLUX=FALSE.
            *   FROZEN_SOIL = FALSE: default method is the approximation of [_Liang et al._ (1999)](../Documentation/References.md#vic). This corresponds to QUICK_FLUX=TRUE.
        *   To change which method is used to compute soil temperatures, set QUICK_FLUX appropriately in the global parameter file. (QUICK_FLUX TRUE = method of [_Liang et al._ (1999)](../Documentation/References.md#vic); FALSE = method of [_Cherkauer and Lettenmaier_ (1999)](../Documentation/References.md))

* **Organic Soil**
    *   VIC now accounts for the effects of organic matter on soil properties.
        *   VIC uses the equations of Farouki (1981) to compute thermal conductivity and heat capacity of soil as a function of organic content.
        *   To use this feature, you must specify ORGANIC_FRACT TRUE in your global parameter file.
        *   If ORGANIC_FRACT is TRUE, VIC assumes that there are 3*Nlayer extra columns in your soil parameter file (immediately following the soil particle density). These extra columns are: 1. Organic fraction of each soil layer, 2. bulk density for each soil layer, and 3. soil particle density for each soil layer.
        *   In this case, VIC assumes that the original quartz content, bulk density, and soil density columns in the soil parameter file pertain to the mineral fraction of the soil. Internally, VIC computes aggregate versions of soil properties (including porosity, thermal conductivity, and heat capacity) as a weighted average of the mineral and organic versions of these properties. The weights are (1-organic_fract)*mineral + organic_fract*organic.
        *   If ORGANIC_FRACT is FALSE, or is omitted entirely from the global parameter file, VIC behaves exactly as before (no extra fields in the soil param file) and organic fraction is assumed to be 0 everywhere.
    *   see [Bohn et al., 2012a](../Documentation/References.md) for more details

* **Other**
    *   Removed **LINK_DEBUG** option, which was rarely used and difficult to maintain in the code base.
    *   Removed the **AR_COMBO** value of the **AERO_RESIST_CANSNOW** option, since it was double-counting aerodynamic resistances.
    *   Removed the **GF_FULL** value of the **GRND_FLUX_TYPE** option, since it was double-counting the attenuation of radiation by the canopy.
    *   Removed **MIN_LIQ** option, which was confusing, unnecessary, and incorrect.

* * *

### Bug Fixes in VIC 4.1.2

* **Forcing estimation/ disaggregation**
    *   Fixed bug (present in all previous releases) in which the annual cycle of daylight hours shifted back by one day in 3 out of every 4 years (i.e., non-leap years). This affected the diurnal cycle of radiation and temperature, most noticeably at high latitudes (north of 66 N, i.e. Arctic Circle, where some days of the year have 0 hours of daylight). For long simulations (e.g., > 70 years), winter day lengths (e.g., a few or even zero hours of daylight) ended up being applied to spring days, possibly interfering with snow melt. Effects south of 66 N were minimal.
    *   Fixed bug (present in all previous releases) in which the diurnal cycle of temperature would occasionally skip a day's tmin and linearly interpolate between one day's tmax and the next day's tmax. This only occurred at high latitudes (North of 60 N) on days having exactly one hour of darkness.
    *   Fixed inconsistencies (present in all previous releases) in timing of sub-daily forcings computed by VIC and sub-daily forcings supplied by the user. Now, we have adopted the following convention:
        *   The parameter **off_gmt** (from the [soil parameter file](../Documentation/SoilParam.md)) determines how VIC interprets sub-daily time steps relative to the model start date and time.
        *   An **off_gmt** value of **0** indicates that the model start date/time is relative to **Greenwich Mean Time (GMT)**.
        *   An **off_gmt** value of **(grid_cell_longitude*24/360)** indicates that the model start date/time is relative to **local time**.
        *   When outputting sub-daily results, VIC's output files are referenced to the model start date/time; therefore they are controlled by **off_gmt** (off_gmt=0 means VIC results are referenced to GMT; off_gmt=(grid_cell_longitude*24/360) means VIC results are referenced to local time).
        *   **Daily** supplied forcings are assumed to start/end at midnight in **local time**; the forcing start date/time is thus in **local time**. When VIC disaggregates daily forcings into sub-daily forcings, **off_gmt** will be used to determine the time lag between the start of the forcing's diurnal cycle and the start of the VIC simulation.
        *   **Sub-daily** supplied forcings are assumed to occur relative to **the time zone indicated by off_gmt**. Therefore, if VIC outputs these sub-daily forcings, they will occur at the exact same time of day as in the input files.
    *   Therefore, if mixing daily and sub-daily forcing inputs, it is important that any sub-daily forcing inputs be shifted as necessary to be in the time zone indicated by off_gmt.

* **Soil moisture and runoff**
    *   Fixed bug (present in 4.1.1 only) in which soil moisture was allowed to fall below residual moisture level.
        *   NOTE: 4.1.1 appears to replace residual moisture with 0 under most modes of operation. The most noticeable consequence of this behavior is that baseflow will exhibit less variability (lower peak values and greater low flows) and soil moisture will be lower than it should be, for a given set of soil parameters. This appears to be the result of the MIN_LIQ option. The MIN_LIQ option was removed in version 4.1.2, so that baseflow and soil moisture behave correctly. But this means that soil parameters calibrated for 4.1.1 may need modification to yield similar results with 4.1.2. A simple fix to force 4.1.2 to behave like 4.1.1 is to set residual moisture to 0 everywhere in the soil parameter file.
    *   Fixed bug (present in all previous releases) in runoff computation when soil column is completely saturated.

* **Errors in lake model**
    *   (present in 4.1.1)
    *   Water balance errors in the lake model when the lake changes area, in particular when ice fraction > 0.
    *   Garbage output when lake fraction goes to 0.
    *   Fixed typo in writing of lake state information in VIC state file.
    *   Fixed typo in writing of number of lake active nodes to ASCII-format state file.
    *   Fixed VIC's inability to handle the case in which some cells have lakes and others don't. Now users can specify cells that have no lakes by setting lake_idx to -1 for those cells in the lake parameter file.

* **Soil Temperature Scheme (primarily for frozen soils)**
    *   Now users are prevented from inadvertently selecting an unstable soil temperature scheme (present in all previous releases)
        *   Now, whenever FROZEN_SOIL is set to TRUE, IMPLICIT is set to TRUE by default. The implicit scheme is stable for any combination of time step and thermal node spacing.
        *   Users may override this by setting IMPLICIT to FALSE in the global parameter file (thereby selecting an explicit scheme). However, VIC now performs a validity check at the beginning of the simulation to determine if the selected model time step and thermal node spacing will lead to instability in the solution; if so, VIC aborts with an error message. (It should be noted that this check is approximate, since stability depends on soil thermal properties, which change during the simulation as soil moisture and temperature change.)
    *   Added a crude HACK to prevent the formation of runaway "cold nose" in soil temperature profile (in presence of ice) (bug is present in all previous releases)
        *   This is only enabled when options.TFALLBACK=TRUE
        *   This hack basically recognizes when one node is colder than its upper and lower neighbors and continues to get colder (un-physical behavior) and sets the node's temperature to the average of its upper and lower neighbors. This results in a small energy balance error but prevents the runaway behavior we have seen in the past.
        *   Preliminary testing indicates that this replaces runaway cold noses with reasonably realistic behaviors, and does not change results for cases when cold noses do not occur.
        *   We hope to implement a more rigorous fix to the numerics of the soil thermal solution in future releases.
    *   Fixed slightly incorrect fallback temperature assignment for finite difference solution. (bug present in 4.1.1 only)
    *   Fixed incorrect calculation of grnd_flux, deltaH, and fusion for EXP_TRANS=TRUE. (bug present in 4.1.1 only)
    *   Added TFALLBACK logic to soil thermal profile solution for case in which max iterations are exceeded. (bug present in 4.1.1 only)
    *   Fixed initialization of various Tfallback counts and flags. (bug present in 4.1.1 only)

* **State files**
    *   Fixed typo in writing of lake state information in VIC state file. (bug present in 4.1.1 only)
    *   Fixed typo in writing of number of lake active nodes to ASCII-format state file. (bug present in 4.1.1 only)
    *   Fixed typo in condition for reading Wdew values. (bug present in all previous releases)

* **Screen output**
    *   Fixed the printing of cumulative water/energy balance errors, numbers of Temperature fallbacks, etc. to the screen; they were ending up out of order when saved to a log file. (bug present in 4.1.1 only)

* **Miscellaneous**
    *   Replaced "assert" statements with "if" statements. This fixes problems with sporadic aborting in some environments. (bug present in all previous releases)
    *   Miscellaneous fixes to how aerodynamic resistance and conductance are handled. (bug present in all previous releases)
    *   Fixed snow albedo decay curves to switch from accumulation to melt at the correct time of year in the southern hemisphere. (bug present in all previous releases)
    *   Fixed failure to handle 5-digit lat/lons correctly. (bug present in all previous releases)

* * *

## VIC 4.1.1

* **Potential ET**
    *   6 Types of PET can now be output:
        *   OUT_PET_SATSOIL: saturated soil
        *   OUT_PET_H2OSURF: water surface (smooth/shallow)
        *   OUT_PET_SHORT: short reference crop (potential transpiration only)
        *   OUT_PET_TALL: tall reference crop (potential transpiration only)
        *   OUT_PET_NATVEG: natural vegetation, no soil moisture limitation (potential transpiration only)
        *   OUT_PET_VEGNOCR: natural vegetation with Rc = 0 (potential transpiration only)

* **T instablilty handling**
    *   In past, soil T solution often failed to converge in arctic climates
    *   Jenny Adam's permafrost enhancements solved most problems
    *   Just in case, VIC now has better error handling: CONTINUEONERROR in global param file
        *   FALSE: VIC exits completely when a grid cell encounters T failure
        *   TRUE: VIC moves to next grid cell when T failure occurs
        *   TFALLBACK: use previous timestep's T value and continue with current grid cell
            *   Applies to surface and canopy T solutions as well
            *   NOTE: by the time the solution fails to converge, T values may have been unreasonable for days to weeks

* **Aerodynamic Resistance in Snow-filled Canopy**
    *   4.0.6 and 4.1.0 handled this issue differently, neither was correct
    *   New global parameter option: AERO_RESIST_CANSNOW
        *   AR_406: multiply Ra by 10 but apply to latent heat
        *   AR_406_LS: multiply by 10 and apply to both latent and sensible heat
        *   AR_406_FULL: same as AR_406_LS but also use canopy Ra (no correction) when there is no snow
        *   AR_410: apply stability correction to both latent and sensible heat
        *   AR_COMBO: apply stability correction and multipy by 10
    *   410 accumulates more canopy snow than 406 due to higher Ra; 411 similar to 406 but with lower sensible heat

* **Ground Flux**
    *   4.0.6 and 4.1.0 handled this differently; 4.1.0 was more correct, but not perfect
    *   New global parameter option: GRND_FLUX_TYPE
        *   GF_406: use (flawed) formulas for grnd_flux, deltaH, and fusion
        *   GF_410: use correct formula for grnd_flux (Liang et al. 1999) but deltaH and fusion don't take surf_atten into account
        *   GF_FULL: use correct formula for grnd_flux and take surf_atten in to account in deltaH and fusion
    *   410 and 411 have larger grnd_flux than 406

* **Snow Albedo**
    *   new global parameter file option: SNOW_ALBEDO
        *   USACE (default): previous algorithm, from US Army Corps; prescribed date of beginning of melt season (optimized for US)
        *   SUN1999: melt condition = function of cold content only (applicable globally)

* **Snow Density**
    *   New global parameter file option: SNOW_DENSITY
        *   DENS_BRAS (default): previous algorithm, taken from Bras (1992)
        *   DENS_SNTHRM: algorithm taken from SNTHRM model
    *   No computations depend on snow density, but it is good for calidation of model results against observations

* **Pressure Lapsing**
    *   Previously, when P and density computed during forcing disaggregation they were set to constant 400 meter elevation values
        *   P and density were lapsed to cell or band elevation for canopy evaporation and transpiration calculations
        *   But sensible heat used un-lapsed values
        *   Cold bias at high elevations, warm bias at low elevations
        *   Magnitude of effect unknown
    *   New global parameter option: PLAPSE
        *   FALSE: previous behavior
        *   TRUE (default): lapse pressure to grid cell average elevation

* **Changes in state file format**
    *  The state files have undergone some small changes in format, in response to bug fixes and new features. Among the changes: 0-area snow/elevation bands are once again stored in state files; some previously-missed state variables associated with lakes and wetlands are now stored. Therefore, VIC 4.1.1 will not be able to read state files saved from previous versions of VIC.

* **Change to value of OUT_LONGWAVE in output files for water balance mode**
    * Previously, in water balance mode, VIC would store the net longwave in the OUT_LONGWAVE variable, rather than incoming longwave. Now, VIC always stores incoming longwave in OUT_LONGWAVE, regardless of the setting of FULL_ENERGY.

* * *

### Differences Between VIC Versions 4.1.1 and 4.1.0

* 24 hour water balance mode
    *   Snow step was broken in 4.1.0 - now fixed

* Treeline Computation
    *   use in conjunction with snow/elevation bands
    *   Uses avg July T of forcings to determine which bands are above treeline
    *   Above treeline, replaces any overstory vegetation tiles with non-overstory vegetation cover
    *   Otherwise, runaway snowpack ("glaciers")
    *   Set COMPUTE_TREELINE to a veg. class ID number in global parameter file (ID number of non-overstory class)
    *   JulyAvgT can be supplies in soil file instead of computed from forcings:
        *   Set JULY_TAVG_SUPPLIED to TRUE in global parameter file

* * *

### Differences Between VIC Versions 4.1.1 and 4.0.6

* Dynamic Lake/Wetland Model
    *   Multi-layer lake model of Hostetler et al. 2000
        *   Energy-balance model
        *   Mixing, radiation attenuation, variable ice cover
    *   Dynamic lake area (taken from topography) allows seasonal inundation of adjacent wetlands
    *   Currently not a part of the channel network
* Permafrost Enhancements
    *   Global Parameter Options
        *   IMPLICIT = implicit solver
        *   EXP_TRANS = exponential node spacing
    *   User_def.h Options
        *   EXCESS_ICE = Excess Ground Ice and Subsidence Algorithm
            *   Excess ice is the concentration of ice in excess of what the soil can hold were it unfrozen
            *   As excess ice in a soil layer melts, the ground subsides
* Soil T heterogeneity: "Spatial Frost"
    *   Linear (uniform) distribution of soil T around mean
    *   Allows some moisture movement in soil when avg T below freezing
    *   Requires extra soil parameter: frost_slope
    *   User_def.h (compile time) option:
        *   Set SPATIAL_FROST to TRUE
        *   Set FROST_SUBAREAS (number of frost "bands")
* Partial Snow Cover: "Spatial Snow"
    *   If melting & avg SWE < depth_full_snow_cover, snow coverage is proportional to SWE
    *   Requires extra soil parameter: depth_full_snow_cover
    *   User_def.h (compile time) option:
        *   Set SPATIAL_SNOW to TRUE
* Canopy Temperatures and Energy Balance
    *   Canopy Surface Temperature
        *   When snow is in the canopy, solves for canopy snow temperature (Tfoilage) that balances canopy snow surface energy budget; otherwise surface temperature equals canopy air temperature (in 4.0.6 and earlier, canopy snow surface temperature did not balance the snow energy budget).
    *   Canopy Air Temperature
        *   By default equals air temperature (as in 4.0.6 and earlier)
        *   If snow is present in the canopy. and CLOSE_ENERGY = TRUE in user_def.h (compile-time option), VIC will determine canopy air temperature by solving the canopy surface energy budget and will iterate between canopy an land surface energy budget solutions until fluxes between the two systems are reconciled.
* Blowing Snow
    *   Blowing snow sublimation
    *   Set BLOWING to TRUE in global parameter file
