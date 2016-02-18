/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine initializes the lake variables for each new grid cell.
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

#include <vic_def.h>
#include <vic_run.h>

/******************************************************************************
 * @brief    This routine initializes the lake variables for each new grid cell
 *****************************************************************************/
int
initialize_lake(lake_var_struct  *lake,
                lake_con_struct   lake_con,
                soil_con_struct  *soil_con,
                cell_data_struct *cell,
                double            airtemp,
                int               skip_hydro)
{
    extern option_struct     options;
    extern parameters_struct param;

    size_t                   i, j;
    int                      k;
    int                      status;
    double                   depth;
    double                   tmp_volume;

    /*  Assume no ice present, lake completely equilibrated with atmosphere. */

    for (i = 0; i < MAX_LAKE_NODES; i++) {
        lake->temp[i] = max(airtemp, 0.0);
    }

    lake->areai = 0.0;
    lake->coldcontent = 0.0;
    lake->hice = 0.0;
    lake->ice_water_eq = 0.0;
    lake->new_ice_area = 0.0;
    lake->pack_temp = 0.0;
    lake->pack_water = 0.0;
    lake->SAlbedo = 0.0;
    lake->sdepth = 0.0;
    lake->surf_temp = 0.0;
    lake->surf_water = 0.0;
    lake->swe = 0.0;
    lake->swe_save = 0.0;
    lake->tempi = 0.0;

    /********************************************************************/
    /* Initialize lake physical parameters.                             */
    /********************************************************************/

    if (!skip_hydro) {
        if (lake_con.depth_in > lake_con.z[0]) {
            lake_con.depth_in = lake_con.z[0];
        }

        lake->ldepth = lake_con.depth_in;

        if (lake->ldepth > param.LAKE_MAX_SURFACE && lake->ldepth < 2 *
            param.LAKE_MAX_SURFACE) {
            /* Not quite enough for two full layers. */
            lake->surfdz = lake->ldepth / 2.;
            lake->dz = lake->ldepth / 2.;
            lake->activenod = 2;
        }
        else if (lake->ldepth >= 2 * param.LAKE_MAX_SURFACE) {
            /* More than two layers. */
            lake->surfdz = param.LAKE_MAX_SURFACE;
            lake->activenod = (int) (lake->ldepth / param.LAKE_MAX_SURFACE);
            if (lake->activenod > MAX_LAKE_NODES) {
                lake->activenod = MAX_LAKE_NODES;
            }
            lake->dz =
                (lake->ldepth -
                 lake->surfdz) / ((double) (lake->activenod - 1));
        }
        else if (lake->ldepth > 0.0) {
            lake->surfdz = lake->ldepth;
            lake->dz = 0.0;
            lake->activenod = 1;
        }
        else {
            lake->surfdz = 0.0;
            lake->dz = 0.0;
            lake->activenod = 0;
            lake->ldepth = 0.0;
        }

        // lake_con.basin equals the surface area at specific depths as input by
        // the user in the lake parameter file or calculated in read_lakeparam(),
        // lake->surface equals the area at the top of each dynamic solution layer

        for (k = 0; k <= lake->activenod; k++) {
            if (k == 0) {
                depth = lake->ldepth;
            }
            else {
                depth = lake->dz * (lake->activenod - k);
            }
            status = get_sarea(lake_con, depth, &(lake->surface[k]));
            if (status < 0) {
                log_err("Error in get_sarea: record = %d, depth = %f, "
                        "sarea = %e", 0, depth, lake->surface[k]);
            }
        }

        lake->sarea = lake->surface[0];
        lake->sarea_save = lake->sarea;
        status = get_volume(lake_con, lake->ldepth, &tmp_volume);
        if (status < 0) {
            log_err("Error in get_volume: record = %d, depth = %f, "
                    "volume = %e", 0, depth, tmp_volume);
            return(status);
        }
        else if (status > 0) {
            log_err("Warning in get_volume: lake depth exceeds maximum; "
                    "setting to maximum; record = %d", 0);
        }
        lake->volume = tmp_volume + lake->ice_water_eq;
        lake->volume_save = lake->volume;

        // Initialize lake moisture fluxes to 0
        lake->baseflow_in = 0.0;
        lake->baseflow_out = 0.0;
        lake->channel_in = 0.0;
        lake->evapw = 0.0;
        lake->prec = 0.0;
        lake->recharge = 0.0;
        lake->runoff_in = 0.0;
        lake->runoff_out = 0.0;
        lake->snowmlt = 0.0;
        lake->vapor_flux = 0.0;
    } // if (!skip_hydro)

    // Initialize other miscellaneous lake properties
    lake->aero_resist = 0;
    for (k = 0; k < lake->activenod; k++) {
        lake->density[k] = CONST_RHOFW;
    }

    // Initialize the snow, energy, and soil components of lake structure
    // If we implement heat flux between lake and underlying soil, we will need to initialize these more correctly
    // Snow state vars
    lake->snow.albedo = 0.0;
    lake->snow.canopy_albedo = 0.0;
    lake->snow.coldcontent = 0.0;
    lake->snow.coverage = 0.0;
    lake->snow.density = 0.0;
    lake->snow.depth = 0.0;
    lake->snow.last_snow = 0;
    lake->snow.max_snow_depth = 0.0;
    lake->snow.MELTING = false;
    lake->snow.pack_temp = 0.0;
    lake->snow.pack_water = 0.0;
    lake->snow.snow = false;
    lake->snow.snow_canopy = 0.0;
    lake->snow.store_coverage = 0.0;
    lake->snow.store_snow = false;
    lake->snow.store_swq = 0.0;
    lake->snow.surf_temp = 0.0;
    lake->snow.surf_temp_fbflag = 0;
    lake->snow.surf_temp_fbcount = 0;
    lake->snow.surf_water = 0.0;
    lake->snow.swq = 0.0;
    lake->snow.snow_distrib_slope = 0.0;
    lake->snow.tmp_int_storage = 0.0;
    // Snow fluxes
    lake->snow.blowing_flux = 0.0;
    lake->snow.canopy_vapor_flux = 0.0;
    lake->snow.mass_error = 0.0;
    lake->snow.melt = 0.0;
    lake->snow.Qnet = 0.0;
    lake->snow.surface_flux = 0.0;
    lake->snow.transport = 0.0;
    lake->snow.vapor_flux = 0.0;
    // Energy state vars
    lake->energy.AlbedoLake = 0.0;
    lake->energy.AlbedoOver = 0.0;
    lake->energy.AlbedoUnder = 0.0;
    lake->energy.frozen = 0.0;
    lake->energy.Nfrost = 0;
    lake->energy.Nthaw = 0;
    lake->energy.T1_index = 0;
    lake->energy.Tcanopy = 0.0;
    lake->energy.Tcanopy_fbflag = 0;
    lake->energy.Tcanopy_fbcount = 0;
    lake->energy.Tfoliage = 0.0;
    lake->energy.Tfoliage_fbflag = 0;
    lake->energy.Tfoliage_fbcount = 0;
    lake->energy.Tsurf = lake->temp[0];
    lake->energy.Tsurf_fbflag = 0;
    lake->energy.Tsurf_fbcount = 0;
    lake->energy.unfrozen = 0.0;
    for (i = 0; i < MAX_FRONTS; i++) {
        lake->energy.fdepth[i] = 0.0;
        lake->energy.tdepth[i] = 0.0;
    }
    for (i = 0; i < 2; i++) {
        lake->energy.Cs[i] = 0.0;
        lake->energy.kappa[i] = 0.0;
    }
    for (i = 0; i < MAX_NODES; i++) {
        lake->energy.Cs_node[i] = 0.0;
        lake->energy.ice[i] = 0.0;
        lake->energy.kappa_node[i] = 0.0;
        lake->energy.moist[i] = 0.0;
        lake->energy.T[i] = lake->temp[0];
        lake->energy.T_fbflag[i] = 0;
        lake->energy.T_fbcount[i] = 0;
    }
    // Energy fluxes
    lake->energy.advected_sensible = 0.0;
    lake->energy.advection = 0.0;
    lake->energy.AtmosError = 0.0;
    lake->energy.AtmosLatent = 0.0;
    lake->energy.AtmosLatentSub = 0.0;
    lake->energy.AtmosSensible = 0.0;
    lake->energy.canopy_advection = 0.0;
    lake->energy.canopy_latent = 0.0;
    lake->energy.canopy_latent_sub = 0.0;
    lake->energy.canopy_refreeze = 0.0;
    lake->energy.canopy_sensible = 0.0;
    lake->energy.deltaCC = 0.0;
    lake->energy.deltaH = 0.0;
    lake->energy.error = 0.0;
    lake->energy.fusion = 0.0;
    lake->energy.grnd_flux = 0.0;
    lake->energy.latent = 0.0;
    lake->energy.latent_sub = 0.0;
    lake->energy.longwave = 0.0;
    lake->energy.LongOverIn = 0.0;
    lake->energy.LongUnderIn = 0.0;
    lake->energy.LongUnderOut = 0.0;
    lake->energy.melt_energy = 0.0;
    lake->energy.NetLongAtmos = 0.0;
    lake->energy.NetLongOver = 0.0;
    lake->energy.NetLongUnder = 0.0;
    lake->energy.NetShortAtmos = 0.0;
    lake->energy.NetShortGrnd = 0.0;
    lake->energy.NetShortOver = 0.0;
    lake->energy.NetShortUnder = 0.0;
    lake->energy.out_long_canopy = 0.0;
    lake->energy.out_long_surface = 0.0;
    lake->energy.refreeze_energy = 0.0;
    lake->energy.sensible = 0.0;
    lake->energy.shortwave = 0.0;
    lake->energy.ShortOverIn = 0.0;
    lake->energy.ShortUnderIn = 0.0;
    lake->energy.snow_flux = 0.0;
    // Soil states and fluxes
    lake->soil.asat = 1.0;
    if (!skip_hydro) {
        lake->soil.baseflow = 0.0;
        lake->soil.inflow = 0.0;
        lake->soil.runoff = 0.0;
    }
    lake->soil.rootmoist = 0.0;
    lake->soil.wetness = 1.0;
    for (i = 0; i < 2; i++) {
        lake->soil.aero_resist[i] = 0.0;
    }
    for (i = 0; i < MAX_LAYERS; i++) {
        lake->soil.layer[i].Cs = cell->layer[i].Cs;
        lake->soil.layer[i].T = lake->temp[0];
        lake->soil.layer[i].evap = 0.0;
        lake->soil.layer[i].kappa = cell->layer[i].kappa;
        lake->soil.layer[i].moist = soil_con->porosity[i] * soil_con->depth[i] *
                                    MM_PER_M;
        lake->soil.layer[i].phi = cell->layer[i].phi;
        for (j = 0; j < options.Nfrost; j++) {
            lake->soil.layer[i].ice[j] = 0.0;
        }
    }
    lake->soil.zwt = 0.0;
    lake->soil.zwt_lumped = 0.0;
    if (!skip_hydro) {
        lake->soil.pot_evap = 0.0;
    }
    if (options.CARBON) {
        lake->soil.RhLitter = 0.0;
        lake->soil.RhLitter2Atm = 0.0;
        lake->soil.RhInter = 0.0;
        lake->soil.RhSlow = 0.0;
        lake->soil.RhTot = 0.0;
        lake->soil.CLitter = 0.0;
        lake->soil.CInter = 0.0;
        lake->soil.CSlow = 0.0;
    }

    return(0);
}

