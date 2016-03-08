/******************************************************************************
 * @section DESCRIPTION
 *
 * Print library.
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

#include <vic_driver_shared_image.h>

/******************************************************************************
 * @brief    Print atmos data structure.
 *****************************************************************************/
void
print_atmos_data(atmos_data_struct *atmos)
{
    extern option_struct options;

    fprintf(LOG_DEST, "atmos_data  :\n");
    fprintf(LOG_DEST, "\tair_temp  : %.4f\n", atmos->air_temp[0]);
    fprintf(LOG_DEST, "\tdensity   : %.4f\n", atmos->density[0]);
    fprintf(LOG_DEST, "\tlongwave  : %.4f\n", atmos->longwave[0]);
    fprintf(LOG_DEST, "\tout_prec  : %.4f\n", atmos->out_prec);
    fprintf(LOG_DEST, "\tout_rain  : %.4f\n", atmos->out_rain);
    fprintf(LOG_DEST, "\tout_snow  : %.4f\n", atmos->out_snow);
    fprintf(LOG_DEST, "\tprec      : %.4f\n", atmos->prec[0]);
    fprintf(LOG_DEST, "\tpressure  : %.4f\n", atmos->pressure[0]);
    fprintf(LOG_DEST, "\tshortwave : %.4f\n", atmos->shortwave[0]);
    fprintf(LOG_DEST, "\tsnowflag  : %d\n", atmos->snowflag[0]);
    fprintf(LOG_DEST, "\tvp        : %.4f\n", atmos->vp[0]);
    fprintf(LOG_DEST, "\tvpd       : %.4f\n", atmos->vpd[0]);
    fprintf(LOG_DEST, "\twind      : %.4f\n", atmos->wind[0]);
    if (options.LAKES) {
        fprintf(LOG_DEST, "\tchannel_in: %.4f\n", atmos->channel_in[0]);
    }
    if (options.CARBON) {
        fprintf(LOG_DEST, "\tCatm      : %.4f\n", atmos->Catm[0]);
        fprintf(LOG_DEST, "\tfdir      : %.4f\n", atmos->fdir[0]);
        fprintf(LOG_DEST, "\tpar       : %.4f\n", atmos->par[0]);
    }
}

/******************************************************************************
 * @brief    Print energy balance structure.
 *****************************************************************************/
void
print_domain(domain_struct *domain,
             bool           print_loc)
{
    size_t i;

    fprintf(LOG_DEST, "domain:\n");
    fprintf(LOG_DEST, "\tncells_total : %zd\n", domain->ncells_total);
    fprintf(LOG_DEST, "\tncells_active: %zd\n", domain->ncells_active);
    fprintf(LOG_DEST, "\tn_nx         : %zd\n", domain->n_nx);
    fprintf(LOG_DEST, "\tn_ny         : %zd\n", domain->n_ny);
    fprintf(LOG_DEST, "\tlocations    : %p\n", domain->locations);
    if (print_loc) {
        for (i = 0; i < domain->ncells_active; i++) {
            print_location(&(domain->locations[i]));
        }
    }
}

/******************************************************************************
 * @brief    Print location structure.
 *****************************************************************************/
void
print_location(location_struct *loc)
{
    fprintf(LOG_DEST, "location:\n");
    fprintf(LOG_DEST, "\tlatitude       : %.4f\n", loc->latitude);
    fprintf(LOG_DEST, "\tlongitude      : %.4f\n", loc->longitude);
    fprintf(LOG_DEST, "\tarea           : %.4f\n", loc->area);
    fprintf(LOG_DEST, "\tfrac           : %.4f\n", loc->frac);
    fprintf(LOG_DEST, "\tnveg           : %zd\n", loc->nveg);
    fprintf(LOG_DEST, "\tglobal_idx     : %zd\n", loc->global_idx);
    fprintf(LOG_DEST, "\tio_idx         : %zd\n", loc->io_idx);
    fprintf(LOG_DEST, "\tlocal_idx      : %zd\n", loc->local_idx);
}

/******************************************************************************
 * @brief    Print location structure as one string.
 *****************************************************************************/
void
sprint_location(char            *str,
                location_struct *loc)
{
    sprintf(str,
            "location:\n"
            "\tlatitude       : %.4f\n"
            "\tlongitude      : %.4f\n"
            "\tarea           : %.4f\n"
            "\tfrac           : %.4f\n"
            "\nveg            : %zd\n"
            "\tglobal_idx     : %zd\n"
            "\tio_idx         : %zd\n"
            "\tlocal_idx      : %zd\n",
            loc->latitude, loc->longitude, loc->area, loc->frac,
            loc->nveg, loc->global_idx, loc->io_idx, loc->local_idx);
}

/******************************************************************************
 * @brief    Print netCDF file structure.
 *****************************************************************************/
void
print_nc_file(nc_file_struct *nc)
{
    fprintf(LOG_DEST, "nc_file:");
    fprintf(LOG_DEST, "\tfname          : %s\n", nc->fname);
    fprintf(LOG_DEST, "\tc_fillvalue    : %d\n", nc->c_fillvalue);
    fprintf(LOG_DEST, "\ti_fillvalue    : %d\n", nc->i_fillvalue);
    fprintf(LOG_DEST, "\td_fillvalue    : %.4f\n", nc->d_fillvalue);
    fprintf(LOG_DEST, "\tf_fillvalue    : %.4f\n", nc->f_fillvalue);
    fprintf(LOG_DEST, "\tnc_id          : %d\n", nc->nc_id);
    fprintf(LOG_DEST, "\tband_dimid     : %d\n", nc->band_dimid);
    fprintf(LOG_DEST, "\tfront_dimid    : %d\n", nc->front_dimid);
    fprintf(LOG_DEST, "\tfrost_dimid    : %d\n", nc->frost_dimid);
    fprintf(LOG_DEST, "\tlayer_dimid    : %d\n", nc->layer_dimid);
    fprintf(LOG_DEST, "\tni_dimid       : %d\n", nc->ni_dimid);
    fprintf(LOG_DEST, "\tnj_dimid       : %d\n", nc->nj_dimid);
    fprintf(LOG_DEST, "\tnode_dimid     : %d\n", nc->node_dimid);
    fprintf(LOG_DEST, "\troot_zone_dimid: %d\n", nc->root_zone_dimid);
    fprintf(LOG_DEST, "\ttime_dimid     : %d\n", nc->time_dimid);
    fprintf(LOG_DEST, "\tveg_dimid      : %d\n", nc->veg_dimid);
    fprintf(LOG_DEST, "\tband_size      : %zd\n", nc->band_size);
    fprintf(LOG_DEST, "\tfront_size     : %zd\n", nc->front_size);
    fprintf(LOG_DEST, "\tfrost_size     : %zd\n", nc->frost_size);
    fprintf(LOG_DEST, "\tlayer_size     : %zd\n", nc->layer_size);
    fprintf(LOG_DEST, "\tni_size        : %zd\n", nc->ni_size);
    fprintf(LOG_DEST, "\tnj_size        : %zd\n", nc->nj_size);
    fprintf(LOG_DEST, "\tnode_size      : %zd\n", nc->node_size);
    fprintf(LOG_DEST, "\troot_zone_size : %zd\n", nc->root_zone_size);
    fprintf(LOG_DEST, "\ttime_size      : %zd\n", nc->time_size);
    fprintf(LOG_DEST, "\tveg_size       : %zd\n", nc->veg_size);
    fprintf(LOG_DEST, "\topen           : %d\n", nc->open);
}

/******************************************************************************
 * @brief    Print netCDF variable structure.
 *****************************************************************************/
void
print_nc_var(nc_var_struct *nc_var,
             size_t         ndims)
{
    size_t i;

    fprintf(LOG_DEST, "nc_var:\n");
    fprintf(LOG_DEST, "\tnc_var_name: %s\n", nc_var->nc_var_name);
    fprintf(LOG_DEST, "\tnc_units: %s\n", nc_var->nc_units);
    fprintf(LOG_DEST, "\tnc_dimids:");
    for (i = 0; i < ndims; i++) {
        fprintf(LOG_DEST, "\t%d", nc_var->nc_dimids[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tnc_counts:");
    for (i = 0; i < ndims; i++) {
        fprintf(LOG_DEST, "\t%d", nc_var->nc_counts[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tnc_type: %d\n", nc_var->nc_type);
    fprintf(LOG_DEST, "\tnc_aggtype: %d\n", nc_var->nc_aggtype);
    fprintf(LOG_DEST, "\tnc_dims: %d\n", nc_var->nc_dims);
    fprintf(LOG_DEST, "\tnc_write: %d\n", nc_var->nc_write);
}

/******************************************************************************
 * @brief    Print veg_con_map structure.
 *****************************************************************************/
void
print_veg_con_map(veg_con_map_struct *veg_con_map)
{
    size_t i;

    fprintf(LOG_DEST, "veg_con_map:\n");
    fprintf(LOG_DEST, "\tnv_types : %zd\n", veg_con_map->nv_types);
    fprintf(LOG_DEST, "\tnv_active: %zd\n", veg_con_map->nv_active);
    for (i = 0; i < veg_con_map->nv_types; i++) {
        fprintf(LOG_DEST, "\t%zd      : %d (vidx) %f (Cv)\n", i,
                veg_con_map->vidx[i],
                veg_con_map->Cv[i]);
    }
}
