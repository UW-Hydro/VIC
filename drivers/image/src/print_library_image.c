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

#include <vic_def.h>
#include <vic_run.h>
#include <vic_driver_image.h>

/******************************************************************************
 * @brief    Print atmos data structure.
 *****************************************************************************/
void
print_atmos_data(atmos_data_struct *atmos)
{
    printf("atmos_data  :\n");
    printf("\tair_temp  : %.4lf\n", atmos->air_temp[0]);
    printf("\tCatm      : %.4lf\n", atmos->Catm[0]);
    printf("\tchannel_in: %.4lf\n", atmos->channel_in[0]);
    printf("\tcoszen    : %.4lf\n", atmos->coszen[0]);
    printf("\tdensity   : %.4lf\n", atmos->density[0]);
    printf("\tfdir      : %.4lf\n", atmos->fdir[0]);
    printf("\tlongwave  : %.4lf\n", atmos->longwave[0]);
    printf("\tout_prec  : %.4lf\n", atmos->out_prec);
    printf("\tout_rain  : %.4lf\n", atmos->out_rain);
    printf("\tout_snow  : %.4lf\n", atmos->out_snow);
    printf("\tpar       : %.4lf\n", atmos->par[0]);
    printf("\tprec      : %.4lf\n", atmos->prec[0]);
    printf("\tpressure  : %.4lf\n", atmos->pressure[0]);
    printf("\tshortwave : %.4lf\n", atmos->shortwave[0]);
    printf("\tsnowflag  : %d\n", atmos->snowflag[0]);
    printf("\ttskc      : %.4lf\n", atmos->tskc[0]);
    printf("\tvp        : %.4lf\n", atmos->vp[0]);
    printf("\tvpd       : %.4lf\n", atmos->vpd[0]);
    printf("\twind      : %.4lf\n", atmos->wind[0]);
}

/******************************************************************************
 * @brief    Print energy balance structure.
 *****************************************************************************/
void
print_domain(domain_struct *domain,
             bool           print_loc)
{
    size_t i;

    printf("domain:\n");
    printf("\tncells_global: %zd\n", domain->ncells_global);
    printf("\tn_nx         : %zd\n", domain->n_nx);
    printf("\tn_ny         : %zd\n", domain->n_ny);
    printf("\tncells_local : %zd\n", domain->ncells_local);
    printf("\tlocations    : %p\n", domain->locations);
    if (print_loc) {
        for (i = 0; i < domain->ncells_global; i++) {
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
    printf("location:\n");
    printf("\tlatitude       : %.4lf\n", loc->latitude);
    printf("\tlongitude      : %.4lf\n", loc->longitude);
    printf("\tarea           : %.4lf\n", loc->area);
    printf("\tfrac           : %.4lf\n", loc->frac);
    printf("\tglobal_cell_idx: %zd\n", loc->global_cell_idx);
    printf("\tglobal_x_idx   : %zd\n", loc->global_x_idx);
    printf("\tglobal_y_idx   : %zd\n", loc->global_y_idx);
    printf("\tlocal_cell_idx : %zd\n", loc->local_cell_idx);
    printf("\tlocal_x_idx    : %zd\n", loc->local_x_idx);
    printf("\tlocal_y_idx    : %zd\n", loc->local_y_idx);
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
            "\tlatitude       : %.4lf\n"
            "\tlongitude      : %.4lf\n"
            "\tarea           : %.4lf\n"
            "\tfrac           : %.4lf\n"
            "\tglobal_cell_idx: %zd\n"
            "\tglobal_x_idx   : %zd\n"
            "\tglobal_y_idx   : %zd\n"
            "\tlocal_cell_idx : %zd\n"
            "\tlocal_x_idx    : %zd\n"
            "\tlocal_y_idx    : %zd\n",
            loc->latitude, loc->longitude, loc->area, loc->frac,
            loc->global_cell_idx, loc->global_x_idx, loc->global_y_idx,
            loc->local_cell_idx, loc->local_x_idx, loc->local_y_idx);
}

/******************************************************************************
 * @brief    Print netCDF file structure.
 *****************************************************************************/
void
print_nc_file(nc_file_struct *nc)
{
    printf("nc_file:");
    printf("\tfname          : %s\n", nc->fname);
    printf("\tc_fillvalue    : %d\n", nc->c_fillvalue);
    printf("\ti_fillvalue    : %d\n", nc->i_fillvalue);
    printf("\td_fillvalue    : %.4lf\n", nc->d_fillvalue);
    printf("\tf_fillvalue    : %.4f\n", nc->f_fillvalue);
    printf("\tnc_id          : %d\n", nc->nc_id);
    printf("\tband_dimid     : %d\n", nc->band_dimid);
    printf("\tfront_dimid    : %d\n", nc->front_dimid);
    printf("\tfrost_dimid    : %d\n", nc->frost_dimid);
    printf("\tlayer_dimid    : %d\n", nc->layer_dimid);
    printf("\tni_dimid       : %d\n", nc->ni_dimid);
    printf("\tnj_dimid       : %d\n", nc->nj_dimid);
    printf("\tnode_dimid     : %d\n", nc->node_dimid);
    printf("\troot_zone_dimid: %d\n", nc->root_zone_dimid);
    printf("\ttime_dimid     : %d\n", nc->time_dimid);
    printf("\tveg_dimid      : %d\n", nc->veg_dimid);
    printf("\tband_size      : %zd\n", nc->band_size);
    printf("\tfront_size     : %zd\n", nc->front_size);
    printf("\tfrost_size     : %zd\n", nc->frost_size);
    printf("\tlayer_size     : %zd\n", nc->layer_size);
    printf("\tni_size        : %zd\n", nc->ni_size);
    printf("\tnj_size        : %zd\n", nc->nj_size);
    printf("\tnode_size      : %zd\n", nc->node_size);
    printf("\troot_zone_size : %zd\n", nc->root_zone_size);
    printf("\ttime_size      : %zd\n", nc->time_size);
    printf("\tveg_size       : %zd\n", nc->veg_size);
    printf("\topen           : %d\n", nc->open);
}

/******************************************************************************
 * @brief    Print netCDF variable structure.
 *****************************************************************************/
void
print_nc_var(nc_var_struct *nc_var,
             size_t         ndims)
{
    size_t i;

    printf("nc_var:\n");
    printf("\tnc_var_name: %s\n", nc_var->nc_var_name);
    printf("\tnc_units: %s\n", nc_var->nc_units);
    printf("\tnc_dimids:");
    for (i = 0; i < ndims; i++) {
        printf("\t%d", nc_var->nc_dimids[i]);
    }
    printf("\n");
    printf("\tnc_counts:");
    for (i = 0; i < ndims; i++) {
        printf("\t%d", nc_var->nc_counts[i]);
    }
    printf("\n");
    printf("\tnc_type: %d\n", nc_var->nc_type);
    printf("\tnc_aggtype: %d\n", nc_var->nc_aggtype);
    printf("\tnc_dims: %d\n", nc_var->nc_dims);
    printf("\tnc_write: %d\n", nc_var->nc_write);
}

/******************************************************************************
 * @brief    Print veg_con_map structure.
 *****************************************************************************/
void
print_veg_con_map(veg_con_map_struct *veg_con_map)
{
    size_t i;

    printf("veg_con_map:\n");
    printf("\tnv_types : %zd\n", veg_con_map->nv_types);
    printf("\tnv_active: %zd\n", veg_con_map->nv_active);
    for (i = 0; i < veg_con_map->nv_types; i++) {
        printf("\t%zd      : %d (vidx) %lf (Cv)\n", i, veg_con_map->vidx[i],
               veg_con_map->Cv[i]);
    }
}
