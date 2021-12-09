/******************************************************************************
 * @section DESCRIPTION
 *
 * This subroutine controls the order and number of forcing variables read from
 * the forcing data files.  Two forcing files are allowed, variables, time step
 * and file format must be defined in the global control file.
 *****************************************************************************/

#include <vic_driver_classic.h>

/******************************************************************************
 * @brief    Control the order and number of forcing variables read from the
 *           forcing data files.
 *****************************************************************************/
double **
read_forcing_data(FILE              **infile,
                  global_param_struct global_param,
                  double          ****veg_hist_data)
{
    extern param_set_struct param_set;
    extern size_t           NF;

    size_t                  i, j;
    double                **forcing_data;

    /** Allocate data arrays for input forcing data **/
    forcing_data = calloc(N_FORCING_TYPES, sizeof(*forcing_data));
    check_alloc_status(forcing_data, "Memory allocation error.");
    (*veg_hist_data) = calloc(N_FORCING_TYPES, sizeof(*(*veg_hist_data)));
    check_alloc_status((*veg_hist_data), "Memory allocation error.");
    for (i = 0; i < N_FORCING_TYPES; i++) {
        if (param_set.TYPE[i].SUPPLIED) {
            if (i != ALBEDO && i != LAI && i != FCANOPY) {
                forcing_data[i] = calloc(global_param.nrecs * NF,
                                         sizeof(*(forcing_data[i])));
                check_alloc_status(forcing_data[i], "Memory allocation error.");
            }
            else {
                (*veg_hist_data)[i] = calloc(param_set.TYPE[i].N_ELEM,
                                             sizeof(*((*veg_hist_data)[i])));
                check_alloc_status((*veg_hist_data)[i],
                                   "Memory allocation error.");
                for (j = 0; j < param_set.TYPE[i].N_ELEM; j++) {
                    (*veg_hist_data)[i][j] = calloc(global_param.nrecs * NF,
                                                    sizeof(*((*veg_hist_data)[i]
                                                             [j])));
                    check_alloc_status((*veg_hist_data)[i][j],
                                       "Memory allocation error.");
                }
            }
        }
    }

    /** Read First Forcing Data File **/
    if (param_set.FORCE_DT[0] > 0) {
        read_atmos_data(infile[0], global_param, 0, global_param.forceskip[0],
                        forcing_data, (*veg_hist_data));
    }
    else {
        log_err("File time step must be defined for at least the first "
                "forcing file (FILE_DT).");
    }

    /** Read Second Forcing Data File **/
    if (param_set.FORCE_DT[1] > 0) {
        read_atmos_data(infile[1], global_param, 1, global_param.forceskip[1],
                        forcing_data, (*veg_hist_data));
    }

    return(forcing_data);
}
