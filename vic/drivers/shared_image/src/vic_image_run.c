/******************************************************************************
 * @section DESCRIPTION
 *
 * Run function for image mode driver.
 *****************************************************************************/

#include <vic_driver_shared_image.h>
#include <rout.h>

/******************************************************************************
 * @brief    Run VIC for one timestep and store output data
 *****************************************************************************/
void
vic_image_run(dmy_struct *dmy_current)
{
    extern size_t              current;
    extern all_vars_struct    *all_vars;
    extern force_data_struct  *force;
    extern domain_struct       local_domain;
    extern option_struct       options;
    extern global_param_struct global_param;
    extern lake_con_struct    *lake_con;
    extern double           ***out_data;
    extern stream_struct      *output_streams;
    extern save_data_struct   *save_data;
    extern soil_con_struct    *soil_con;
    extern veg_con_struct    **veg_con;
    extern veg_hist_struct   **veg_hist;
    extern veg_lib_struct    **veg_lib;

    char                       dmy_str[MAXSTRING];
    size_t                     i;
    timer_struct               timer;

    // Print the current timestep info before running vic_run
    sprint_dmy(dmy_str, dmy_current);
    debug("Running timestep %zu: %s", current, dmy_str);

    // If running with OpenMP, run this for loop using multiple threads
    #pragma omp parallel for default(shared) private(i, timer, vic_run_ref_str)
    for (i = 0; i < local_domain.ncells_active; i++) {
        // Set global reference string (for debugging inside vic_run)
        sprintf(vic_run_ref_str, "Gridcell io_idx: %zu, timestep info: %s",
                local_domain.locations[i].io_idx, dmy_str);

        update_step_vars(&(all_vars[i]), veg_con[i], veg_hist[i]);

        timer_start(&timer);
        vic_run(&(force[i]), &(all_vars[i]), dmy_current, &global_param,
                &(lake_con[i]), &(soil_con[i]), veg_con[i], veg_lib[i]);
        timer_stop(&timer);

        put_data(&(all_vars[i]), &(force[i]), &(soil_con[i]), veg_con[i],
                 veg_lib[i], &(lake_con[i]), out_data[i], &(save_data[i]),
                 &timer);
    }

    // run routing over the domain
    rout_run();     // Routing routine (extension)

    for (i = 0; i < options.Noutstreams; i++) {
        agg_stream_data(&(output_streams[i]), dmy_current, out_data);
    }
}
