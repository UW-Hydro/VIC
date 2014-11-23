#include <vic_def.h>
#include <vic_run.h>
#include <vic_driver_image.h>

char               *version = "VIC 5.0 image";

int                 NF, NR;
size_t              current;
all_vars_struct    *all_vars = NULL;
atmos_data_struct  *atmos = NULL;
dmy_struct         *dmy = NULL;
filenames_struct    filenames;
filep_struct        filep;
domain_struct       global_domain;
global_param_struct global_param;
lake_con_struct     lake_con;
nc_file_struct      nc_hist_file;
nc_var_struct       nc_vars[N_OUTVAR_TYPES];
option_struct       options;
out_data_struct   **out_data;
save_data_struct   *save_data;
param_set_struct    param_set;
soil_con_struct    *soil_con = NULL;
veg_con_map_struct *veg_con_map = NULL;
veg_con_struct    **veg_con = NULL;
veg_hist_struct   **veg_hist = NULL;
veg_lib_struct    **veg_lib = NULL;

int
main(int    argc,
     char **argv)
{
    // process command line arguments
    cmd_proc(argc, argv, filenames.global);

    // read global parameters
    vic_start();

    // allocate memory
    vic_alloc();

    // initialize model parameters from parameter files
    vic_init();

    // restore model state, either using a cold start or from a restart file
    vic_restore();

    // initialize output structures
    vic_init_output();

    // loop over all timesteps
    for (current = 0; current < global_param.nrecs; current++) {
        // read forcing data
        vic_force();

        // run vic over the domain
        vic_image_run();

        // if output:
        vic_write();

        // if save:
        if (current == global_param.nrecs - 1) {
            vic_store();
        }
    }

    // clean up
    vic_finalize();

    return EXIT_SUCCESS;
}
