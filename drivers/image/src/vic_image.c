#include <vic_def.h>
#include <vic_run.h>
#include <vic_driver_image.h>

char               *version = "VIC 5.0 image";

int                 NF, NR;
filenames_struct    filenames;
filep_struct        filep;
global_param_struct global_param;
option_struct       options;
param_set_struct    param_set;
domain_struct       global_domain;

int
main(int    argc,
     char **argv)
{
    cmd_proc(argc, argv, filenames.global);

    vic_start();
    vic_alloc();
    // vic_init()
    // vic_restore()
    // foreach timestep:
    // foreach gridcell:
    // vic_run()
    // if output:
    // vic_write()
    // if save:
    // vic_save()
    // vic_finalize()

    return EXIT_SUCCESS;
}
