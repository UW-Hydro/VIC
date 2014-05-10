#include <vic_def.h>
#include <vic_run.h>
#include <vic_driver_image.h>

void
vic_start()
{
    extern filenames_struct    filenames;
    extern filep_struct        filep;
    extern global_param_struct global_param;
    extern domain_struct       global_domain;

    // read global settings
    initialize_global();
    filep.globalparam = open_file(filenames.global, "r");
    global_param = get_global_param(&filenames, filep.globalparam);

    // read domain info
    get_global_domain(filenames.domain, &global_domain);

    // decompose the mask
}
