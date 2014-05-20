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
    extern option_struct       options;

    // read global settings
    initialize_global();
    filep.globalparam = open_file(filenames.global, "r");
    global_param = get_global_param(&filenames, filep.globalparam);

    // read domain info
    get_global_domain(filenames.domain, &global_domain);
    // print_domain(&global_domain, true);

    // decompose the mask

    // get dimensions (number of vegetation types, soil zones, etc)
    options.ROOT_ZONES = get_nc_dimension(filenames.soil, "root_zone");
    options.Nlayer = get_nc_dimension(filenames.soil, "nlayer");
    options.NVEGTYPES = get_nc_dimension(filenames.veg, "veg_class");
    options.SNOW_BAND = get_nc_dimension(filenames.snowband, "snow_band");
}
