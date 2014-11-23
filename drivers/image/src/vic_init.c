#include <vic_def.h>
#include <vic_run.h>
#include <vic_driver_image.h>

void
vic_init(void)
{
    extern all_vars_struct    *all_vars;
    extern size_t              current;
    extern dmy_struct         *dmy;
    extern domain_struct       global_domain;
    extern option_struct       options;
    extern filenames_struct    filenames;
    extern global_param_struct global_param;
    extern soil_con_struct    *soil_con;
    extern veg_con_map_struct *veg_con_map;
    extern veg_con_struct    **veg_con;
    // extern veg_hist_struct   **veg_hist;
    extern veg_lib_struct    **veg_lib;

    bool                       found;
    bool                       no_overstory;
    char                       errstr[MAXSTRING];
    char                       locstr[MAXSTRING];
    double                     mean;
    double                     sum;
    double                    *dvar = NULL;
    size_t                     i;
    size_t                     j;
    size_t                     k;
    size_t                    *idx;
    size_t                     nveg;
    int                        veg_class;
    int                        vidx;
    size_t                     d2count[2];
    size_t                     d2start[2];
    size_t                     d3count[3];
    size_t                     d3start[3];
    size_t                     d4count[4];
    size_t                     d4start[4];

    // make_dmy()
    dmy = make_dmy(&global_param);

    // allocate memory for variables to be read
    dvar = (double *) malloc(global_domain.n_ny * global_domain.n_nx *
                             sizeof(double));
    if (dvar == NULL) {
        nrerror("Memory allocation error in vic_init().");
    }

    // The method used to convert the NetCDF fields to VIC structures for
    // individual grid cells is to read a 2D slice and then loop over the
    // domain cells to assign the values to the VIC structures

    d2start[0] = 0;
    d2start[1] = 0;
    d2count[0] = global_domain.n_ny;
    d2count[1] = global_domain.n_nx;

    d3start[0] = 0;
    d3start[1] = 0;
    d3start[2] = 0;
    d3count[0] = 1;
    d3count[1] = global_domain.n_ny;
    d3count[2] = global_domain.n_nx;

    d4start[0] = 0;
    d4start[1] = 0;
    d4start[2] = 0;
    d4start[3] = 0;
    d4count[0] = 1;
    d4count[1] = 1;
    d4count[2] = global_domain.n_ny;
    d4count[3] = global_domain.n_nx;

    // start the clock
    current = 0;

    // get 1D indices used in mapping the netcdf fields to the locations
    idx = (size_t *) malloc(global_domain.ncells_global *
                            sizeof(size_t));
    if (idx == NULL) {
        nrerror("Memory allocation error in vic_init().");
    }
    for (i = 0; i < global_domain.ncells_global; i++) {
        idx[i] = get_global_idx(&global_domain, i);
    }

    // read_veglib()

    // TBD: Check that options.NVEGTYPES is the right number. Bare soil is
    // the complicating factor here
    for (i = 0; i < global_domain.ncells_global; i++) {
        for (j = 0; j < options.NVEGTYPES; j++) {
            veg_lib[i][j].NVegLibTypes = options.NVEGTYPES;
            veg_lib[i][j].veg_class = (int) j;
        }
    }

    // overstory
    for (j = 0; j < options.NVEGTYPES; j++) {
        d3start[0] = j;
        get_nc_field_double(filenames.veglib, "overstory",
                            d3start, d3count, dvar);
        for (i = 0; i < global_domain.ncells_global; i++) {
            veg_lib[i][j].overstory = (int) dvar[idx[i]];
        }
    }

    // rarc
    for (j = 0; j < options.NVEGTYPES; j++) {
        d3start[0] = j;
        get_nc_field_double(filenames.veglib, "rarc",
                            d3start, d3count, dvar);
        for (i = 0; i < global_domain.ncells_global; i++) {
            veg_lib[i][j].rarc = (double) dvar[idx[i]];
        }
    }

    // rmin
    for (j = 0; j < options.NVEGTYPES; j++) {
        d3start[0] = j;
        get_nc_field_double(filenames.veglib, "rmin",
                            d3start, d3count, dvar);
        for (i = 0; i < global_domain.ncells_global; i++) {
            veg_lib[i][j].rmin = (double) dvar[idx[i]];
        }
    }

    // wind height
    for (j = 0; j < options.NVEGTYPES; j++) {
        d3start[0] = j;
        get_nc_field_double(filenames.veglib, "wind_h",
                            d3start, d3count, dvar);
        for (i = 0; i < global_domain.ncells_global; i++) {
            veg_lib[i][j].wind_h = (double) dvar[idx[i]];
        }
    }

    // RGL
    for (j = 0; j < options.NVEGTYPES; j++) {
        d3start[0] = j;
        get_nc_field_double(filenames.veglib, "RGL",
                            d3start, d3count, dvar);
        for (i = 0; i < global_domain.ncells_global; i++) {
            veg_lib[i][j].RGL = (float) dvar[idx[i]];
        }
    }

    // rad_atten
    for (j = 0; j < options.NVEGTYPES; j++) {
        d3start[0] = j;
        get_nc_field_double(filenames.veglib, "rad_atten",
                            d3start, d3count, dvar);
        for (i = 0; i < global_domain.ncells_global; i++) {
            veg_lib[i][j].rad_atten = (double) dvar[idx[i]];
        }
    }

    // wind_atten
    for (j = 0; j < options.NVEGTYPES; j++) {
        d3start[0] = j;
        get_nc_field_double(filenames.veglib, "wind_atten",
                            d3start, d3count, dvar);
        for (i = 0; i < global_domain.ncells_global; i++) {
            veg_lib[i][j].wind_atten = (double) dvar[idx[i]];
        }
    }

    // trunk_ratio
    for (j = 0; j < options.NVEGTYPES; j++) {
        d3start[0] = j;
        get_nc_field_double(filenames.veglib, "trunk_ratio",
                            d3start, d3count, dvar);
        for (i = 0; i < global_domain.ncells_global; i++) {
            veg_lib[i][j].trunk_ratio = (double) dvar[idx[i]];
        }
    }

    // LAI and Wdmax
    for (j = 0; j < options.NVEGTYPES; j++) {
        d4start[0] = j;
        for (k = 0; k < MONTHSPERYEAR; k++) {
            d4start[1] = k;
            get_nc_field_double(filenames.veglib, "LAI",
                                d4start, d4count, dvar);
            for (i = 0; i < global_domain.ncells_global; i++) {
                veg_lib[i][j].LAI[k] = (double) dvar[idx[i]];
                veg_lib[i][j].Wdmax[k] = LAI_WATER_FACTOR *
                                         veg_lib[i][j].LAI[k];
            }
        }
    }

    // albedo
    for (j = 0; j < options.NVEGTYPES; j++) {
        d4start[0] = j;
        for (k = 0; k < MONTHSPERYEAR; k++) {
            d4start[1] = k;
            get_nc_field_double(filenames.veglib, "albedo",
                                d4start, d4count, dvar);
            for (i = 0; i < global_domain.ncells_global; i++) {
                veg_lib[i][j].albedo[k] = (double) dvar[idx[i]];
            }
        }
    }

    // veg_rough
    for (j = 0; j < options.NVEGTYPES; j++) {
        d4start[0] = j;
        for (k = 0; k < MONTHSPERYEAR; k++) {
            d4start[1] = k;
            get_nc_field_double(filenames.veglib, "veg_rough",
                                d4start, d4count, dvar);
            for (i = 0; i < global_domain.ncells_global; i++) {
                veg_lib[i][j].roughness[k] = (double) dvar[idx[i]];
            }
        }
    }

    // displacement
    for (j = 0; j < options.NVEGTYPES; j++) {
        d4start[0] = j;
        for (k = 0; k < MONTHSPERYEAR; k++) {
            d4start[1] = k;
            get_nc_field_double(filenames.veglib, "displacement",
                                d4start, d4count, dvar);
            for (i = 0; i < global_domain.ncells_global; i++) {
                veg_lib[i][j].displacement[k] = (double) dvar[idx[i]];
            }
        }
    }

    // vegcover not implemented in image model
    for (j = 0; j < options.NVEGTYPES; j++) {
        for (k = 0; k < MONTHSPERYEAR; k++) {
            for (i = 0; i < global_domain.ncells_global; i++) {
                veg_lib[i][j].vegcover[k] = 1.0;
            }
        }
    }

    // read_soilparam()

    // b_infilt
    get_nc_field_double(filenames.soil, "infilt",
                        d2start, d2count, dvar);
    for (i = 0; i < global_domain.ncells_global; i++) {
        soil_con[i].b_infilt = (double) dvar[idx[i]];
    }

    // Ds
    get_nc_field_double(filenames.soil, "Ds",
                        d2start, d2count, dvar);
    for (i = 0; i < global_domain.ncells_global; i++) {
        soil_con[i].Ds = (double) dvar[idx[i]];
    }

    // Dsmax
    get_nc_field_double(filenames.soil, "Dsmax",
                        d2start, d2count, dvar);
    for (i = 0; i < global_domain.ncells_global; i++) {
        soil_con[i].Dsmax = (double) dvar[idx[i]];
    }

    // Ws
    get_nc_field_double(filenames.soil, "Ws",
                        d2start, d2count, dvar);
    for (i = 0; i < global_domain.ncells_global; i++) {
        soil_con[i].Ws = (double) dvar[idx[i]];
    }

    // c
    get_nc_field_double(filenames.soil, "c",
                        d2start, d2count, dvar);
    for (i = 0; i < global_domain.ncells_global; i++) {
        soil_con[i].c = (double) dvar[idx[i]];
    }

    // expt: unsaturated hydraulic conductivity exponent for each layer
    for (j = 0; j < options.Nlayer; j++) {
        d3start[0] = j;
        get_nc_field_double(filenames.soil, "expt",
                            d3start, d3count, dvar);
        for (i = 0; i < global_domain.ncells_global; i++) {
            soil_con[i].expt[j] = (double) dvar[idx[i]];
        }
    }

    // Ksat: saturated hydraulic conductivity for each layer
    for (j = 0; j < options.Nlayer; j++) {
        d3start[0] = j;
        get_nc_field_double(filenames.soil, "Ksat",
                            d3start, d3count, dvar);
        for (i = 0; i < global_domain.ncells_global; i++) {
            soil_con[i].Ksat[j] = (double) dvar[idx[i]];
        }
    }

    // init_moist: initial soil moisture for cold start
    for (j = 0; j < options.Nlayer; j++) {
        d3start[0] = j;
        get_nc_field_double(filenames.soil, "init_moist",
                            d3start, d3count, dvar);
        for (i = 0; i < global_domain.ncells_global; i++) {
            soil_con[i].init_moist[j] = (double) dvar[idx[i]];
        }
    }

    // phi_s
    for (j = 0; j < options.Nlayer; j++) {
        d3start[0] = j;
        get_nc_field_double(filenames.soil, "phi_s",
                            d3start, d3count, dvar);
        for (i = 0; i < global_domain.ncells_global; i++) {
            soil_con[i].phi_s[j] = (double) dvar[idx[i]];
        }
    }

    // elevation: mean grid cell elevation
    get_nc_field_double(filenames.soil, "elev",
                        d2start, d2count, dvar);
    for (i = 0; i < global_domain.ncells_global; i++) {
        soil_con[i].elevation = (double) dvar[idx[i]];
    }

    // depth: thickness for each soil layer
    for (j = 0; j < options.Nlayer; j++) {
        d3start[0] = j;
        get_nc_field_double(filenames.soil, "depth",
                            d3start, d3count, dvar);
        for (i = 0; i < global_domain.ncells_global; i++) {
            soil_con[i].depth[j] = (double) dvar[idx[i]];
        }
    }

    // avg_temp: mean grid temperature
    get_nc_field_double(filenames.soil, "avg_T",
                        d2start, d2count, dvar);
    for (i = 0; i < global_domain.ncells_global; i++) {
        soil_con[i].avg_temp = (double) dvar[idx[i]];
    }

    // dp: damping depth
    get_nc_field_double(filenames.soil, "dp",
                        d2start, d2count, dvar);
    for (i = 0; i < global_domain.ncells_global; i++) {
        soil_con[i].dp = (double) dvar[idx[i]];
    }

    // bubble: bubbling pressure for each soil layer
    for (j = 0; j < options.Nlayer; j++) {
        d3start[0] = j;
        get_nc_field_double(filenames.soil, "bubble",
                            d3start, d3count, dvar);
        for (i = 0; i < global_domain.ncells_global; i++) {
            soil_con[i].bubble[j] = (double) dvar[idx[i]];
        }
    }

    // quartz: quartz content for each soil layer
    for (j = 0; j < options.Nlayer; j++) {
        d3start[0] = j;
        get_nc_field_double(filenames.soil, "quartz",
                            d3start, d3count, dvar);
        for (i = 0; i < global_domain.ncells_global; i++) {
            soil_con[i].quartz[j] = (double) dvar[idx[i]];
        }
    }

    // bulk_dens_min: mineral bulk density for each soil layer
    for (j = 0; j < options.Nlayer; j++) {
        d3start[0] = j;
        get_nc_field_double(filenames.soil, "bulk_density",
                            d3start, d3count, dvar);
        for (i = 0; i < global_domain.ncells_global; i++) {
            soil_con[i].bulk_dens_min[j] = (double) dvar[idx[i]];
        }
    }

    // soil_dens_min: mineral soil density for each soil layer
    for (j = 0; j < options.Nlayer; j++) {
        d3start[0] = j;
        get_nc_field_double(filenames.soil, "soil_density",
                            d3start, d3count, dvar);
        for (i = 0; i < global_domain.ncells_global; i++) {
            soil_con[i].soil_dens_min[j] = (double) dvar[idx[i]];
        }
    }


    // organic soils
    if (options.ORGANIC_FRACT) {
        // organic
        for (j = 0; j < options.Nlayer; j++) {
            d3start[0] = j;
            get_nc_field_double(filenames.soil, "organic",
                                d3start, d3count, dvar);
            for (i = 0; i < global_domain.ncells_global; i++) {
                soil_con[i].organic[j] = (double) dvar[idx[i]];
            }
        }

        // bulk_dens_org: organic bulk density for each soil layer
        for (j = 0; j < options.Nlayer; j++) {
            d3start[0] = j;
            get_nc_field_double(filenames.soil, "bulk_density_org",
                                d3start, d3count, dvar);
            for (i = 0; i < global_domain.ncells_global; i++) {
                soil_con[i].bulk_dens_org[j] = (double) dvar[idx[i]];
            }
        }

        // soil_dens_org: organic soil density for each soil layer
        for (j = 0; j < options.Nlayer; j++) {
            d3start[0] = j;
            get_nc_field_double(filenames.soil, "soil_density_org",
                                d3start, d3count, dvar);
            for (i = 0; i < global_domain.ncells_global; i++) {
                soil_con[i].soil_dens_org[j] = (double) dvar[idx[i]];
            }
        }
    }

    // Wcr: critical point for each layer
    // Note this value is  multiplied with the maximum moisture in each layer
    for (j = 0; j < options.Nlayer; j++) {
        d3start[0] = j;
        get_nc_field_double(filenames.soil, "Wcr_FRACT",
                            d3start, d3count, dvar);
        for (i = 0; i < global_domain.ncells_global; i++) {
            soil_con[i].Wcr[j] = (double) dvar[idx[i]];
        }
    }

    // Wpwp: wilting point for each layer
    // Note this value is  multiplied with the maximum moisture in each layer
    for (j = 0; j < options.Nlayer; j++) {
        d3start[0] = j;
        get_nc_field_double(filenames.soil, "Wpwp_FRACT",
                            d3start, d3count, dvar);
        for (i = 0; i < global_domain.ncells_global; i++) {
            soil_con[i].Wpwp[j] = (double) dvar[idx[i]];
        }
    }

    // rough: soil roughness
    get_nc_field_double(filenames.soil, "rough",
                        d2start, d2count, dvar);
    for (i = 0; i < global_domain.ncells_global; i++) {
        soil_con[i].rough = (double) dvar[idx[i]];
    }

    // snow_rough: snow roughness
    get_nc_field_double(filenames.soil, "snow_rough",
                        d2start, d2count, dvar);
    for (i = 0; i < global_domain.ncells_global; i++) {
        soil_con[i].snow_rough = (double) dvar[idx[i]];
    }

    // annual_prec: annual precipitation
    get_nc_field_double(filenames.soil, "annual_prec",
                        d2start, d2count, dvar);
    for (i = 0; i < global_domain.ncells_global; i++) {
        soil_con[i].annual_prec = (double) dvar[idx[i]];
    }

    // resid_moist: residual moisture content for each layer
    for (j = 0; j < options.Nlayer; j++) {
        d3start[0] = j;
        get_nc_field_double(filenames.soil, "resid_moist",
                            d3start, d3count, dvar);
        for (i = 0; i < global_domain.ncells_global; i++) {
            soil_con[i].resid_moist[j] = (double) dvar[idx[i]];
        }
    }

    // fs_active: frozen soil active flag
    get_nc_field_double(filenames.soil, "fs_active",
                        d2start, d2count, dvar);
    for (i = 0; i < global_domain.ncells_global; i++) {
        soil_con[i].FS_ACTIVE = (char) dvar[idx[i]];
    }

    // spatial snow
    if (options.SPATIAL_SNOW) {
        // max_snow_distrib_slope
        get_nc_field_double(filenames.soil, "max_snow_distrib_slope",
                            d2start, d2count, dvar);
        for (i = 0; i < global_domain.ncells_global; i++) {
            soil_con[i].max_snow_distrib_slope = (double) dvar[idx[i]];
        }

        // frost_slope: slope of frozen soil distribution
        get_nc_field_double(filenames.soil, "frost_slope",
                            d2start, d2count, dvar);
        for (i = 0; i < global_domain.ncells_global; i++) {
            soil_con[i].frost_slope = (double) dvar[idx[i]];
        }
    }

    if (options.COMPUTE_TREELINE) {
        // avgJulyAirTemp: average July air temperature
        get_nc_field_double(filenames.soil, "avgJulyAirTemp",
                            d2start, d2count, dvar);
        for (i = 0; i < global_domain.ncells_global; i++) {
            soil_con[i].avgJulyAirTemp = (double) dvar[idx[i]];
        }
    }

    // Additional processing of the soil variables
    for (i = 0; i < global_domain.ncells_global; i++) {
        for (j = 0; j < options.Nlayer; j++) {
            // compute layer properties
            soil_con[i].bulk_density[j] =
                (1 - soil_con[i].organic[j]) * soil_con[i].bulk_dens_min[j] +
                soil_con[i].organic[j] * soil_con[i].bulk_dens_org[j];
            soil_con[i].soil_density[j] =
                (1 - soil_con[i].organic[j]) * soil_con[i].soil_dens_min[j] +
                soil_con[i].organic[j] * soil_con[i].soil_dens_org[j];
            if (soil_con[i].resid_moist[j] == MISSING) {
                soil_con[i].resid_moist[j] = RESID_MOIST;
            }
            soil_con[i].porosity[j] = 1 - soil_con[i].bulk_density[j] /
                                      soil_con[i].soil_density[j];
            soil_con[i].max_moist[j] = soil_con[i].depth[j] *
                                       soil_con[i].porosity[j] * 1000.;
            // check layer thicknesses
            if (soil_con[i].depth[j] < MINSOILDEPTH) {
                sprint_location(locstr, &(global_domain.locations[i]));
                sprintf(errstr, "ERROR: Model will not function with layer "
                        "%zd depth %f < %f m.\n%s", j, soil_con[i].depth[j],
                        MINSOILDEPTH, locstr);
                nrerror(errstr);
            }
        }
        // check relative thickness of top two layers
        if (soil_con[i].depth[0] > soil_con[i].depth[1]) {
            sprint_location(locstr, &(global_domain.locations[i]));
            sprintf(errstr, "ERROR: Model will not function with layer 0 "
                    "thicker than layer 1 (%f m > %f m).\n%s",
                    soil_con[i].depth[0], soil_con[i].depth[1], locstr);
            nrerror(errstr);
        }
        // compute maximum infiltration for upper layers
        if (options.Nlayer == 2) {
            soil_con[i].max_infil = (1.0 + soil_con[i].b_infilt) *
                                    soil_con[i].max_moist[0];
        }
        else {
            soil_con[i].max_infil = (1.0 + soil_con[i].b_infilt) *
                                    (soil_con[i].max_moist[0] +
                                     soil_con[i].max_moist[1]);
        }

        // compute soil layer critical and wilting point moisture contents
        for (j = 0; j < options.Nlayer; j++) {
            soil_con[i].Wcr[j] *= soil_con[i].max_moist[j];
            soil_con[i].Wpwp[j] *= soil_con[i].max_moist[j];
            if (soil_con[i].Wpwp[j] > soil_con[i].Wcr[j]) {
                sprint_location(locstr, &(global_domain.locations[i]));
                sprintf(errstr, "Calculated wilting point moisture (%f mm) is "
                        "greater than calculated critical point moisture "
                        "(%f mm) for layer %zd."
                        "\n\tIn the soil parameter file, "
                        "Wpwp_FRACT MUST be <= Wcr_FRACT.\n%s",
                        soil_con[i].Wpwp[j], soil_con[i].Wcr[j], j, locstr);
                nrerror(errstr);
            }
            if (soil_con[i].Wpwp[j] < soil_con[i].resid_moist[j] *
                soil_con[i].depth[j] * 1000.) {
                sprint_location(locstr, &(global_domain.locations[i]));
                sprintf(errstr, "Calculated wilting point moisture (%f mm) is "
                        "less than calculated residual moisture (%f mm) for "
                        "layer %zd.\n\tIn the soil parameter file, "
                        "Wpwp_FRACT MUST be >= resid_moist / "
                        "(1.0 - bulk_density/soil_density).\n%s",
                        soil_con[i].Wpwp[j], soil_con[i].resid_moist[j] *
                        soil_con[i].depth[j] * 1000., j, locstr);
                nrerror(errstr);
            }
        }

        // validate spatial snow/frost params
        if (options.SPATIAL_SNOW) {
            if (soil_con[i].max_snow_distrib_slope < 0.0) {
                sprint_location(locstr, &(global_domain.locations[i]));
                sprintf(errstr, "max_snow_distrib_slope (%f) must be "
                        "positive.\n%s", soil_con[i].max_snow_distrib_slope,
                        locstr);
                nrerror(errstr);
            }
        }

        if (options.SPATIAL_FROST) {
            if (soil_con[i].frost_slope < 0.0) {
                sprint_location(locstr, &(global_domain.locations[i]));
                sprintf(errstr, "frost_slope (%f) must be positive.\n%s",
                        soil_con[i].frost_slope, locstr);
                nrerror(errstr);
            }
        }

        // If BASEFLOW = NIJSSEN2001 then convert NIJSSEN2001
        // parameters d1, d2, d3, and d4 to ARNO baseflow
        // parameters Ds, Dsmax, Ws, and c
        if (options.BASEFLOW == NIJSSEN2001) {
            j = options.Nlayer - 1;
            soil_con[i].Dsmax = soil_con[i].Dsmax *
                                pow((double)(1. /
                                             (soil_con[i].max_moist[j] -
                                              soil_con[i].Ws)),
                                    -soil_con[i].c) +
                                soil_con[i].Ds * soil_con[i].max_moist[j];
            soil_con[i].Ds = soil_con[i].Ds *
                             soil_con[i].Ws / soil_con[i].Dsmax;
            soil_con[i].Ws = soil_con[i].Ws / soil_con[i].max_moist[j];
        }

        soil_moisture_from_water_table(&(soil_con[i]), options.Nlayer);

        if (options.CARBON) {
            // TBD Remove hardcoded parameter values
            soil_con[i].AlbedoPar = 0.92 * BARE_SOIL_ALBEDO - 0.015;
            if (soil_con[i].AlbedoPar < AlbSoiParMin) {
                soil_con[i].AlbedoPar = AlbSoiParMin;
            }
        }
    }

    // read_snowband()
    if (options.SNOW_BAND == 1) {
        for (i = 0; i < global_domain.ncells_global; i++) {
            soil_con[i].AreaFract[0] = 1.;
            soil_con[i].BandElev[0] = soil_con[i].elevation;
            soil_con[i].Pfactor[0] = 1.;
            soil_con[i].Tfactor[0] = 0.;
        }
    }
    else {
        // AreaFract: fraction of grid cell in each snow band
        for (j = 0; j < options.SNOW_BAND; j++) {
            d3start[0] = j;
            get_nc_field_double(filenames.snowband, "AreaFract",
                                d3start, d3count, dvar);
            for (i = 0; i < global_domain.ncells_global; i++) {
                soil_con[i].AreaFract[j] = (double) dvar[idx[i]];
            }
        }
        // elevation: elevation of each snow band
        for (j = 0; j < options.SNOW_BAND; j++) {
            d3start[0] = j;
            get_nc_field_double(filenames.snowband, "elevation",
                                d3start, d3count, dvar);
            for (i = 0; i < global_domain.ncells_global; i++) {
                soil_con[i].BandElev[j] = (double) dvar[idx[i]];
            }
        }
        // Pfactor: precipitation multiplier for each snow band
        for (j = 0; j < options.SNOW_BAND; j++) {
            d3start[0] = j;
            get_nc_field_double(filenames.snowband, "Pfactor",
                                d3start, d3count, dvar);
            for (i = 0; i < global_domain.ncells_global; i++) {
                soil_con[i].BandElev[j] = (double) dvar[idx[i]];
            }
        }
        // Run some checks and corrections for soil
        for (i = 0; i < global_domain.ncells_global; i++) {
            // Make sure area fractions are positive and add to 1
            sum = 0.;
            for (j = 0; j < options.SNOW_BAND; j++) {
                if (soil_con[i].AreaFract[j] < 0) {
                    sprint_location(locstr, &(global_domain.locations[i]));
                    sprintf(errstr, "Negative snow band area fraction "
                            "(%f) read from file\n%s",
                            soil_con[i].AreaFract[j], locstr);
                    nrerror(errstr);
                }
                sum += soil_con[i].AreaFract[j];
            }
            // TBD: Need better check for equal to 1.
            if (sum != 1.) {
                sprint_location(locstr, &(global_domain.locations[i]));
                fprintf(stderr, "WARNING: Sum of the snow band area "
                        "fractions does not equal 1 (%f), dividing "
                        "each fraction by the sum\n%s", sum, locstr);
                for (j = 0; j < options.SNOW_BAND; j++) {
                    soil_con[i].AreaFract[j] /= sum;
                }
            }
            // check that the mean elevation from the snow bands matches the
            // grid cell mean elevation. If not reset mean
            mean = 0.;
            for (j = 0; j < options.SNOW_BAND; j++) {
                mean += soil_con[i].BandElev[j];
            }
            mean /= options.SNOW_BAND;
            if (fabs(soil_con[i].elevation - soil_con[i].BandElev[j]) > 1.0) {
                sprint_location(locstr, &(global_domain.locations[i]));
                fprintf(stderr, "WARNING: average band elevation %f not "
                        "equal to grid_cell average elevation %f; "
                        "setting grid cell elevation to average "
                        "band elevation.\n%s",
                        mean, soil_con[i].elevation, locstr);
                soil_con[i].elevation = (float) mean;
            }
            // Tfactor: calculate the temperature factor
            for (j = 0; j < options.SNOW_BAND; j++) {
                // TBD: Ensure that Tlapse is implemented consistently
                soil_con[i].Tfactor[j] = (soil_con[i].elevation -
                                          soil_con[i].BandElev[j]) * T_LAPSE;
            }
            // Pfactor: calculate Pfactor from the precipitation fraction read
            // from file
            // TBD: Ensure that netCDF variable is appropriately named
            sum = 0.;
            for (j = 0; j < options.SNOW_BAND; j++) {
                if (soil_con[i].Pfactor[j] < 0.) {
                    sprint_location(locstr, &(global_domain.locations[i]));
                    sprintf(errstr, "Snow band precipitation fraction (%lf) "
                            "must be between 0 and 1.\n%s",
                            soil_con[i].Pfactor[j], locstr);
                    nrerror(errstr);
                }
                if (soil_con[i].Pfactor[j] > 0. &&
                    soil_con[i].AreaFract[j] == 0) {
                    // TBD: Check to make sure whether this check is actually
                    // needed
                    sprint_location(locstr, &(global_domain.locations[i]));
                    sprintf(errstr, "Snow band precipitation fraction (%lf) "
                            "should be 0 when the area fraction is "
                            "0. (band = %zd).\n%s",
                            soil_con[i].AreaFract[j], j, locstr);
                    nrerror(errstr);
                }
                sum += soil_con[i].Pfactor[j];
            }
            // TBD: Need better check for equal to 1.
            if (sum != 1.) {
                sprint_location(locstr, &(global_domain.locations[i]));
                fprintf(stderr, "WARNING: Sum of the snow band precipitation "
                        "fractions does not equal 1 (%f), dividing "
                        "each fraction by the sum\n%s", sum, locstr);
                for (j = 0; j < options.SNOW_BAND; j++) {
                    soil_con[i].Pfactor[j] /= sum;
                }
            }
            // Pfactor: convert precipitation fraction to Pfactor
            for (j = 0; j < options.SNOW_BAND; j++) {
                if (soil_con[i].AreaFract[j] > 0) {
                    soil_con[i].Pfactor[j] /= soil_con[i].AreaFract[j];
                }
                else {
                    soil_con[i].Pfactor[j] = 0.;
                }
            }
        }
    }

    // logic from compute_treeline()
    for (i = 0; i < global_domain.ncells_global; i++) {
        for (j = 0; j < options.SNOW_BAND; j++) {
            // Lapse average annual July air temperature
            if (soil_con[i].avgJulyAirTemp + soil_con[i].Tfactor[j] <=
                TREELINE_TEMPERATURE) {
                // Snow band is above treeline
                soil_con[i].AboveTreeLine[j] = TRUE;
            }
            else {
                soil_con[i].AboveTreeLine[j] = FALSE;
            }
        }
    }

    // read_vegparam()

    // reading the vegetation parameters is slightly more complicated because
    // VIC allocates memory for veg_con only if the vegetation type exists in
    // the grid cell. The veg_con_map_struct is used to provide some of this
    // mapping

    // number of vegetation types - in vic this is defined without the bare soil
    // and the vegetation above the treeline
    for (i = 0; i < global_domain.ncells_global; i++) {
        nveg = veg_con_map[i].nv_active - 1;
        if (options.AboveTreelineVeg >= 0) {
            nveg -= 1;
        }
        for (j = 0; j < veg_con_map[i].nv_active; j++) {
            veg_con[i][j].vegetat_type_num = (int) nveg;
        }
    }

    // Cv: for each vegetation type, read the cover fraction into the mapping
    // structure. Then assign only the ones with a fraction greater than 0 to
    // the veg_con structure

    for (j = 0; j < options.NVEGTYPES; j++) {
        d3start[0] = j;
        get_nc_field_double(filenames.veg, "Cv",
                            d3start, d3count, dvar);
        for (i = 0; i < global_domain.ncells_global; i++) {
            veg_con_map[i].Cv[j] = (double) dvar[idx[i]];
        }
    }

    // do the mapping
    for (i = 0; i < global_domain.ncells_global; i++) {
        k = 0;
        for (j = 0; j < options.NVEGTYPES; j++) {
            if (veg_con_map[i].Cv[j] > 0) {
                veg_con_map[i].vidx[j] = k;
                veg_con[i][k].Cv = veg_con_map[i].Cv[j];
                veg_con[i][k].veg_class = j;
                k++;
            }
            else {
                veg_con_map[i].vidx[j] = -1;
            }
        }
    }

    // zone_depth: root zone depths
    for (j = 0; j < options.NVEGTYPES; j++) {
        d4start[0] = j;
        for (k = 0; k < options.ROOT_ZONES; k++) {
            d4start[1] = k;
            get_nc_field_double(filenames.veg, "root_depth",
                                d4start, d4count, dvar);
            for (i = 0; i < global_domain.ncells_global; i++) {
                vidx = veg_con_map[i].vidx[j];
                if (vidx != -1) {
                    veg_con[i][vidx].zone_depth[k] = (double) dvar[idx[i]];
                }
            }
        }
    }

    // zone_fract: root fractions
    for (j = 0; j < options.NVEGTYPES; j++) {
        d4start[0] = j;
        for (k = 0; k < options.ROOT_ZONES; k++) {
            d4start[1] = k;
            get_nc_field_double(filenames.veg, "root_fract",
                                d4start, d4count, dvar);
            for (i = 0; i < global_domain.ncells_global; i++) {
                vidx = veg_con_map[i].vidx[j];
                if (vidx != -1) {
                    veg_con[i][vidx].zone_fract[k] = (double) dvar[idx[i]];
                }
            }
        }
    }

    // calculate root fractions
    for (i = 0; i < global_domain.ncells_global; i++) {
        calc_root_fractions(veg_con[i], &(soil_con[i]));
    }

    // Run some checks and corrections for vegetation
    for (i = 0; i < global_domain.ncells_global; i++) {
        no_overstory = FALSE;
        // Only run to options.NVEGTYPES - 1, since bare soil is the last type
        for (j = 0; j < options.NVEGTYPES - 1; j++) {
            vidx = veg_con_map[i].vidx[j];
            if (vidx != -1) {
                sum = 0;
                for (k = 0; k < options.ROOT_ZONES; k++) {
                    sum += veg_con[i][vidx].zone_depth[k];
                }
                if (sum <= 0) {
                    sprint_location(locstr, &(global_domain.locations[i]));
                    sprintf(errstr,
                            "Root zone depths must sum to a value greater "
                            "than 0 (sum = %.2lf) - Type: %zd.\n%s", sum, j,
                            locstr);
                    nrerror(errstr);
                }
                sum = 0;
                for (k = 0; k < options.ROOT_ZONES; k++) {
                    sum += veg_con[i][vidx].zone_fract[k];
                }
                // TBD: Need better test for not equal to 1.
                if (sum != 1.) {
                    sprint_location(locstr, &(global_domain.locations[i]));
                    fprintf(stderr,
                            "WARNING: Root zone fractions sum to more than 1 "
                            "(%f), normalizing fractions.  If the sum is "
                            "large, check your vegetation parameter file.\n%s",
                            sum, locstr);
                    for (k = 0; k < options.ROOT_ZONES; k++) {
                        veg_con[i][vidx].zone_fract[k] /= sum;
                    }
                }
                // check that the vegetation type is defined in the vegetation
                // library
                found = FALSE;
                for (k = 0; k < options.NVEGTYPES; k++) {
                    if (veg_con[i][vidx].veg_class == veg_lib[i][k].veg_class) {
                        found = TRUE;
                        break;
                    }
                }
                if (!found) {
                    sprint_location(locstr, &(global_domain.locations[i]));
                    sprintf(errstr,
                            "The vegetation class id %i in vegetation tile %i "
                            "from cell %zd is not defined in the vegetation "
                            "library\n%s", veg_con[i][vidx].veg_class, vidx, i,
                            locstr);
                    nrerror(errstr);
                }
                // bad use of indexing -- Cv_sum should not be part of the
                // structure. Simply maintained for backward compatibility with
                // classic mode
                veg_con[i][0].Cv_sum += veg_con[i][vidx].Cv;

                // check for overstory
                if (!veg_lib[i][j].overstory) {
                    no_overstory = TRUE;
                }
            }
        }

        // handle the vegetation for the treeline option. This is somewhat
        // confusingly handled in VIC. If I am not mistaken, in VIC classic
        // this is handled in the following way:
        //
        // The treeline option is only active if there is more than one snow
        // band and options.COMPUTE_TREELINE is explicitly set in the global
        // file. If the treeline option is active, then there a few cases:
        //
        // 1. The grid cell contains one or more vegetation types that
        // do not have an overstory (either bare soil or vegetation). Nothing
        // further needs to be done to the input. For the elevation bands above
        // the treeline, the values from vegetation with an overstory are simply
        // ignored and the understory and bare ground values are scaled so they
        // cover the entire band. This scaling is done in put_data()
        //
        // 2. The grid cell contains only vegetation with an overstory.
        // In that case a small area of bare soil or vegetation without an
        // overstory must be created.  This will have almost no effect
        // on the results for most elevation bands, but above the treeline, the
        // elevation band will consists entirely of bare soil or the understory
        // vegetation (because of the scaling in put_data(). There are two
        // cases:
        //
        // 2.a. options.AboveTreelineVeg < 0. In that case a small amount of
        // bare soil is created (fraction is 0.001).
        //
        // 2.b. options.AboveTreelineVeg > 0. In that case a small amount of
        // the new vegetation is created (fraction is 0.001). This vegetation
        // should not have an overstory.
        //
        // The tricky parts are:
        //
        // Ensure that the correct number of vegetation types are reflected
        // for each cell.
        //
        // Ensure that bare soil remains the last vegetation type (the one with
        // the highest number). This will seem odd, but that is how it is
        // handled within VIC.
        //
        // Only case 2 needs to be handled explicitly

        if (options.SNOW_BAND > 1 && options.COMPUTE_TREELINE &&
            !no_overstory && veg_con[i][0].Cv_sum == 1.) {
            // Use bare soil above treeline
            if (options.AboveTreelineVeg < 0) {
                for (j = 0; j < options.NVEGTYPES; j++) {
                    vidx = veg_con_map[i].vidx[j];
                    if (vidx != -1) {
                        veg_con[i][vidx].Cv -=
                            0.001 / veg_con[i][vidx].vegetat_type_num;
                    }
                }
                veg_con[i][0].Cv_sum -= 0.001;
            }
            // Use defined vegetation type above treeline
            else {
                for (j = 0; j < options.NVEGTYPES; j++) {
                    vidx = veg_con_map[i].vidx[j];
                    if (vidx != -1) {
                        veg_con[i][vidx].Cv -=
                            0.001 / veg_con[i][vidx].vegetat_type_num;
                        veg_con[i][vidx].vegetat_type_num += 1;
                    }
                }
                veg_con[i][options.NVEGTYPES - 1].Cv = 0.001;
                veg_con[i][options.NVEGTYPES - 1].veg_class =
                    options.AboveTreelineVeg;
                veg_con[i][options.NVEGTYPES - 1].Cv_sum =
                    veg_con[i][0].Cv_sum;
                veg_con[i][options.NVEGTYPES - 1].vegetat_type_num =
                    veg_con[i][0].vegetat_type_num;
                // Since root zones are not defined they are copied from another
                // vegetation type.
                for (j = 0; j < options.ROOT_ZONES; j++) {
                    veg_con[i][options.NVEGTYPES - 1].zone_depth[j] =
                        veg_con[i][0].zone_depth[j];
                    veg_con[i][options.NVEGTYPES - 1].zone_fract[j] =
                        veg_con[i][0].zone_fract[j];
                }
                // redo the mapping to ensure that the veg type is active
                k = 0;
                for (j = 0; j < options.NVEGTYPES; j++) {
                    if (veg_con_map[i].Cv[j] > 0) {
                        veg_con_map[i].vidx[j] = k;
                        veg_con[i][k].Cv = veg_con_map[i].Cv[j];
                        veg_con[i][k].veg_class = j;
                        k++;
                    }
                    else {
                        veg_con_map[i].vidx[j] = -1;
                    }
                }
                // check that the vegetation type is defined in the vegetation
                // library
                found = FALSE;
                for (k = 0; k < options.NVEGTYPES; k++) {
                    if (veg_con[i][vidx].veg_class == veg_lib[i][k].veg_class) {
                        found = TRUE;
                        break;
                    }
                }
                if (!found) {
                    sprint_location(locstr, &(global_domain.locations[i]));
                    sprintf(errstr,
                            "The vegetation class id %i in vegetation tile %i "
                            "from cell %zd is not defined in the vegetation "
                            "library\n%s", veg_con[i][vidx].veg_class, vidx, i,
                            locstr);
                    nrerror(errstr);
                }
                // make sure it has no overstory
                veg_class = veg_con[i][options.NVEGTYPES - 1].veg_class;
                if (veg_lib[i][veg_class].overstory) {
                    sprint_location(locstr, &(global_domain.locations[i]));
                    sprintf(errstr,
                            "Vegetation class %i is defined to have overstory, "
                            "so it cannot be used as the default vegetation "
                            "type for above canopy snow bands.\n%s", veg_class,
                            locstr);
                    nrerror(errstr);
                }
            }
        }

        // Bare soil is now read in as the "last" (highest index) vegetation
        // class
        // rescale vegetation classes to 1.0 if their sum is greater than 0.99
        // otherwise throw an error
        // TBD: Need better check for equal to 1.
        if (veg_con[i][0].Cv_sum != 1.) {
            sprint_location(locstr, &(global_domain.locations[i]));
            fprintf(stderr,
                    "Cv !=  1.0 (%f) at grid cell %zd. Rescaling ...\n%s",
                    veg_con[i][0].Cv_sum, i, locstr);
            for (j = 0; j < options.NVEGTYPES; j++) {
                vidx = veg_con_map[i].vidx[j];
                if (vidx != -1) {
                    veg_con[i][vidx].Cv /= veg_con[i][0].Cv_sum;
                }
            }
            veg_con[i][0].Cv_sum = 1.;
        }
    }

    // TBD: implement the blowing snow option
    if (options.BLOWING) {
        sprintf(errstr,
                "BLOWING option not yet implemented in vic_init()");
        nrerror(errstr);
    }

    // read_lakeparam()
    // TBD: read lake parameters
    if (options.LAKES) {
        sprintf(errstr,
                "LAKES option not yet implemented in vic_init()");
        nrerror(errstr);
    }

    // initialize structures with default values
    for (i = 0; i < global_domain.ncells_global; i++) {
        nveg = veg_con[i][0].vegetat_type_num;
        initialize_snow(all_vars[i].snow, nveg);
        initialize_soil(all_vars[i].cell, &(soil_con[i]), veg_con[i], nveg);
        initialize_veg(all_vars[i].veg_var, veg_con[i], nveg);
        if (options.LAKES) {
            sprintf(errstr,
                    "LAKES option not yet implemented in vic_init()");
            nrerror(errstr);
        }
        initialize_energy(all_vars[i].energy, &(soil_con[i]), nveg);
    }

    // TBD: handle decomposed domain

    // cleanup
    free(dvar);
    free(idx);
}
