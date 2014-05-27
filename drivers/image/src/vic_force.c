#include <vic_def.h>
#include <vic_run.h>
#include <vic_driver_image.h>

void
vic_force()
{
    extern int                 NF;
    extern int                 NR;
    extern size_t              current;
    extern atmos_data_struct  *atmos;
    extern dmy_struct         *dmy;
    extern domain_struct       global_domain;
    extern filenames_struct    filenames;
    extern global_param_struct global_param;
    extern option_struct       options;

    double                     t_offset;
    float                     *fvar = NULL;
    size_t                     i;
    size_t                     j;
    size_t                    *idx = NULL;
    size_t                     d3count[3];
    size_t                     d3start[3];

    // allocate memory for variables to be read
    fvar = (float *) malloc(global_domain.n_ny * global_domain.n_nx *
                            sizeof(float));
    if (fvar == NULL) {
        nrerror("Memory allocation error in vic_force().");
    }

    // get 1D indices used in mapping the netcdf fields to the locations
    idx = (size_t *) malloc(global_domain.ncells_global *
                            sizeof(size_t));
    if (idx == NULL) {
        nrerror("Memory allocation error in vic_force().");
    }
    for (i = 0; i < global_domain.ncells_global; i++) {
        idx[i] = global_domain.locations[i].global_x_idx * global_domain.n_ny +
                 global_domain.locations[i].global_y_idx;
    }

    printf("%04d-%02d-%02d\n", dmy[current].year, dmy[current].month,
           dmy[current].day);

    // for now forcing file is determined by the year
    sprintf(filenames.forcing[0], "%s%4d.nc", filenames.f_path_pfx[0],
            dmy[current].year);

    // TBD: implement the correct offset for the data (d3start[0])
    d3start[1] = 0;
    d3start[2] = 0;
    d3count[0] = 1;
    d3count[1] = global_domain.n_ny;
    d3count[2] = global_domain.n_nx;

    // Air temperature: tas
    for (j = 0; j < NF; j++) {
        d3start[0] = current * NF + j;
        for (i = 0; i < global_domain.ncells_global; i++) {
            get_nc_field_float(filenames.forcing[0], "tas",
                               d3start, d3count, fvar);
            atmos[i].air_temp[j] = (double) fvar[idx[i]];
        }
    }

    // Precipitation: prcp
    for (j = 0; j < NF; j++) {
        d3start[0] = current * NF + j;
        for (i = 0; i < global_domain.ncells_global; i++) {
            get_nc_field_float(filenames.forcing[0], "prcp",
                               d3start, d3count, fvar);
            atmos[i].prec[j] = (double) fvar[idx[i]];
        }
    }

    // Downward solar radiation: dswrf
    for (j = 0; j < NF; j++) {
        d3start[0] = current * NF + j;
        for (i = 0; i < global_domain.ncells_global; i++) {
            get_nc_field_float(filenames.forcing[0], "dswrf",
                               d3start, d3count, fvar);
            atmos[i].shortwave[j] = (double) fvar[idx[i]];
        }
    }

    // Downward longwave radiation: dlwrf
    for (j = 0; j < NF; j++) {
        d3start[0] = current * NF + j;
        for (i = 0; i < global_domain.ncells_global; i++) {
            get_nc_field_float(filenames.forcing[0], "dlwrf",
                               d3start, d3count, fvar);
            atmos[i].longwave[j] = (double) fvar[idx[i]];
        }
    }

    // Wind speed: wind
    for (j = 0; j < NF; j++) {
        d3start[0] = current * NF + j;
        for (i = 0; i < global_domain.ncells_global; i++) {
            get_nc_field_float(filenames.forcing[0], "wind",
                               d3start, d3count, fvar);
            atmos[i].wind[j] = (double) fvar[idx[i]];
        }
    }

    // Specific humidity: shum
    for (j = 0; j < NF; j++) {
        d3start[0] = current * NF + j;
        for (i = 0; i < global_domain.ncells_global; i++) {
            get_nc_field_float(filenames.forcing[0], "shum",
                               d3start, d3count, fvar);
            atmos[i].vp[j] = (double) fvar[idx[i]];
        }
    }

    // Pressure: pressure
    for (j = 0; j < NF; j++) {
        d3start[0] = current * NF + j;
        for (i = 0; i < global_domain.ncells_global; i++) {
            get_nc_field_float(filenames.forcing[0], "pres",
                               d3start, d3count, fvar);
            atmos[i].pressure[j] = (double) fvar[idx[i]];
        }
    }

    if (options.SNOW_BAND > 1) {
        nrerror("SNOW_BAND not implemented in vic_force()");
    }
    else {
        t_offset = 0;
    }
    // Convert forcings into what we need and calculate missing ones
    for (i = 0; i < global_domain.ncells_global; i++) {
        for (j = 0; j < NF; j++) {
            // temperature in Kelvin
            atmos[i].air_temp[j] -= KELVIN;
            // precipitation in mm/period
            atmos[i].prec[j] *= (double) (options.SNOW_STEP * SECPHOUR);
            // pressure in kPa
            atmos[i].pressure[j] /= 1000.;
            // vapor pressure in kPa (we read specific humidity in kg/kg)
            atmos[i].vp[j] = q_to_vp(atmos[i].vp[j], atmos[i].pressure[j]);
            // vapor pressure deficit
            atmos[i].vpd[j] = svp(atmos[i].air_temp[j]) - atmos[i].vp[j];
            // photosynthetically active radiation
            atmos[i].par[j] = SW2PAR * atmos[i].shortwave[j];
            // air density
            atmos[i].density[j] = air_density(atmos[i].air_temp[j],
                                              atmos[i].pressure[j],
                                              atmos[i].vp[j]);
            // snow flag
            atmos[i].snowflag[j] = will_it_snow(&(atmos[i].air_temp[j]),
                                                t_offset,
                                                global_param.MAX_SNOW_TEMP,
                                                &(atmos[i].prec[j]), 1);
        }
    }


    // Put average value in NR field
    for (i = 0; i < global_domain.ncells_global; i++) {
        atmos[i].air_temp[NR] = average(atmos[i].air_temp, NF);
        // For precipitation put total
        atmos[i].prec[NR] = average(atmos[i].prec, NF) * NF;
        atmos[i].shortwave[NR] = average(atmos[i].shortwave, NF);
        atmos[i].longwave[NR] = average(atmos[i].longwave, NF);
        atmos[i].pressure[NR] = average(atmos[i].pressure, NF);
        atmos[i].wind[NR] = average(atmos[i].wind, NF);
        atmos[i].vp[NR] = average(atmos[i].vp, NF);
        atmos[i].vpd[NR] = (svp(atmos[i].air_temp[NR]) - atmos[i].vp[NR]);
        atmos[i].density[NR] = air_density(atmos[i].air_temp[NR],
                                           atmos[i].pressure[NR],
                                           atmos[i].vp[NR]);
        atmos[i].snowflag[NR] = will_it_snow(atmos[i].air_temp, t_offset,
                                             global_param.MAX_SNOW_TEMP,
                                             atmos[i].prec, NF);
    }

    // TBD: coszen (used for some of the carbon functions), fdir (if needed)
    // Catm, fdir (not used as far as I can tell)

    // cleanup
    free(fvar);
    free(idx);
}

double
average(double *ar,
        size_t  n)
{
    size_t i;
    double sum = 0.;

    if (n <= 0) {
        nrerror("Error in calc_average: divide by zero or negative");
    }
    else if (n == 1) {
        return ar[0];
    }
    else {
        for (i = 0; i < n; i++) {
            sum += ar[i];
        }
    }

    return sum / n;
}

// convert specific humidity (q) to vapor pressure (vp) based on pressure (p)
// returned units are those of p
double
q_to_vp(double q,
        double p)
{
    double vp;

    // full equation
    // vp = q/(q+EPS*(1-q))*p;

    // approximation used in VIC
    vp = q * p / EPS;

    return vp;
}

// convert surface pressure (kPa) to density (kg/m3) based on pressure (p),
// vapor pressure (vp), and temperature
double
air_density(double t,
            double p,
            double vp)
{
    double rho;

    // full equation
    // rho = (p*1000)/(Rd * *t+KELVIN) + (pv*1000)/(Rv * *t+KELVIN);

    // approximation used in VIC
    rho = 0.003486 * p / (275.0 + t);

    return rho;
}

char
will_it_snow(double *t,
             double  t_offset,
             double  max_snow_temp,
             double *prcp,
             int     n)
{
    size_t i;

    for (i = 0; i < n; i++) {
        if ((t[i] + t_offset) < max_snow_temp && prcp[i] > 0.) {
            return 1;
        }
    }

    return 0;
}
