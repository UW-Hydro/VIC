/******************************************************************************
 * @section DESCRIPTION
 *
 * This subroutine redistributes soil properties based on the thermal solutions
 * found for the current time step.
 *
 * @section LICENSE
 *
 * The Variable Infiltration Capacity (VIC) macroscale hydrological model
 * Copyright (C) 2016 The Computational Hydrology Group, Department of Civil
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

#include <vic_run.h>

/******************************************************************************
 * @brief    This subroutine redistributes soil properties based on the thermal
 *           solutions found for the current time step.
 *****************************************************************************/
int
calc_layer_average_thermal_props(energy_bal_struct *energy,
                                 layer_data_struct *layer,
                                 soil_con_struct   *soil_con,
                                 size_t             Nnodes,
                                 double            *T)
{
    extern option_struct options;

    size_t               i;
    int                  ErrorFlag;
    size_t               tmpTshape[] = {
        options.Nlayer, Nnodes,
        options.Nfrost + 1
    };
    size_t               tmpZshape[] = {
        options.Nlayer, Nnodes
    };
    double            ***tmpT;
    double             **tmpZ;

    // allocate memory for tmpT and tmpZ
    malloc_3d_double(tmpTshape, &tmpT);
    malloc_2d_double(tmpZshape, &tmpZ);

    if (options.FROZEN_SOIL && soil_con->FS_ACTIVE) {
        find_0_degree_fronts(energy, soil_con->Zsum_node, T, Nnodes);
    }
    else {
        energy->Nfrost = 0;
    }

    /** Store Layer Temperature Values **/
    for (i = 0; i < Nnodes; i++) {
        energy->T[i] = T[i];
    }

    if (energy->Nfrost > 0) {
        energy->frozen = true;
    }
    else {
        energy->frozen = false;
    }

    /** Compute Soil Layer average  properties **/
    if (options.QUICK_FLUX) {
        ErrorFlag = estimate_layer_temperature_quick_flux(layer,
                                                          soil_con->depth,
                                                          soil_con->dp,
                                                          energy->T[0],
                                                          energy->T[1],
                                                          soil_con->avg_temp);
        if (ErrorFlag == ERROR) {
            return (ERROR);
        }
        ErrorFlag = estimate_layer_ice_content_quick_flux(layer,
                                                          soil_con->depth,
                                                          soil_con->max_moist,
                                                          soil_con->expt,
                                                          soil_con->bubble,
                                                          soil_con->frost_fract,
                                                          soil_con->frost_slope,
                                                          soil_con->FS_ACTIVE);
        if (ErrorFlag == ERROR) {
            return (ERROR);
        }
    }
    else {
        estimate_frost_temperature_and_depth(tmpT,
                                             tmpZ,
                                             soil_con->Zsum_node,
                                             energy->T,
                                             soil_con->depth,
                                             soil_con->frost_fract,
                                             soil_con->frost_slope,
                                             Nnodes,
                                             options.Nlayer);
        ErrorFlag = estimate_layer_temperature(layer,
                                               tmpT,
                                               tmpZ,
                                               soil_con->Zsum_node,
                                               soil_con->depth,
                                               Nnodes,
                                               options.Nlayer);
        if (ErrorFlag == ERROR) {
            return (ERROR);
        }
        ErrorFlag = estimate_layer_ice_content(layer,
                                               tmpT,
                                               tmpZ,
                                               soil_con->Zsum_node,
                                               soil_con->depth,
                                               soil_con->max_moist,
                                               soil_con->expt,
                                               soil_con->bubble,
                                               Nnodes,
                                               options.Nlayer,
                                               soil_con->FS_ACTIVE);
        if (ErrorFlag == ERROR) {
            return (ERROR);
        }
    }

    // free memory for tmpT and tmpZ
    free_3d_double(tmpTshape, tmpT);
    free_2d_double(tmpZshape, tmpZ);

    return (0);
}

/******************************************************************************
 * @brief    Iteratively solve the soil temperature profile using a numerical
 *           difference equation.  The solution equation is second order in
 *           space, and first order in time.
 *****************************************************************************/
int
solve_T_profile(double   *T,
                double   *T0,
                char     *Tfbflag,
                unsigned *Tfbcount,
                double   *Zsum,
                double   *kappa,
                double   *Cs,
                double   *moist,
                double    deltat,
                double   *max_moist,
                double   *bubble,
                double   *expt,
                double   *ice,
                double   *alpha,
                double   *beta,
                double   *gamma,
                double    Dp,
                int       Nnodes,
                int      *FIRST_SOLN,
                int       FS_ACTIVE,
                int       NOFLUX,
                int       EXP_TRANS)
{
    double *aa, *bb, *cc, *dd, *ee, Bexp;
    int     Error;
    int     j;

    // TODO: remove use of static variables (see GH #735), for now:
    // make static variables thread safe
    static double A[MAX_NODES];
    static double B[MAX_NODES];
    static double C[MAX_NODES];
    static double D[MAX_NODES];
    static double E[MAX_NODES];
    #pragma omp threadprivate(A, B, C, D, E)

    if (FIRST_SOLN[0]) {
        if (EXP_TRANS) {
            Bexp = logf(Dp + 1.) / (double) (Nnodes - 1);
        }

        FIRST_SOLN[0] = false;
        if (!EXP_TRANS) {
            for (j = 1; j < Nnodes - 1; j++) {
                A[j] = Cs[j] * alpha[j - 1] * alpha[j - 1];
                B[j] = (kappa[j + 1] - kappa[j - 1]) * deltat;
                C[j] = 2 * deltat * kappa[j] * alpha[j - 1] / gamma[j - 1];
                D[j] = 2 * deltat * kappa[j] * alpha[j - 1] / beta[j - 1];
                E[j] = CONST_RHOICE * CONST_LATICE *
                       alpha[j - 1] * alpha[j - 1];
            }
            if (NOFLUX) {
                j = Nnodes - 1;
                A[j] = Cs[j] * alpha[j - 1] * alpha[j - 1];
                B[j] = (kappa[j] - kappa[j - 1]) * deltat;
                C[j] = 2 * deltat * kappa[j] * alpha[j - 1] / gamma[j - 1];
                D[j] = 2 * deltat * kappa[j] * alpha[j - 1] / beta[j - 1];
                E[j] = CONST_RHOICE * CONST_LATICE *
                       alpha[j - 1] * alpha[j - 1];
            }
        }
        else { // grid transformation terms
            for (j = 1; j < Nnodes - 1; j++) {
                A[j] = 4 * Bexp * Bexp * Cs[j] * (Zsum[j] + 1) * (Zsum[j] + 1);
                B[j] = (kappa[j + 1] - kappa[j - 1]) * deltat;
                C[j] = 4 * deltat * kappa[j];
                D[j] = 2 * deltat * kappa[j] * Bexp;
                E[j] = 4 * Bexp * Bexp * CONST_RHOICE * CONST_LATICE *
                       (Zsum[j] + 1) * (Zsum[j] + 1);
            }
            if (NOFLUX) {
                j = Nnodes - 1;
                A[j] = 4 * Bexp * Bexp * Cs[j] * (Zsum[j] + 1) * (Zsum[j] + 1);
                B[j] = (kappa[j] - kappa[j - 1]) * deltat;
                C[j] = 4 * deltat * kappa[j];
                D[j] = 2 * deltat * kappa[j] * Bexp;
                E[j] = 4 * Bexp * Bexp * CONST_RHOICE * CONST_LATICE *
                       (Zsum[j] + 1) * (Zsum[j] + 1);
            }
        }
    }

    aa = &A[0];
    bb = &B[0];
    cc = &C[0];
    dd = &D[0];
    ee = &E[0];

    for (j = 0; j < Nnodes; j++) {
        T[j] = T0[j];
    }

    Error = calc_soil_thermal_fluxes(Nnodes, T, T0, Tfbflag, Tfbcount, moist,
                                     max_moist, ice, bubble, expt, gamma, aa,
                                     bb, cc, dd, ee, FS_ACTIVE, NOFLUX,
                                     EXP_TRANS);

    return (Error);
}

/******************************************************************************
 * @brief    Iteratively solve the soil temperature profile using a numerical
 *           difference equation.  The solution equation is second order in
 *           space, and first order in time.
 *****************************************************************************/
int
solve_T_profile_implicit(double   *T,                               // update
                         double   *T0,                        // keep
                         char     *Tfbflag,
                         unsigned *Tfbcount,
                         double   *Zsum,                      // soil parameter
                         double   *kappa,                     // update if necessary
                         double   *Cs,                        // update if necessary
                         double   *moist,                     // keep
                         double    deltat,                    // model parameter
                         double   *max_moist,                 // soil parameter
                         double   *bubble,                    // soil parameter
                         double   *expt,                      // soil parameter
                         double   *ice,                       // update if necessary
                         double   *alpha,                     // soil parameter
                         double   *beta,                      // soil parameter
                         double   *gamma,                     // soil parameter
                         double    Dp,                        // soil parameter
                         int       Nnodes,                   // model parameter
                         int      *FIRST_SOLN,               // update
                         int       NOFLUX,
                         int       EXP_TRANS,
                         double   *bulk_dens_min,              // soil parameter
                         double   *soil_dens_min,              // soil parameter
                         double   *quartz,                    // soil parameter
                         double   *bulk_density,              // soil parameter
                         double   *soil_density,              // soil parameter
                         double   *organic,                    // soil parameter
                         double   *depth)                     // soil parameter
{
    extern option_struct options;
    int                  n, Error;
    double               res[MAX_NODES];
    void                 (*vecfunc)(double *, double *, int, int, ...);
    int                  j;

    if (FIRST_SOLN[0]) {
        FIRST_SOLN[0] = false;
    }

    // initialize fda_heat_eqn:
    // pass model parameters, initial states, and soil parameters
    // it MUST be initialized before Newton-Raphson searching
    if (!NOFLUX) {
        n = Nnodes - 2;
    }
    else {
        n = Nnodes - 1;
    }

    fda_heat_eqn(&T[1], res, n, 1, deltat, NOFLUX, EXP_TRANS, T0,
                 moist, ice, kappa, Cs, max_moist, bubble, expt,
                 alpha, beta, gamma, Zsum, Dp, bulk_dens_min, soil_dens_min,
                 quartz, bulk_density, soil_density, organic, depth,
                 options.Nlayer);

    // modified Newton-Raphson to solve for new T
    vecfunc = &(fda_heat_eqn);
    Error = newt_raph(vecfunc, &T[1], n);

    // update temperature boundaries
    if (Error == 0) {
        T[0] = T0[0]; // surface
        if (!NOFLUX) {
            T[Nnodes - 1] = T0[Nnodes - 1]; // bottom boundary
        }
        if (options.TFALLBACK) {
            // HACK to prevent runaway cold nose
            // Handle the case in which the a node was colder than both the nodes above and below
            // in the last time step, and that both differences have increased between the last
            // time step and the current one.
            for (j = 1; j < Nnodes - 1; j++) {
                if ((T0[j - 1] - T0[j] > 0 && T0[j + 1] - T0[j] > 0 &&
                     (T[j - 1] - T[j]) - (T0[j - 1] - T0[j]) > 0 &&
                     (T[j + 1] - T[j]) - (T0[j + 1] - T0[j]) > 0) ||
                    (T0[j - 1] - T0[j] < 0 && T0[j + 1] - T0[j] < 0 &&
                     (T[j - 1] - T[j]) - (T0[j - 1] - T0[j]) < 0 &&
                     (T[j + 1] - T[j]) - (T0[j + 1] - T0[j]) < 0)) {
                    T[j] = 0.5 * (T[j - 1] + T[j + 1]); // crude fix for now; just average the T's without taking distance, conductivities into account
                    Tfbflag[j] = 1;
                    Tfbcount[j]++;
                }
            }
        }
    }

    return (Error);
}

/******************************************************************************
 * @brief    Calculate soil thermal fluxes
 *****************************************************************************/
int
calc_soil_thermal_fluxes(int       Nnodes,
                         double   *T,
                         double   *T0,
                         char     *Tfbflag,
                         unsigned *Tfbcount,
                         double   *moist,
                         double   *max_moist,
                         double   *ice,
                         double   *bubble,
                         double   *expt,
                         double   *gamma,
                         double   *A,
                         double   *B,
                         double   *C,
                         double   *D,
                         double   *E,
                         int       FS_ACTIVE,
                         int       NOFLUX,
                         int       EXP_TRANS)
{
    extern option_struct     options;
    extern parameters_struct param;

    int                      Error;
    char                     Done;
    int                      j;
    int                      ItCount;
    double                   threshold = 1.e-2; /* temperature profile iteration threshold */
    double                   maxdiff;
    double                   diff;
    double                   oldT;
    double                   Tlast[MAX_NODES];

    Error = 0;
    Done = false;
    ItCount = 0;

    /* initialize Tlast */
    for (j = 0; j < Nnodes; j++) {
        Tlast[j] = T[j];
    }

    /* initialize Tfbflag, Tfbcount */
    for (j = 0; j < Nnodes; j++) {
        Tfbflag[j] = 0;
        Tfbcount[j] = 0;
    }

    while (!Done && Error == 0 && ItCount < param.FROZEN_MAXITER) {
        ItCount++;
        maxdiff = threshold;
        for (j = 1; j < Nnodes - 1; j++) {
            oldT = T[j];

            /**	2nd order variable kappa equation **/

            if (T[j] >= 0 || !FS_ACTIVE || !options.FROZEN_SOIL) {
                if (!EXP_TRANS) {
                    T[j] = (A[j] * T0[j] +
                            B[j] * (T[j + 1] - T[j - 1]) +
                            C[j] * T[j + 1] +
                            D[j] * T[j - 1] +
                            E[j] * (0. - ice[j])) / (A[j] + C[j] + D[j]);
                }
                else {
                    T[j] = (A[j] * T0[j] +
                            B[j] * (T[j + 1] - T[j - 1]) +
                            C[j] * (T[j + 1] + T[j - 1]) -
                            D[j] * (T[j + 1] - T[j - 1]) +
                            E[j] * (0. - ice[j])) / (A[j] + 2. * C[j]);
                }
            }
            else {
                T[j] =
                    root_brent(T0[j] - (param.SOIL_DT), T0[j] + (param.SOIL_DT),
                               soil_thermal_eqn,
                               T[j + 1], T[j - 1], T0[j], moist[j],
                               max_moist[j], bubble[j], expt[j], ice[j],
                               A[j], B[j], C[j], D[j], E[j], EXP_TRANS, j);
                if (T[j] <= -998) {
                    if (options.TFALLBACK) {
                        T[j] = T0[j];
                        Tfbflag[j] = 1;
                        Tfbcount[j]++;
                    }
                    else {
                        error_solve_T_profile(T[j], T[j + 1], T[j - 1], T0[j],
                                              moist[j],
                                              max_moist[j], bubble[j], expt[j],
                                              ice[j],
                                              gamma[j - 1], A[j], B[j], C[j],
                                              D[j],
                                              E[j]);
                        return (ERROR);
                    }
                }
            }

            diff = fabs(oldT - T[j]);
            if (diff > maxdiff) {
                maxdiff = diff;
            }
        }

        if (NOFLUX) {
            /** Solve for bottom temperature if using no flux lower boundary **/
            j = Nnodes - 1;
            oldT = T[j];

            if (T[j] >= 0 || !FS_ACTIVE || !options.FROZEN_SOIL) {
                if (!EXP_TRANS) {
                    T[j] =
                        (A[j] * T0[j] + B[j] * (T[j] - T[j - 1]) +
                         C[j] * T[j] + D[j] * T[j - 1] +
                         E[j] * (0. - ice[j])) / (A[j] + C[j] + D[j]);
                }
                else {
                    T[j] = (A[j] * T0[j] +
                            B[j] * (T[j] - T[j - 1]) +
                            C[j] * (T[j] + T[j - 1]) -
                            D[j] * (T[j] - T[j - 1]) +
                            E[j] * (0. - ice[j])) / (A[j] + 2. * C[j]);
                }
            }
            else {
                T[Nnodes - 1] = root_brent(T0[Nnodes - 1] - param.SOIL_DT,
                                           T0[Nnodes - 1] + param.SOIL_DT,
                                           soil_thermal_eqn,
                                           T[Nnodes - 1],
                                           T[Nnodes - 2], T0[Nnodes - 1],
                                           moist[Nnodes - 1],
                                           max_moist[Nnodes - 1],
                                           bubble[j], expt[Nnodes - 1],
                                           ice[Nnodes - 1],
                                           A[j], B[j], C[j], D[j], E[j],
                                           EXP_TRANS, j);
                if (T[j] <= -998) {
                    if (options.TFALLBACK) {
                        T[j] = T0[j];
                        Tfbflag[j] = 1;
                        Tfbcount[j]++;
                    }
                    else {
                        error_solve_T_profile(T[Nnodes - 1], T[Nnodes - 1],
                                              T[Nnodes - 2], T0[Nnodes - 1],
                                              moist[Nnodes - 1],
                                              max_moist[Nnodes - 1],
                                              bubble[Nnodes - 1],
                                              expt[Nnodes - 1], ice[Nnodes - 1],
                                              gamma[Nnodes - 2],
                                              A[j], B[j], C[j], D[j], E[j]);
                        return (ERROR);
                    }
                }
            }

            diff = fabs(oldT - T[Nnodes - 1]);
            if (diff > maxdiff) {
                maxdiff = diff;
            }
        }

        if (maxdiff <= threshold) {
            Done = true;
        }
    }

    if (options.TFALLBACK) {
        // HACK to prevent runaway cold nose
        // Handle the case in which the a node was colder than both the nodes above and below
        // in the last time step, and that both differences have increased between the last
        // time step and the current one.
        for (j = 1; j < Nnodes - 1; j++) {
            if ((Tlast[j - 1] - Tlast[j] > 0 && Tlast[j + 1] - Tlast[j] > 0 &&
                 (T[j - 1] - T[j]) - (Tlast[j - 1] - Tlast[j]) > 0 &&
                 (T[j + 1] - T[j]) - (Tlast[j + 1] - Tlast[j]) > 0) ||
                (Tlast[j - 1] - Tlast[j] < 0 && Tlast[j + 1] - Tlast[j] < 0 &&
                 (T[j - 1] - T[j]) - (Tlast[j - 1] - Tlast[j]) < 0 &&
                 (T[j + 1] - T[j]) - (Tlast[j + 1] - Tlast[j]) < 0)) {
                T[j] = 0.5 * (T[j - 1] + T[j + 1]); // crude fix for now; just average the T's without taking distance, conductivities into account
                Tfbflag[j] = 1;
                Tfbcount[j]++;
            }
        }
    }

    if (!Done && !Error) {
        if (options.TFALLBACK) {
            for (j = 0; j < Nnodes; j++) {
                T[j] = T0[j];
                Tfbflag[j] = 1;
                Tfbcount[j]++;
            }
        }
        else {
            fprintf(LOG_DEST,
                    "ERROR: Temperature Profile Unable to Converge!!!\n");
            fprintf(LOG_DEST, "Dumping Profile Temperatures (last, new).\n");
            for (j = 0; j < Nnodes; j++) {
                fprintf(LOG_DEST, "%f\t%f\n", T0[j], T[j]);
            }
            log_err("Cannot solve temperature profile:\n"
                    "\tToo Many Iterations in solve_T_profile");
            return (ERROR);
        }
    }

    return (Error);
}

/******************************************************************************
 * @brief    Dummy function to allow calling of error_print_solve_T_profile()
 *****************************************************************************/
double
error_solve_T_profile(double Tj,
                      ...)
{
    va_list ap;

    double  error;

    va_start(ap, Tj);
    error = error_print_solve_T_profile(Tj, ap);
    va_end(ap);

    return error;
}

/******************************************************************************
 * @brief    Print soil temperature terms.
 *****************************************************************************/
double
error_print_solve_T_profile(double  T,
                            va_list ap)
{
    double TL;
    double TU;
    double T0;
    double moist;
    double max_moist;
    double bubble;
    double expt;
    double ice0;
    double gamma;
    double A;
    double B;
    double C;
    double D;
    double E;

    TL = (double) va_arg(ap, double);
    TU = (double) va_arg(ap, double);
    T0 = (double) va_arg(ap, double);
    moist = (double) va_arg(ap, double);
    max_moist = (double) va_arg(ap, double);
    bubble = (double) va_arg(ap, double);
    expt = (double) va_arg(ap, double);
    ice0 = (double) va_arg(ap, double);
    gamma = (double) va_arg(ap, double);
    A = (double) va_arg(ap, double);
    B = (double) va_arg(ap, double);
    C = (double) va_arg(ap, double);
    D = (double) va_arg(ap, double);
    E = (double) va_arg(ap, double);

    log_warn("solve_T_profile failed to converge to a solution "
             "in root_brent.  Variable values will be dumped to the "
             "screen, check for invalid values.");

    fprintf(LOG_DEST, "T\t%f\n", T);
    fprintf(LOG_DEST, "TL\t%f\n", TL);
    fprintf(LOG_DEST, "TU\t%f\n", TU);
    fprintf(LOG_DEST, "T0\t%f\n", T0);
    fprintf(LOG_DEST, "moist\t%f\n", moist);
    fprintf(LOG_DEST, "max_moist\t%f\n", max_moist);
    fprintf(LOG_DEST, "bubble\t%f\n", bubble);
    fprintf(LOG_DEST, "expt\t%f\n", expt);
    fprintf(LOG_DEST, "ice0\t%f\n", ice0);
    fprintf(LOG_DEST, "gamma\t%f\n", gamma);
    fprintf(LOG_DEST, "A\t%f\n", A);
    fprintf(LOG_DEST, "B\t%f\n", B);
    fprintf(LOG_DEST, "C\t%f\n", C);
    fprintf(LOG_DEST, "D\t%f\n", D);
    fprintf(LOG_DEST, "E\t%f\n", E);

    log_warn("Finished dumping values for solve_T_profile.\n"
             "Try increasing SOIL_DT to get model to complete cell.\n"
             "Then check output for instabilities.");

    return(ERROR);
}

/******************************************************************************
 * @brief    Heat Equation for implicit scheme (used to calculate residual of
 *           the heat equation) passed from solve_T_profile_implicit
 *****************************************************************************/
void
fda_heat_eqn(double T_2[],
             double res[],
             int    n,
             int    init,
             ...)
{
    char    PAST_BOTTOM;
    double  storage_term, flux_term, phase_term, flux_term1, flux_term2;
    double  Lsum;
    int     i;
    size_t  lidx;
    int     focus, left, right;

    // argument list handling
    va_list arg_addr;

    // TODO: remove use of static variables (see GH #735), for now:
    // make static variables thread safe
    static double  deltat;
    static int     NOFLUX;
    static int     EXP_TRANS;
    static double *T0;
    static double *moist;
    static double *ice;
    static double *kappa;
    static double *Cs;
    static double *max_moist;
    static double *bubble;
    static double *expt;
    static double *alpha;
    static double *beta;
    static double *gamma;
    static double *Zsum;
    static double  Dp;
    static double *bulk_dens_min;
    static double *soil_dens_min;
    static double *quartz;
    static double *bulk_density;
    static double *soil_density;
    static double *organic;
    static double *depth;
    static size_t  Nlayers;

    // variables used to calculate residual of the heat equation
    // defined here
    static double Ts;
    static double Tb;

    // locally used variables
    static double ice_new[MAX_NODES], Cs_new[MAX_NODES], kappa_new[MAX_NODES];
    static double DT[MAX_NODES], DT_down[MAX_NODES], DT_up[MAX_NODES];
    static double Dkappa[MAX_NODES];
    static double Bexp;

    #pragma omp threadprivate(deltat, NOFLUX, EXP_TRANS, T0, moist, ice, \
    kappa, Cs, max_moist, bubble, expt, alpha, beta, gamma, Zsum, Dp, \
    bulk_dens_min, soil_dens_min, quartz, bulk_density, soil_density, organic, \
    depth, Nlayers, Ts, Tb, ice_new, Cs_new, kappa_new, DT, DT_down, DT_up, \
    Dkappa, Bexp)

    // initialize variables if init==1
    if (init == 1) {
        va_start(arg_addr, init);
        deltat = va_arg(arg_addr, double);
        NOFLUX = va_arg(arg_addr, int);
        EXP_TRANS = va_arg(arg_addr, int);
        T0 = va_arg(arg_addr, double *);
        moist = va_arg(arg_addr, double *);
        ice = va_arg(arg_addr, double *);
        kappa = va_arg(arg_addr, double *);
        Cs = va_arg(arg_addr, double *);
        max_moist = va_arg(arg_addr, double *);
        bubble = va_arg(arg_addr, double *);
        expt = va_arg(arg_addr, double *);
        alpha = va_arg(arg_addr, double *);
        beta = va_arg(arg_addr, double *);
        gamma = va_arg(arg_addr, double *);
        Zsum = va_arg(arg_addr, double *);
        Dp = va_arg(arg_addr, double);
        bulk_dens_min = va_arg(arg_addr, double *);
        soil_dens_min = va_arg(arg_addr, double *);
        quartz = va_arg(arg_addr, double *);
        bulk_density = va_arg(arg_addr, double *);
        soil_density = va_arg(arg_addr, double *);
        organic = va_arg(arg_addr, double *);
        depth = va_arg(arg_addr, double *);
        Nlayers = va_arg(arg_addr, size_t);

        if (EXP_TRANS) {
            if (!NOFLUX) {
                Bexp = logf(Dp + 1.) / (double)(n + 1);
            }
            else {
                Bexp = logf(Dp + 1.) / (double)(n);
            }
        }

        Ts = T0[0];
        if (!NOFLUX) {
            Tb = T0[n + 1];
        }
        else {
            Tb = T0[n];
        }
        for (i = 0; i < n; i++) {
            T_2[i] = T0[i + 1];
        }
    }
    // calculate residuals if init==0
    else {
        // get the range of columns to calculate
        va_start(arg_addr, init);
        focus = va_arg(arg_addr, int);

        // calculate all entries if focus == -1
        if (focus == -1) {
            lidx = 0;
            Lsum = 0.;
            PAST_BOTTOM = false;

            for (i = 0; i < n + 1; i++) {
                kappa_new[i] = kappa[i];
                if (i >= 1) { // all but surface node
                    // update ice contents
                    if (T_2[i - 1] < 0) {
                        ice_new[i] = moist[i] - maximum_unfrozen_water(
                            T_2[i - 1],
                            max_moist[
                                i], bubble[i], expt[i]);
                        if (ice_new[i] < 0) {
                            ice_new[i] = 0;
                        }
                    }
                    else {
                        ice_new[i] = 0;
                    }
                    Cs_new[i] = Cs[i];

                    // update other states due to ice content change
                    /***********************************************/
                    if (ice_new[i] != ice[i]) {
                        kappa_new[i] = soil_conductivity(moist[i],
                                                         moist[i] - ice_new[i],
                                                         soil_dens_min[lidx],
                                                         bulk_dens_min[lidx],
                                                         quartz[lidx],
                                                         soil_density[lidx],
                                                         bulk_density[lidx],
                                                         organic[lidx]);
                        Cs_new[i] = volumetric_heat_capacity(
                            bulk_density[lidx] / soil_density[lidx],
                            moist[i] - ice_new[i], ice_new[i], organic[lidx]);
                    }
                    /************************************************/
                }

                if (Zsum[i] > Lsum + depth[lidx] && !PAST_BOTTOM) {
                    Lsum += depth[lidx];
                    lidx++;
                    if (lidx == Nlayers) {
                        PAST_BOTTOM = true;
                        lidx = Nlayers - 1;
                    }
                }
            }

            // constants used in fda equation
            for (i = 0; i < n; i++) {
                if (i == 0) {
                    DT[i] = T_2[i + 1] - Ts;
                    DT_up[i] = T_2[i] - Ts;
                    DT_down[i] = T_2[i + 1] - T_2[i];
                }
                else if (i == n - 1) {
                    DT[i] = Tb - T_2[i - 1];
                    DT_up[i] = T_2[i] - T_2[i - 1];
                    DT_down[i] = Tb - T_2[i];
                }
                else {
                    DT[i] = T_2[i + 1] - T_2[i - 1];
                    DT_up[i] = T_2[i] - T_2[i - 1];
                    DT_down[i] = T_2[i + 1] - T_2[i];
                }
                if (i < n - 1) {
                    Dkappa[i] = kappa_new[i + 2] - kappa_new[i];
                }
                else if (!NOFLUX) {
                    Dkappa[i] = kappa_new[i + 2] - kappa_new[i];
                }
                else {
                    Dkappa[i] = kappa_new[i + 1] - kappa_new[i];
                }
            }

            for (i = 0; i < n; i++) {
                storage_term =
                    Cs_new[i +
                           1] *
                    (T_2[i] -
                     T0[i +
                        1]) / deltat + T_2[i] *
                    (Cs_new[i + 1] - Cs[i + 1]) / deltat;
                if (!EXP_TRANS) {
                    flux_term1 = Dkappa[i] / alpha[i] * DT[i] / alpha[i];
                    flux_term2 =
                        kappa_new[i +
                                  1] *
                        (DT_down[i] / gamma[i] - DT_up[i] /
                         beta[i]) / (0.5 * alpha[i]);
                }
                else { // grid transformation
                    flux_term1 = Dkappa[i] / 2. * DT[i] / 2. /
                                 (Bexp *
                                  (Zsum[i +
                                        1] + 1.)) / (Bexp * (Zsum[i + 1] + 1.));
                    flux_term2 =
                        kappa_new[i +
                                  1] *
                        ((DT_down[i] -
                          DT_up[i]) /
                         (Bexp *
                          (Zsum[i +
                                1] +
                           1.)) /
                         (Bexp *
                          (Zsum[i +
                                1] +
                           1.)) - DT[i] / 2. /
                         (Bexp * (Zsum[i + 1] + 1.) * (Zsum[i + 1] + 1.)));
                }
                // inelegant fix for "cold nose" problem - when a very cold node skates off to
                // much colder and breaks the second law of thermodynamics (because
                // flux_term1 exceeds flux_term2 in absolute magnitude) - therefore, don't let
                // that node get any colder.  This only seems to happen in the first and
                // second near-surface nodes.
                flux_term = flux_term1 + flux_term2;
                phase_term = CONST_RHOICE * CONST_LATICE *
                             (ice_new[i + 1] - ice[i + 1]) / deltat;
                res[i] = flux_term + phase_term - storage_term;
            }
        }
        // only calculate entries focus-1, focus, and focus+1 if focus has a value>=0
        else {
            if (focus == 0) {
                left = 0;
            }
            else {
                left = focus - 1;
            }
            if (focus == n - 1) {
                right = n - 1;
            }
            else {
                right = focus + 1;
            }

            // update ice content for node focus and its adjacents
            for (i = left; i <= right; i++) {
                if (T_2[i] < 0) {
                    ice_new[i + 1] = moist[i + 1] - maximum_unfrozen_water(
                        T_2[i],
                        max_moist[
                            i + 1], bubble[i + 1], expt[i + 1]);
                    if (ice_new[i + 1] < 0) {
                        ice_new[i + 1] = 0;
                    }
                }
                else {
                    ice_new[i + 1] = 0;
                }
            }

            // update other parameters due to ice content change
            /********************************************************/
            lidx = 0;
            Lsum = 0.;
            PAST_BOTTOM = false;
            for (i = 0; i <= right + 1; i++) {
                if (i >= left + 1) {
                    if (ice_new[i] != ice[i]) {
                        kappa_new[i] = soil_conductivity(moist[i],
                                                         moist[i] - ice_new[i],
                                                         soil_dens_min[lidx],
                                                         bulk_dens_min[lidx],
                                                         quartz[lidx],
                                                         soil_density[lidx],
                                                         bulk_density[lidx],
                                                         organic[lidx]);
                        Cs_new[i] = volumetric_heat_capacity(
                            bulk_density[lidx] / soil_density[lidx],
                            moist[i] - ice_new[i], ice_new[i], organic[lidx]);
                    }
                }
                if (Zsum[i] > Lsum + depth[lidx] && !PAST_BOTTOM) {
                    Lsum += depth[lidx];
                    lidx++;
                    if (lidx == Nlayers) {
                        PAST_BOTTOM = true;
                        lidx = Nlayers - 1;
                    }
                }
            }
            /*********************************************************/

            // update other states due to ice content change
            for (i = left; i <= right; i++) {
                if (i == 0) {
                    DT[i] = T_2[i + 1] - Ts;
                    DT_up[i] = T_2[i] - Ts;
                    DT_down[i] = T_2[i + 1] - T_2[i];
                }
                else if (i == n - 1) {
                    DT[i] = Tb - T_2[i - 1];
                    DT_up[i] = T_2[i] - T_2[i - 1];
                    DT_down[i] = Tb - T_2[i];
                }
                else {
                    DT[i] = T_2[i + 1] - T_2[i - 1];
                    DT_up[i] = T_2[i] - T_2[i - 1];
                    DT_down[i] = T_2[i + 1] - T_2[i];
                }
                // update Dkappa due to ice content change
                /*******************************************/
                if (i < n - 1) {
                    Dkappa[i] = kappa_new[i + 2] - kappa_new[i];
                }
                else if (!NOFLUX) {
                    Dkappa[i] = kappa_new[i + 2] - kappa_new[i];
                }
                else {
                    Dkappa[i] = kappa_new[i + 1] - kappa_new[i];
                }
                /********************************************/
            }

            for (i = left; i <= right; i++) {
                storage_term =
                    Cs_new[i +
                           1] *
                    (T_2[i] -
                     T0[i +
                        1]) / deltat + T_2[i] *
                    (Cs_new[i + 1] - Cs[i + 1]) / deltat;
                if (!EXP_TRANS) {
                    flux_term1 = Dkappa[i] / alpha[i] * DT[i] / alpha[i];
                    flux_term2 =
                        kappa_new[i +
                                  1] *
                        (DT_down[i] / gamma[i] - DT_up[i] /
                         beta[i]) / (0.5 * alpha[i]);
                }
                else { // grid transformation
                    flux_term1 = Dkappa[i] / 2. * DT[i] / 2. /
                                 (Bexp *
                                  (Zsum[i +
                                        1] + 1.)) / (Bexp * (Zsum[i + 1] + 1.));
                    flux_term2 =
                        kappa_new[i +
                                  1] *
                        ((DT_down[i] -
                          DT_up[i]) /
                         (Bexp *
                          (Zsum[i +
                                1] +
                           1.)) /
                         (Bexp *
                          (Zsum[i +
                                1] +
                           1.)) - DT[i] / 2. /
                         (Bexp * (Zsum[i + 1] + 1.) * (Zsum[i + 1] + 1.)));
                }
                // inelegant fix for "cold nose" problem - when a very cold node skates off to
                // much colder and breaks the second law of thermodynamics (because
                // flux_term1 exceeds flux_term2 in absolute magnitude) - therefore, don't let
                // that node get any colder.  This only seems to happen in the first and
                // second near-surface nodes.
                flux_term = flux_term1 + flux_term2;
                phase_term = CONST_RHOICE * CONST_LATICE *
                             (ice_new[i + 1] - ice[i + 1]) / deltat;
                res[i] = flux_term + phase_term - storage_term;
            }
        } // end of calculation of focus node only
    } // end of non-init
}
