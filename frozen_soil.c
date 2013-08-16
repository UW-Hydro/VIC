#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <vicNl.h>

#define MAXIT 1000

static char vcid[] = "$Id: frozen_soil.c,v 4.2.2.3 2004/09/22 00:40:33 vicadmin Exp $";

void setup_frozen_soil(soil_con_struct   *soil_con,
		       layer_data_struct *layer_wet,
		       layer_data_struct *layer_dry,
		       layer_data_struct *layer,
		       energy_bal_struct  energy,
		       int                rec,
		       int                veg,
		       int                Nnodes,
		       double             mu,
		       double            *kappa,
		       double            *Cs,
		       double            *moist) {
/**********************************************************************
  setup_frozen_soil	Keith Cherkauer		July 27, 1998

  This subroutine prepares the data arrays needed to solve the soil
  thermal fluxes through all soil thermal nodes.

  soil_con_struct    soil_con   soil parameter structure
  layer_data_struct *layer_wet  soil variables for wet fraction
  layer_data_struct *layer_dry  soil variables for dry fraction
  layer_data_struct *layer      average soil variables
  energy_bal_struct  energy     energy balance structure
  int                rec        record number
  int                veg        vegetation type number
  int                Nnodes     number of soil thermal nodes 
  double             mu         fraction of grid cell that received precip
  double            *kappa      soil layer thermal conductivity (W/m/K)
  double            *Cs         soil layer heat capacity (J/m^3/K)
  double            *moist      soil layer moisture (mm)
  
  Modifications
  04-Jun-04 Added descriptive error message to beginning of screen dump
	    in error_print_solve_T_profile.			TJB
  21-Sep-04 Added ErrorString to store error messages from
	    root_brent.						TJB

**********************************************************************/
  extern option_struct options;
#if LINK_DEBUG
  extern debug_struct  debug;
#endif

  int nidx;

  for(nidx=0;nidx<Nnodes;nidx++) {
    moist[nidx] = energy.moist[nidx];
    kappa[nidx] = energy.kappa_node[nidx];
    Cs[nidx]    = energy.Cs_node[nidx];
  }

}

void finish_frozen_soil_calcs(energy_bal_struct *energy,
			      layer_data_struct *layer_wet,
			      layer_data_struct *layer_dry,
			      layer_data_struct *layer,
			      soil_con_struct   *soil_con,
			      int                Nnodes,
			      int                veg,
			      double             mu,
			      double            *T,
			      double            *kappa,
			      double            *Cs,
			      double            *moist) {
/******************************************************************
  finish_frozen_soil_calcs      Keith Cherkauer      July 27, 1998

  This subroutine redistributes soil properties based on the 
  thermal solutions found for the current time step.

  Modifications:
  3-12-03 Modified so that soil layer ice content is only 
          calculated if the frozen soil algorithm is implemented 
          and active in the current grid cell.               KAC

******************************************************************/

  extern option_struct options;
#if LINK_DEBUG
  extern debug_struct  debug;
#endif

  char    ErrorString[MAXSTRING];
  int     i, j;
  int     index;
  double  fdepth;
  double  tdepth;
  double  unfrozen;
  double  frozen;
  double  old_fdepth[MAX_LAYERS];
  double  old_tdepth[MAX_LAYERS];
  double *tmp_ptr;

  find_0_degree_fronts(energy, soil_con->dz_node, T, Nnodes);

  /** Store Layer Temperature Values **/
  for(i=0;i<Nnodes;i++) energy->T[i] = T[i];

  if(energy->Nfrost>0) energy->frozen = TRUE;
  else energy->frozen = FALSE;

  /** Redistribute Soil Properties for New Frozen Soil Layer Size **/
  if(soil_con->FS_ACTIVE && options.FROZEN_SOIL)
    estimate_layer_ice_content(layer_wet, soil_con->dz_node, energy->T,
			       soil_con->max_moist_node, 
#if QUICK_FS
			       soil_con->ufwc_table_node,
#else
			       soil_con->expt_node, soil_con->bubble_node, 
#endif
			       soil_con->depth, soil_con->max_moist, 
#if QUICK_FS
			       soil_con->ufwc_table_layer,
#else
			       soil_con->expt, soil_con->bubble, 
#endif
			       soil_con->bulk_density,
			       soil_con->soil_density, soil_con->quartz,
			       soil_con->layer_node_fract, Nnodes, 
			       options.Nlayer, soil_con->FS_ACTIVE);
  if(options.DIST_PRCP && soil_con->FS_ACTIVE && options.FROZEN_SOIL)
    estimate_layer_ice_content(layer_dry, soil_con->dz_node, energy->T,
			       soil_con->max_moist_node, 
#if QUICK_FS
			       soil_con->ufwc_table_node,
#else
			       soil_con->expt_node, soil_con->bubble_node, 
#endif
			       soil_con->depth, soil_con->max_moist, 
#if QUICK_FS
			       soil_con->ufwc_table_layer,
#else
			       soil_con->expt, soil_con->bubble, 
#endif
			       soil_con->bulk_density, soil_con->soil_density, 
			       soil_con->quartz, soil_con->layer_node_fract, 
			       Nnodes, options.Nlayer, soil_con->FS_ACTIVE);

#if LINK_DEBUG
  if(debug.PRT_BALANCE && debug.DEBUG) {
    printf("After Moisture Redistribution\n");
    write_layer(layer,veg,options.Nlayer,soil_con->depth);
  } 
#endif
  
}

void solve_T_profile(double *T,
		     double *T0,
		     double *dz,
		     double *kappa,
		     double *Cs,
		     double *moist,
		     double  deltat,
		     double *max_moist,
		     double *bubble,
		     double *expt,
		     double *ice,
		     double *alpha,
		     double *beta,
		     double *gamma,
#if QUICK_FS
		     double ***ufwc_table_node,
#endif
		     int     Nnodes,
		     char   *FIRST_SOLN,
		     char    FIRST_TIME, 
		     int     FS_ACTIVE) {
/**********************************************************************
  This subroutine was written to iteratively solve the soil temperature
  profile using a numerical difference equation.  The solution equation
  is second order in space, and first order in time.
**********************************************************************/

  extern option_struct options;
#if LINK_DEBUG
  extern debug_struct  debug;
#endif
  
  static double A[MAX_NODES];
  static double B[MAX_NODES];
  static double C[MAX_NODES];
  static double D[MAX_NODES];
  static double E[MAX_NODES];

  double *aa, *bb, *cc, *dd, *ee;

  char   Error;
  int    try;
  double maxdiff, diff;
  double oldT;
  double fprime;
  int    j, Done, ItCount;

  if(FIRST_SOLN[0]) {
    FIRST_SOLN[0] = FALSE;
    for(j=1;j<Nnodes-1;j++) {
      A[j] = beta[j-1]*deltat*(kappa[j+1]-kappa[j-1]);
      B[j] = 2.*alpha[j-1]*alpha[j-1]*deltat*kappa[j];
      C[j] = alpha[j-1]*alpha[j-1]*beta[j-1]*Cs[j]*T0[j];
      D[j] = alpha[j-1]*alpha[j-1]*beta[j-1]*ice_density*Lf;
      E[j] = alpha[j-1]*alpha[j-1]*beta[j-1]*Cs[j] 
	+ 4.*kappa[j]*alpha[j-1]*alpha[j-1]*deltat;
    }
    if(options.NOFLUX) {
      j = Nnodes-1;
      A[j] = beta[j-1]*deltat*(kappa[j]-kappa[j-1]);
      B[j] = 2.*alpha[j-1]*alpha[j-1]*deltat*kappa[j];
      C[j] = alpha[j-1]*alpha[j-1]*beta[j-1]*Cs[j]*T0[j];
      D[j] = alpha[j-1]*alpha[j-1]*beta[j-1]*ice_density*Lf;
      E[j] = alpha[j-1]*alpha[j-1]*beta[j-1]*Cs[j] 
	+ 4.*kappa[j]*alpha[j-1]*alpha[j-1]*deltat;
    }
  }
    
  aa = &A[0];
  bb = &B[0];
  cc = &C[0];
  dd = &D[0];
  ee = &E[0];

  for(j=0;j<Nnodes;j++) T[j]=T0[j];

#if QUICK_FS
  Error = calc_soil_thermal_fluxes(Nnodes, T, T0, moist, max_moist, ice, 
				   bubble, expt, alpha, gamma, aa, bb, cc, 
				   dd, ee, ufwc_table_node, FS_ACTIVE);
#else
  Error = calc_soil_thermal_fluxes(Nnodes, T, T0, moist, max_moist, ice, 
				   bubble, expt, alpha, gamma, aa, bb, cc, 
				   dd, ee, FS_ACTIVE);
#endif
  
}
  

int calc_soil_thermal_fluxes(int     Nnodes,
			     double *T,
			     double *T0,
			     double *moist,
			     double *max_moist,
			     double *ice,
			     double *bubble,
			     double *expt,
			     double *alpha,
			     double *gamma,
			     double *A, 
			     double *B, 
			     double *C, 
			     double *D, 
			     double *E,
#if QUICK_FS
			     double ***ufwc_table_node,
#endif
			     char    FS_ACTIVE) {
  
  /** Eventually the nodal ice contents will also have to be updated **/

  extern option_struct options;

  int    Error;
  char   Done;
  int    j;
  int    ItCount;
  double threshold = 1.e-2;	/* temperature profile iteration threshold */
  double maxdiff;
  double diff;
  double oldT;
  double fprime;
  char ErrorString[MAXSTRING];

  Error = 0;
  Done = FALSE;
  ItCount = 0;
  
  while(!Done && Error==0 && ItCount<MAXIT) {
    ItCount++;
    maxdiff=threshold;
    for(j=1;j<Nnodes-1;j++) {
      oldT=T[j];
      
      /**	2nd order variable kappa equation **/
      fprime = (T[j+1]-T[j-1])/alpha[j-1];
      
      if(T[j] >= 0 || !FS_ACTIVE || !options.FROZEN_SOIL) {
	T[j] = (A[j]*(T[j+1]-T[j-1])
		+ B[j]*(T[j+1]+T[j-1]-gamma[j-1]*fprime)
		+ C[j] + D[j]*(0.-ice[j])) / (E[j]);
      }
      else {
#if QUICK_FS
	T[j] = root_brent(T0[j]-(SOIL_DT), T0[j]+(SOIL_DT),
			  ErrorString, soil_thermal_eqn, 
			  T[j+1], T[j-1], T0[j], moist[j], max_moist[j], 
			  ufwc_table_node[j], ice[j], gamma[j-1], fprime, 
			  A[j], B[j], C[j], D[j], E[j]);
#else
	T[j] = root_brent(T0[j]-(SOIL_DT), T0[j]+(SOIL_DT),
			  ErrorString, soil_thermal_eqn, 
			  T[j+1], T[j-1], T0[j], moist[j], max_moist[j], 
			  bubble[j], expt[j], ice[j], gamma[j-1], fprime, 
			  A[j], B[j], C[j], D[j], E[j]);
#endif
	
	if(T[j] <= -9998) 
	  error_solve_T_profile(T[j], T[j+1], T[j-1], T0[j], moist[j], 
				max_moist[j], bubble[j], expt[j], ice[j], 
				gamma[j-1], fprime, A[j], B[j], C[j], D[j], 
				E[j], ErrorString);
      }
      
      diff=fabs(oldT-T[j]);
      if(diff > maxdiff) maxdiff=diff;
    }
    
    if(options.NOFLUX) { 
      /** Solve for bottom temperature if using no flux lower boundary **/
      oldT=T[Nnodes-1];
      
      fprime = (T[Nnodes-1]-T[Nnodes-2])/alpha[Nnodes-2];
      
      j = Nnodes-1;
      
      if(T[j] >= 0 || !FS_ACTIVE || !options.FROZEN_SOIL) {
	T[j] = (A[j]*(T[j]-T[j-1]) + B[j]*(T[j] + T[j-1] - gamma[j-1]*fprime)
		+ C[j] + D[j]*(0.-ice[j])) / E[j];
      }
      else {
#if QUICK_FS
	T[Nnodes-1] = root_brent(T0[Nnodes-1]-SOIL_DT,T0[Nnodes-1]+SOIL_DT,
				 ErrorString, soil_thermal_eqn, T[Nnodes-1],
				 T[Nnodes-2], T0[Nnodes-1], 
				 moist[Nnodes-1], max_moist[Nnodes-1], 
				 ufwc_table_node[Nnodes-1], 
				 ice[Nnodes-1], 
				 gamma[Nnodes-2], fprime, 
				 A[j], B[j], C[j], D[j], E[j]);
#else
	T[Nnodes-1] = root_brent(T0[Nnodes-1]-SOIL_DT,T0[Nnodes-1]+SOIL_DT,
				 ErrorString, soil_thermal_eqn, T[Nnodes-1],
				 T[Nnodes-2], T0[Nnodes-1], 
				 moist[Nnodes-1], max_moist[Nnodes-1], 
				 bubble[j], expt[Nnodes-1], ice[Nnodes-1], 
				 gamma[Nnodes-2], fprime, 
				 A[j], B[j], C[j], D[j], E[j]);
#endif
	
	if(T[j] <= -9998) 
	  error_solve_T_profile(T[Nnodes-1], T[Nnodes-1],
				T[Nnodes-2], T0[Nnodes-1], 
				moist[Nnodes-1], max_moist[Nnodes-1], 
				bubble[Nnodes-1], 
				expt[Nnodes-1], ice[Nnodes-1], 
				gamma[Nnodes-2], fprime, 
				A[j], B[j], C[j], D[j], E[j], ErrorString);
      }
      
      diff=fabs(oldT-T[Nnodes-1]);
      if(diff>maxdiff) maxdiff=diff;
    }
    
    if(maxdiff <= threshold) Done=TRUE;
    
  }
  
  if(!Done && !Error) {
    fprintf(stderr,"ERROR: Temperature Profile Unable to Converge!!!\n");
    fprintf(stderr,"Dumping Profile Temperatures (last, new).\n");
    for(j=0;j<Nnodes;j++) fprintf(stderr,"%f\t%f\n",T0[j],T[j]);
    vicerror("ERROR: Cannot solve temperature profile:\n\tToo Many Iterations in solve_T_profile");
  }

  return (Error);

}

double error_solve_T_profile (double Tj, ...) {

  va_list ap;

  double error;

  va_start(ap,Tj);
  error = error_print_solve_T_profile(Tj, ap);
  va_end(ap);

  return error;

}

double error_print_solve_T_profile(double T, va_list ap) {

  double TL;
  double TU;
  double T0;
  double moist;
  double max_moist;
  double bubble;
  double expt;
  double ice0;
  double gamma;
  double fprime;
  double A;
  double B;
  double C;
  double D;
  double E;
  char *ErrorString;

  TL        = (double) va_arg(ap, double);
  TU        = (double) va_arg(ap, double);
  T0        = (double) va_arg(ap, double);
  moist     = (double) va_arg(ap, double);
  max_moist = (double) va_arg(ap, double);
  bubble    = (double) va_arg(ap, double);
  expt      = (double) va_arg(ap, double);
  ice0      = (double) va_arg(ap, double);
  gamma     = (double) va_arg(ap, double);
  fprime    = (double) va_arg(ap, double);
  A         = (double) va_arg(ap, double);
  B         = (double) va_arg(ap, double);
  C         = (double) va_arg(ap, double);
  D         = (double) va_arg(ap, double);
  E         = (double) va_arg(ap, double);
  ErrorString = (char*) va_arg(ap, char *);
  
  fprintf(stderr, "%s", ErrorString);
  fprintf(stderr, "ERROR: solve_T_profile failed to converge to a solution in root_brent.  Variable values will be dumped to the screen, check for invalid values.\n");

  fprintf(stderr,"TL\t%f\n",TL);
  fprintf(stderr,"TU\t%f\n",TU);
  fprintf(stderr,"T0\t%f\n",T0);
  fprintf(stderr,"moist\t%f\n",moist);
  fprintf(stderr,"max_moist\t%f\n",max_moist);
  fprintf(stderr,"bubble\t%f\n",bubble);
  fprintf(stderr,"expt\t%f\n",expt);
  fprintf(stderr,"ice0\t%f\n",ice0);
  fprintf(stderr,"gamma\t%f\n",gamma);
  fprintf(stderr,"fprime\t%f\n",fprime);
  fprintf(stderr,"A\t%f\n",A);
  fprintf(stderr,"B\t%f\n",B);
  fprintf(stderr,"C\t%f\n",C);
  fprintf(stderr,"D\t%f\n",D);
  fprintf(stderr,"E\t%f\n",E);

  vicerror("Finished dumping values for solve_T_profile.\nTry increasing SOIL_DT to get model to complete cell.\nThen check output for instabilities.\n");
  
  return(0.0);

}

#undef MAXIT
