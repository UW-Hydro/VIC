#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <vicNl.h>

#define MAXIT 1000

static char vcid[] = "$Id$";

void setup_frozen_soil(soil_con_struct    soil_con,
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

**********************************************************************/
  extern option_struct options;
  extern debug_struct  debug;

  char    ErrorString[MAXSTRING];
  int     i;

  for(i=0;i<options.Nlayer;i++)
    layer[i] = find_average_layer(layer_wet[i],layer_dry[i],
				   soil_con.depth[i],mu);

  soil_thermal_calc(soil_con, layer, energy, kappa, Cs, moist, 
		    options.Nlayer, Nnodes);

}

void finish_frozen_soil_calcs(energy_bal_struct *energy,
			      layer_data_struct *layer_wet,
			      layer_data_struct *layer_dry,
			      layer_data_struct *layer,
			      soil_con_struct    soil_con,
			      int                Nnodes,
			      int                veg,
			      double             mu,
			      double            *T,
			      double            *dz,
			      double            *kappa,
			      double            *Cs,
			      double            *moist,
			      double            *expt,
			      double            *max_moist) {
/******************************************************************
  finish_frozen_soil_calcs      Keith Cherkauer      July 27, 1998

  This subroutine redistributes soil properties based on the 
  thermal solutions found for the current time step.

******************************************************************/

  extern option_struct options;
  extern debug_struct  debug;

  char    ErrorString[MAXSTRING];
  int     i, j;
  int     index;
  double  fdepth;
  double  tdepth;
  double  unfrozen;
  double  frozen;
  double *old_fdepth;
  double *old_tdepth;
  double *tmp_ptr;

  fdepth = energy->fdepth[0];
  tdepth = energy->fdepth[1];

  /** Calculate New Layer Depths **/
  old_fdepth = (double *)calloc(options.Nlayer,sizeof(double));
  old_tdepth = (double *)calloc(options.Nlayer,sizeof(double));
  for(i=0;i<options.Nlayer;i++) {
    old_fdepth[i] = layer[i].fdepth;
    old_tdepth[i] = layer[i].tdepth;
  }
  find_0_degree_fronts(energy, layer, soil_con.dp, soil_con.depth,
      T, options.Nlayer, Nnodes);

  /***** Find New Layer Temperatures *****/
  find_sublayer_temperatures(layer,T,dz,soil_con.depth,
      energy->fdepth[0],energy->fdepth[1],options.Nlayer,Nnodes);

  /** Store Layer Temperature Values **/
  for(i=0;i<Nnodes;i++) energy->T[i] = T[i];

  if(energy->fdepth[0]>0) energy->frozen = TRUE;
  else energy->frozen = FALSE;

  /** Redistribute Soil Properties for New Frozen Soil Layer Size **/

  if(energy->frozen || fdepth>0.0) {
    for(i=0;i<options.Nlayer;i++) {
      layer_wet[i].fdepth = layer[i].fdepth;
      layer_wet[i].tdepth = layer[i].tdepth;
      layer_wet[i].T      = layer[i].T;
      layer_wet[i].T      = layer[i].T;
      layer_wet[i].T_froz = layer[i].T_froz;
      layer_wet[i].T_froz = layer[i].T_froz;
      layer_wet[i].T_thaw = layer[i].T_thaw;
      layer_wet[i].T_thaw = layer[i].T_thaw;
    }
    redistribute_moisture(layer_wet, energy->fdepth,
			  soil_con.max_moist, old_fdepth, old_tdepth, 
			  soil_con.depth, options.Nlayer);
    if(options.DIST_PRCP) {
      for(i=0;i<options.Nlayer;i++) {
	layer_dry[i].fdepth = layer[i].fdepth;
	layer_dry[i].tdepth = layer[i].tdepth;
	layer_dry[i].T      = layer[i].T;
	layer_dry[i].T      = layer[i].T;
	layer_dry[i].T_froz = layer[i].T_froz;
	layer_dry[i].T_froz = layer[i].T_froz;
	layer_dry[i].T_thaw = layer[i].T_thaw;
	layer_dry[i].T_thaw = layer[i].T_thaw;
      }
      redistribute_moisture(layer_dry, energy->fdepth,
			    soil_con.max_moist, old_fdepth, old_tdepth, 
			    soil_con.depth, options.Nlayer);
    }
  }
  free((char*)old_fdepth);
  free((char*)old_tdepth);

  if(debug.PRT_BALANCE && debug.DEBUG) {
    printf("After Moisture Redistribution\n");
    write_layer(layer,veg,options.Nlayer,soil_con.depth);
  } 

  /** Compute Amount of Unfrozen Moisture in Frozen Layer **/
  for(index=0;index<options.Nlayer;index++) {
    if(layer[index].fdepth > 0.0) {
      if(layer[index].T_froz<0. && layer[index].T_froz != -999.) {
	unfrozen = maximum_unfrozen_water(layer[index].T_froz,
	           soil_con.max_moist[index], soil_con.bubble,
                   soil_con.expt[index]);
	if(unfrozen>soil_con.max_moist[index] || unfrozen<0.) 
	  unfrozen = soil_con.max_moist[index];

	layer_wet[index].unfrozen = unfrozen;
	frozen = layer_wet[index].moist_froz - unfrozen;
	if(frozen < 0.0) {
	  frozen = 0.0;
	  unfrozen = layer_wet[index].moist_froz;
	}
	layer_wet[index].ice = frozen;
	layer_wet[index].moist_froz = unfrozen;
	if(options.DIST_PRCP) {
	  layer_dry[index].unfrozen = unfrozen;
	  frozen = layer_dry[index].moist_froz - unfrozen;
	  if(frozen < 0.0) {
	    frozen = 0.0;
	    unfrozen = layer_dry[index].moist_froz;
	  }
	  layer_dry[index].ice = frozen;
	  layer_dry[index].moist_froz = unfrozen;
	}
      }
      else if(layer[index].T_froz == 0.) {
	layer_wet[index].unfrozen = soil_con.max_moist[index];
	layer_wet[index].ice = 0.;
	if(options.DIST_PRCP) {
	  layer_dry[index].unfrozen = soil_con.max_moist[index];
	  layer_dry[index].ice = 0.;
	}
      }
      else if(layer[index].T_froz != -999.) {
        sprintf(ErrorString,"ERROR: Frozen Layer Temperature > 0C (%lf)",
            layer[index].T_froz);
	vicerror(ErrorString);
      }
    }
    else {
      layer_wet[index].ice=0.0;
      if(options.DIST_PRCP) layer_dry[index].ice=0.0;
    }
  }

  for(j=0;j<Nnodes;j++) {
    if(T[j]<0.) {
      energy->ice[j] = moist[j] - maximum_unfrozen_water(T[j],max_moist[j],
               soil_con.bubble,expt[j]);
      if(energy->ice[j]<0.) energy->ice[j]=0.;
      if(energy->ice[j]>max_moist[j]) energy->ice[j]=max_moist[j];
    }
    else energy->ice[j]=0.;
  }

  for(i=0;i<options.Nlayer;i++)
    layer[i] = find_average_layer(layer_wet[i],layer_dry[i],
				  soil_con.depth[i],mu);

  if(debug.PRT_BALANCE && debug.DEBUG) {
    printf("After Refreezing Moisture\n");
    write_layer(layer,veg,options.Nlayer,soil_con.depth);
  } 

  tmp_ptr = NULL;

  soil_thermal_calc(soil_con, layer, *energy, tmp_ptr,
      tmp_ptr, tmp_ptr, options.Nlayer, Nnodes);

  /** Save Thermal Conductivities for Energy Balance **/
  for(i=0;i<options.Nlayer;i++) {
    layer_wet[i].kappa = layer[i].kappa;
    layer_wet[i].Cs    = layer[i].Cs;
    if(options.DIST_PRCP) {
      layer_dry[i].kappa = layer[i].kappa;
      layer_dry[i].Cs    = layer[i].Cs;
    }
  }
  energy->kappa[0] = layer[0].kappa;
  energy->Cs[0] = layer[0].Cs;
  energy->kappa[1] = layer[1].kappa;
  energy->Cs[1] = layer[1].Cs;
  
}

void solve_T_profile(double *T,
		     double *T0,
		     double *dz,
		     double *kappa,
		     double *Cs,
		     double *moist,
		     double  deltat,
		     double *max_moist,
		     double  bubble,
		     double *expt,
		     double *ice,
		     double *alpha,
		     double *beta,
		     double *gamma,
		     int     Nnodes,
		     char   *FIRST_SOLN,
		     char   FIRST_TIME) {
/**********************************************************************
  This subroutine was written to iteratively solve the soil temperature
  profile using a numerical difference equation.  The solution equation
  is second order in space, and first order in time.
**********************************************************************/

  extern option_struct options;
  extern debug_struct  debug;

  static double *A, *B, *C, *D, *E;

  double maxdiff, diff;
  double threshold = 1.e-3;	/* temperature profile iteration threshold */
  double oldT;
  double fprime;
  int j, Done, ItCount;

  if(FIRST_TIME) {
    A = (double *)calloc(Nnodes,sizeof(double));
    B = (double *)calloc(Nnodes,sizeof(double));
    C = (double *)calloc(Nnodes,sizeof(double));
    D = (double *)calloc(Nnodes,sizeof(double));
    E = (double *)calloc(Nnodes,sizeof(double));
  }
  else {

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
    
    for(j=0;j<Nnodes;j++) T[j]=T0[j];
    
    Done = FALSE;
    ItCount = 0;
    
    while(!Done && ItCount<MAXIT) {
      ItCount++;
      maxdiff=threshold;
      for(j=1;j<Nnodes-1;j++) {
	oldT=T[j];
  
	/**	2nd order variable kappa equation **/
	fprime = (T[j+1]-T[j-1])/alpha[j-1];
	
	if(T[j]>=0) {
	  T[j] = (A[j]*(T[j+1]-T[j-1])
		  + B[j]*(T[j+1]+T[j-1]-gamma[j-1]*fprime)
		  + C[j] + D[j]*(0.-ice[j])) / (E[j]);
	}
	else {
	  T[j] = root_brent(T0[j]-SOIL_DT,T0[j]+SOIL_DT,soil_thermal_eqn, 
			    T[j+1], T[j-1], T0[j], moist[j], max_moist[j], 
			    bubble, expt[j], ice[j], gamma[j-1], fprime, 
			    A[j], B[j], C[j], D[j], E[j]);
	  
	  if(T[j] <= -9998) {
	    error_solve_T_profile(T[j], T[j+1], T[j-1], T0[j], moist[j], 
				  max_moist[j], bubble, expt[j], ice[j], 
				  gamma[j-1], fprime, A[j], B[j], C[j], D[j], 
				  E[j]);
	  }
	}
	
	diff=fabs(oldT-T[j]);
	if(diff>maxdiff) maxdiff=diff;
      }
      
      if(options.NOFLUX) { 
	/** Solve for bottom temperature if using no flux lower boundary **/
	oldT=T[Nnodes-1];
	
	fprime = (T[Nnodes-1]-T[Nnodes-2])/alpha[Nnodes-2];
	
	j = Nnodes-1;
	
	if(T[j]>=0) {
	  T[j] = (A[j]*(T[j]-T[j-1]) + B[j]*(T[j] + T[j-1] - gamma[j-1]*fprime)
		  + C[j] + D[j]*(0.-ice[j])) / E[j];
	}
	else {
	  T[Nnodes-1] = root_brent(T0[Nnodes-1]-SOIL_DT,T0[Nnodes-1]
				   +SOIL_DT,soil_thermal_eqn, T[Nnodes-1],
				   T[Nnodes-2], T0[Nnodes-1], 
				   moist[Nnodes-1], max_moist[Nnodes-1], 
				   bubble, expt[Nnodes-1], ice[Nnodes-1], 
				   gamma[Nnodes-2], fprime, 
				   A[j], B[j], C[j], D[j], E[j]);
	  
	  if(T[Nnodes-1] <= -9998) {
	    error_solve_T_profile(T[Nnodes-1], T[Nnodes-1],
				  T[Nnodes-2], T0[Nnodes-1], 
				  moist[Nnodes-1], max_moist[Nnodes-1], bubble, 
				  expt[Nnodes-1], ice[Nnodes-1], 
				  gamma[Nnodes-2], fprime, 
				  A[j], B[j], C[j], D[j], E[j]);
	  }
	}
	
	diff=fabs(oldT-T[Nnodes-1]);
	if(diff>maxdiff) maxdiff=diff;
      }
    
      if(maxdiff<=threshold) Done=TRUE;
  
    }
  
    if(!Done) {
      fprintf(stderr,"ERROR: Temperature Profile Unable to Converge!!!\n");
      fprintf(stderr,"Dumping Profile Temperatures (last, new).\n");
      for(j=0;j<Nnodes;j++) fprintf(stderr,"%lf\t%lf\n",T0[j],T[j]);
      vicerror("ERROR: Cannot solve temperature profile:\n\tToo Many Iterations in solve_T_profile");
    }

  }
  
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
  
  fprintf(stderr,"TL\t%lf\n",TL);
  fprintf(stderr,"TU\t%lf\n",TU);
  fprintf(stderr,"T0\t%lf\n",T0);
  fprintf(stderr,"moist\t%lf\n",moist);
  fprintf(stderr,"max_moist\t%lf\n",max_moist);
  fprintf(stderr,"bubble\t%lf\n",bubble);
  fprintf(stderr,"expt\t%lf\n",expt);
  fprintf(stderr,"ice0\t%lf\n",ice0);
  fprintf(stderr,"gamma\t%lf\n",gamma);
  fprintf(stderr,"fprime\t%lf\n",fprime);
  fprintf(stderr,"A\t%lf\n",A);
  fprintf(stderr,"B\t%lf\n",B);
  fprintf(stderr,"C\t%lf\n",C);
  fprintf(stderr,"D\t%lf\n",D);
  fprintf(stderr,"E\t%lf\n",E);

  vicerror("Finished dumping values for solve_T_profile.\nTry increasing SOIL_DT to get model to complete cell.\nThen check output for instabilities.");

}

#undef MAXIT
