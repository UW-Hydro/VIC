#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <vicNl.h>
#include <stdarg.h>

#define MAXIT 1000

static char vcid[] = "$Id$";

int calc_layer_average_thermal_props(energy_bal_struct *energy,
				     layer_data_struct *layer,
				     soil_con_struct   *soil_con,
				     int                Nnodes,
				     int                veg,
				     double            *T) {
  /******************************************************************
  calc_layer_average_thermal_props      Keith Cherkauer      July 27, 1998

  This subroutine redistributes soil properties based on the 
  thermal solutions found for the current time step.

  Modifications:
  3-12-03 Modified so that soil layer ice content is only 
          calculated if the frozen soil algorithm is implemented 
          and active in the current grid cell.               KAC
  2007-Apr-24 Removed setup_frozen_soil function (was above this one).	JCA
  2007-Apr-24 Added functionality for EXP_TRANS option.			JCA
              (including passing Zsum_node to find_0_degree_fronts)
  2007-Apr-24 For IMPLICIT option, added Ming Pan's new functions
	      (solve_T_profile_implicit and fda_heat_eqn).		JCA
  2007-Aug-09 Added features for EXCESS_ICE option.			JCA
  2009-Feb-09 Removed dz_node from call to find_0_degree_front.		KAC via TJB
  2009-Feb-09 Modified to handle error flags and pass Zsum_node instead 
              of dz_node to esimate_layer_ice_content.			KAC via TJB
  2009-Mar-16 Added resid_moist to argument list of
	      estimate_layer_ice_content().  This allows computation
	      of min_liq, the minimum allowable liquid water content
	      in each layer as a function of temperature.		TJB
  2009-Jun-10 Fixed incorrect placement of checks on ErrorFlag.		TJB
  2009-Jul-31 Removed unused layer_node_fract array from call to
	      estimate_layer_ice_content().				TJB
  2009-Dec-11 Removed min_liq and options.MIN_LIQ.  As a result, no
	      longer need to include resid_moist in the arg list for
	      estimate_layer_ice_content().				TJB
  2011-May-24 Added calls to estimate_layer_ice_content_quick_flux to
	      handle the cases for QUICK_FLUX = TRUE (or any other case
	      in which the full finite-element T profile is not used).
	      Changed name from finish_frozen_soil_calcs to
	      calc_layer_average_thermal_props.  This function now can
	      be used to compute soil layer average T and ice, regardless
	      of the settings of FROZEN_SOIL, QUICK_FLUX, etc.		TJB
  2012-Jan-16 Removed LINK_DEBUG code					BN
  2013-Dec-26 Removed EXCESS_ICE option.				TJB
  2013-Dec-27 Moved SPATIAL_FROST to options_struct.			TJB
  2013-Dec-27 Removed QUICK_FS option.					TJB
  2014-Mar-28 Removed DIST_PRCP option.					TJB
******************************************************************/

  extern option_struct options;
  int     i, ErrorFlag;

  if (options.FROZEN_SOIL && soil_con->FS_ACTIVE)
    find_0_degree_fronts(energy, soil_con->Zsum_node, T, Nnodes);
  else
    energy->Nfrost = 0;

  /** Store Layer Temperature Values **/
  for(i=0;i<Nnodes;i++) energy->T[i] = T[i];

  if(energy->Nfrost>0) energy->frozen = TRUE;
  else energy->frozen = FALSE;

  /** Compute Soil Layer average  properties **/
  if (options.QUICK_FLUX) {
    ErrorFlag = estimate_layer_ice_content_quick_flux(layer, soil_con->depth, soil_con->dp,
					   energy->T[0], energy->T[1], soil_con->avg_temp,
					   soil_con->max_moist, 
					   soil_con->expt, soil_con->bubble, 
					   soil_con->frost_fract, soil_con->frost_slope, soil_con->FS_ACTIVE);
    if ( ErrorFlag == ERROR ) return (ERROR);
  }
  else {
    ErrorFlag = estimate_layer_ice_content(layer, soil_con->Zsum_node, energy->T,
					   soil_con->max_moist_node, 
					   soil_con->expt_node, soil_con->bubble_node, 
					   soil_con->depth, soil_con->max_moist, 
					   soil_con->expt, soil_con->bubble, 
					   soil_con->frost_fract, soil_con->frost_slope, 
					   Nnodes, options.Nlayer, soil_con->FS_ACTIVE);
    if ( ErrorFlag == ERROR ) return (ERROR);
  }
  
  return (0);
  
}

int  solve_T_profile(double *T,
		     double *T0,
		     char   *Tfbflag,
		     int    *Tfbcount,
		     double *Zsum,
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
		     double Dp,
		     double *depth,
		     int     Nnodes,
		     int    *FIRST_SOLN,
		     int     FS_ACTIVE,
		     int     NOFLUX,
		     int     EXP_TRANS, 
		     int veg_class) {
/**********************************************************************
  This subroutine was written to iteratively solve the soil temperature
  profile using a numerical difference equation.  The solution equation
  is second order in space, and first order in time.

  Modifications:
  2007-Apr-11 Changed type of Error from char to int.				GCT
  2007-Apr-24 Added option for EXP_TRANS.					JCA
              (including passing in dz_node, Zsum, Dp, depth, EXP_TRANS,
              veg_class; and removing FIRST_TIME and fprime)
  2007-Apr-24 Rearranged terms in finite-difference heat equation (equation 8
              of Cherkauer et al. (1999)); Therefore the constants (A-E)
              are calculated in a new way.  These constants are equal to the
              constants in each of the terms in equation 8 multiplied by
              alpha^2*deltat.  This was done to make EXP_TRANS option
              easier to code.							JCA
  2007-Apr-24 Replaced second term of heat flux with alternate derivative
              approximation (a form found in most text books).			JCA
  2007-Aug-08 Added option for EXCESS_ICE.					JCA
  2007-Oct-08 Fixed error in EXP_TRANS formulation.				JCA
  2007-Oct-11 Fixed error in EXP_TRANS formulation.				JCA
  2009-Feb-09 Removed dz_node from call to solve_T_profile and 
              solve_T_profile_implicit.                                         KAC
  2009-Jun-19 Added T fbflag to indicate whether TFALLBACK occurred.		TJB
  2009-Sep-19 Added T fbcount to count TFALLBACK occurrences.			TJB
  2012-Jan-16 Removed LINK_DEBUG code						BN
  2013-Dec-26 Removed EXCESS_ICE option.				TJB
  2013-Dec-27 Moved SPATIAL_FROST to options_struct.			TJB
  2013-Dec-27 Removed QUICK_FS option.					TJB
**********************************************************************/

  extern option_struct options;
  static double A[MAX_NODES];
  static double B[MAX_NODES];
  static double C[MAX_NODES];
  static double D[MAX_NODES];
  static double E[MAX_NODES];

  double *aa, *bb, *cc, *dd, *ee, Bexp;

  int    Error;
  int    j;

  if(FIRST_SOLN[0]) {
    //fprintf(stderr,"*************EXPLICIT SOLUTION***********\n");
    
    if(EXP_TRANS)
      Bexp = logf(Dp+1.)/(double)(Nnodes-1); 

    FIRST_SOLN[0] = FALSE;
    if(!EXP_TRANS) {
      for(j=1;j<Nnodes-1;j++) {
	A[j] = Cs[j]*alpha[j-1]*alpha[j-1];
	B[j] = (kappa[j+1]-kappa[j-1])*deltat;

	//C[j] = 2*deltat*kappa[j]*powf(alpha[j-1],2.)/(powf(gamma[j-1],2.)+powf(beta[j-1],2.)); // old formulation
	//D[j] = 2*deltat*kappa[j]*(gamma[j-1]-beta[j-1])/(powf(gamma[j-1],2.)+powf(beta[j-1],2.));  // old formulation

	C[j] = 2*deltat*kappa[j]*alpha[j-1]/gamma[j-1]; // new formulation
	D[j] = 2*deltat*kappa[j]*alpha[j-1]/beta[j-1];  // new formulation

	E[j] = ice_density*Lf*alpha[j-1]*alpha[j-1];
      }
      if(NOFLUX) {
	j = Nnodes-1;
	A[j] = Cs[j]*alpha[j-1]*alpha[j-1];
	B[j] = (kappa[j]-kappa[j-1])*deltat;


	//C[j] = 2*deltat*kappa[j]*powf(alpha[j-1],2.)/(powf(gamma[j-1],2.)+powf(beta[j-1],2.)); // old formulation
	//D[j] = 2*deltat*kappa[j]*(gamma[j-1]-beta[j-1])/(powf(gamma[j-1],2.)+powf(beta[j-1],2.));  // old formulation

	C[j] = 2*deltat*kappa[j]*alpha[j-1]/gamma[j-1]; //new formulation
	D[j] = 2*deltat*kappa[j]*alpha[j-1]/beta[j-1]; //new formulation

	E[j] = ice_density*Lf*alpha[j-1]*alpha[j-1];
      }
    }
    else { //grid transformation terms
      for(j=1;j<Nnodes-1;j++) {
	A[j] = 4*Bexp*Bexp*Cs[j]*(Zsum[j]+1)*(Zsum[j]+1);
	B[j] = (kappa[j+1]-kappa[j-1])*deltat;
	C[j] = 4*deltat*kappa[j];
	D[j] = 2*deltat*kappa[j]*Bexp;
	E[j] = 4*Bexp*Bexp*ice_density*Lf*(Zsum[j]+1)*(Zsum[j]+1);
      }
      if(NOFLUX) {
	j = Nnodes-1;
	A[j] = 4*Bexp*Bexp*Cs[j]*(Zsum[j]+1)*(Zsum[j]+1);
	B[j] = (kappa[j]-kappa[j-1])*deltat;
	C[j] = 4*deltat*kappa[j];
	D[j] = 2*deltat*kappa[j]*Bexp;
	E[j] = 4*Bexp*Bexp*ice_density*Lf*(Zsum[j]+1)*(Zsum[j]+1);
      }
    }
  }
  
  aa = &A[0];
  bb = &B[0];
  cc = &C[0];
  dd = &D[0];
  ee = &E[0];
  
  for(j=0;j<Nnodes;j++) T[j]=T0[j];

  Error = calc_soil_thermal_fluxes(Nnodes, T, T0, Tfbflag, Tfbcount, moist, max_moist, ice, 
				   bubble, expt, alpha, gamma, aa, bb, cc, 
				   dd, ee, 
				   FS_ACTIVE, NOFLUX, EXP_TRANS, veg_class);

  return ( Error );
  
}
 


int solve_T_profile_implicit(double *T,                           // update
			     double *T0,                    // keep
			     char   *Tfbflag,
			     int    *Tfbcount,
			     double *Zsum,                  // soil parameter
			     double *kappa,                 // update if necessary
			     double *Cs,                    // update if necessary
			     double *moist,                 // keep
			     double  deltat,                // model parameter
			     double *max_moist,             // soil parameter
			     double *bubble,                // soil parameter
			     double *expt,                  // soil parameter
			     double *ice,                   // update if necessary
			     double *alpha,                 // soil parameter
			     double *beta,                  // soil parameter
			     double *gamma,                 // soil parameter
			     double Dp,                     // soil parameter
			     int     Nnodes,               // model parameter
			     int   *FIRST_SOLN,            // update
			     int     FS_ACTIVE,
			     int  NOFLUX,
			     int EXP_TRANS,
			     int veg_class,                // model parameter
			     double *bulk_dens_min,          // soil parameter
			     double *soil_dens_min,          // soil parameter
			     double *quartz,                // soil parameter
			     double *bulk_density,          // soil parameter
			     double *soil_density,          // soil parameter
			     double *organic,                // soil parameter
			     double *depth)                 // soil parameter
{    
  /**********************************************************************
  This subroutine was written to iteratively solve the soil temperature
  profile using a numerical difference equation.  The solution equation
  is second order in space, and first order in time.
  By Ming Pan, mpan@Princeton.EDU

  Modifications:
  2006-Aug-08 Integrated with 4.1.0 (from 4.0.3).			JCA
  2006-Aug-08 Added NOFLUX option.					JCA
  2006-Aug-08 Added EXP_TRANS option - this allows the heat flux equations
	      to be solved on a transformed grid using an exponential
	      distribution.						JCA
  2006-Aug-09 Included terms needed for Cs and kappa nodal updating for
	      new ice.							JCA
  2007-Aug-08 Added EXCESS_ICE OPTION.					JCA
  2009-Feb-09 Removed dz_node from call to solve_T_profile and 
              solve_T_profile_implicit.					KAC
  2011-Jun-03 Added ORGANIC_FRACT option.				TJB
  2011-Jun-03 Added options.ORGANIC_FRACT.  Soil properties now take
	      organic fraction into account.				TJB
  2012-Jan-16 Removed LINK_DEBUG code					BN
  2013-Dec-26 Removed EXCESS_ICE option.				TJB
  2014-Jan-14 Modified cold nose hack to also cover warm nose case.
  **********************************************************************/
  
  extern option_struct options;
  int  n, Error;
  double res[MAX_NODES];
  void (*vecfunc)(double *, double *, int, int, ...);
  int j;

  if(FIRST_SOLN[0]) 
    FIRST_SOLN[0] = FALSE;
  
  // initialize fda_heat_eqn:
  // pass model parameters, initial states, and soil parameters
  // it MUST be initialized before Newton-Raphson searching  
  if(!NOFLUX)
    n = Nnodes-2;
  else
    n = Nnodes-1;
  
  fda_heat_eqn(&T[1], res, n, 1, deltat, FS_ACTIVE, NOFLUX, EXP_TRANS, T0, moist, ice, kappa, Cs, max_moist, bubble, expt, 
	       alpha, beta, gamma, Zsum, Dp, bulk_dens_min, soil_dens_min, quartz, bulk_density, soil_density, organic, depth, options.Nlayer);
  
  // modified Newton-Raphson to solve for new T
  vecfunc = &(fda_heat_eqn);
  Error = newt_raph(vecfunc, &T[1], n);
 
  // update temperature boundaries
  if(Error == 0 ){
    T[0] = T0[0]; //surface
    if(!NOFLUX)
      T[Nnodes-1] = T0[Nnodes-1]; //bottom boundary
    if (options.TFALLBACK) {
      // HACK to prevent runaway cold nose
      // Handle the case in which the a node was colder than both the nodes above and below
      // in the last time step, and that both differences have increased between the last
      // time step and the current one.
      for(j=1;j<Nnodes-1;j++) {
        if ( ( T0[j-1]-T0[j] > 0 && T0[j+1]-T0[j] > 0
               && (T[j-1]-T[j]) - (T0[j-1]-T0[j]) > 0
               && (T[j+1]-T[j]) - (T0[j+1]-T0[j]) > 0 )
          || ( T0[j-1]-T0[j] < 0 && T0[j+1]-T0[j] < 0
               && (T[j-1]-T[j]) - (T0[j-1]-T0[j]) < 0
               && (T[j+1]-T[j]) - (T0[j+1]-T0[j]) < 0 ) ) {
          T[j] = 0.5*(T[j-1]+T[j+1]); // crude fix for now; just average the T's without taking distance, conductivities into account
          Tfbflag[j] = 1;
          Tfbcount[j]++;
        }
      }
    }
  }

  return (Error);

}

 

int calc_soil_thermal_fluxes(int     Nnodes,
			     double *T,
			     double *T0,
			     char   *Tfbflag,
			     int    *Tfbcount,
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
			     int    FS_ACTIVE, 
			     int    NOFLUX,
			     int EXP_TRANS,
			     int veg_class) {
  
  /**********************************************************************
  Modifications:
  2007-Apr-24 Added EXP_TRANS option.						JCA
	      (including passing in EXP_TRANS and veg_class; and removing fprime)
  2007-Apr-24 Rearranged terms in finite-difference heat equation (equation 8
	      of Cherkauer et al. (1999)).  see note in solve_T_profile.
	      This affects the equation for T[j].				JCA
  2007-Apr-24 Passed j to soil_thermal_eqn for "cold nose" problem in
	      explicit solution.						JCA
  2007-Aug-08 Added EXCESS_ICE option.						JCA
  2007-Aug-31 Checked root_brent return value against -998 rather than -9998.	JCA
  2009-May-22 Added TFALLBACK value to options.CONTINUEONERROR.  This
	      allows simulation to continue when energy balance fails
	      to converge by using previous T value.				TJB
  2009-Jun-19 Added T fbflag to indicate whether TFALLBACK occurred.		TJB
  2009-Sep-19 Added T fbcount to count TFALLBACK occurrences.			TJB
  2009-Nov-11 Changed the value of T for TFALLBACK from oldT to T0.		TJB
  2010-Feb-03 Corrected typo in initialization of Tfbflag.			TJB
  2010-Mar-08 Added TFallback logic for case in which max iterations exceeded.	TJB
  2010-Apr-24 Added initialization of Tfbcount.					TJB
  2010-Apr-24 Added hack to prevent cold nose.  Only active when TFALLBACK
	      is TRUE.								TJB
  2013-Dec-26 Removed EXCESS_ICE option.				TJB
  2013-Dec-27 Moved SPATIAL_FROST to options_struct.			TJB
  2013-Dec-27 Removed QUICK_FS option.					TJB
  2014-Mar-28 Modified cold nose hack to also cover warm nose case.	TJB
  **********************************************************************/

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
  char ErrorString[MAXSTRING];
  double Tlast[MAX_NODES];

  Error = 0;
  Done = FALSE;
  ItCount = 0;
 
  /* initialize Tlast */
  for(j=0;j<Nnodes;j++)
    Tlast[j] = T[j];

  /* initialize Tfbflag, Tfbcount */
  for(j=0;j<Nnodes;j++) {
    Tfbflag[j] = 0;
    Tfbcount[j] = 0;
  }

  while(!Done && Error==0 && ItCount<MAXIT) {
    ItCount++;
    maxdiff=threshold;
    for(j=1;j<Nnodes-1;j++) {
      oldT=T[j];
      
      /**	2nd order variable kappa equation **/
      
      if(T[j] >= 0 || !FS_ACTIVE || !options.FROZEN_SOIL) {
	if(!EXP_TRANS)
	  T[j] = (A[j]*T0[j]+B[j]*(T[j+1]-T[j-1])+C[j]*T[j+1]+D[j]*T[j-1]+E[j]*(0.-ice[j]))/(A[j]+C[j]+D[j]);
	else
	  T[j] = (A[j]*T0[j]+B[j]*(T[j+1]-T[j-1])+C[j]*(T[j+1]+T[j-1])-D[j]*(T[j+1]-T[j-1])+E[j]*(0.-ice[j]))/(A[j]+2.*C[j]);
      }
      else {
	T[j] = root_brent(T0[j]-(SOIL_DT), T0[j]+(SOIL_DT),
			  ErrorString, soil_thermal_eqn, 
			  T[j+1], T[j-1], T0[j], moist[j], max_moist[j], 
			  bubble[j], expt[j], 
			  ice[j], gamma[j-1], 
			  A[j], B[j], C[j], D[j], E[j], EXP_TRANS, j);
	if(T[j] <= -998 ) {
          if (options.TFALLBACK) {
            T[j] = T0[j];
            Tfbflag[j] = 1;
            Tfbcount[j]++;
          }
          else {
	    error_solve_T_profile(T[j], T[j+1], T[j-1], T0[j], moist[j], 
				  max_moist[j], bubble[j], expt[j], ice[j], 
				  gamma[j-1], A[j], B[j], C[j], D[j], 
				  E[j], ErrorString);
            return ( ERROR );
	  }
	}
      }
      
      diff=fabs(oldT-T[j]);
      if(diff > maxdiff) maxdiff=diff;
    }
    
    if(NOFLUX) { 
      /** Solve for bottom temperature if using no flux lower boundary **/
      j = Nnodes-1;
      oldT=T[j];
      
      if(T[j] >= 0 || !FS_ACTIVE || !options.FROZEN_SOIL) {
	if(!EXP_TRANS )
	  T[j] = (A[j]*T0[j]+B[j]*(T[j]-T[j-1])+C[j]*T[j]+D[j]*T[j-1]+E[j]*(0.-ice[j]))/(A[j]+C[j]+D[j]);
	else
	  T[j] = (A[j]*T0[j]+B[j]*(T[j]-T[j-1])+C[j]*(T[j]+T[j-1])-D[j]*(T[j]-T[j-1])+E[j]*(0.-ice[j]))/(A[j]+2.*C[j]);
      }
      else {
	T[Nnodes-1] = root_brent(T0[Nnodes-1]-SOIL_DT, T0[Nnodes-1]+SOIL_DT,
				 ErrorString, soil_thermal_eqn, T[Nnodes-1],
				 T[Nnodes-2], T0[Nnodes-1], 
				 moist[Nnodes-1], max_moist[Nnodes-1], 
				 bubble[j], expt[Nnodes-1], 
				 ice[Nnodes-1], 
				 gamma[Nnodes-2], 
				 A[j], B[j], C[j], D[j], E[j], EXP_TRANS, j);
	if(T[j] <= -998 ) {
          if (options.TFALLBACK) {
            T[j] = T0[j];
            Tfbflag[j] = 1;
            Tfbcount[j]++;
          }
          else {
	    error_solve_T_profile(T[Nnodes-1], T[Nnodes-1],
				  T[Nnodes-2], T0[Nnodes-1], 
				  moist[Nnodes-1], max_moist[Nnodes-1], 
				  bubble[Nnodes-1], 
				  expt[Nnodes-1], ice[Nnodes-1], 
				  gamma[Nnodes-2], 
				  A[j], B[j], C[j], D[j], E[j], ErrorString);
            return ( ERROR );
          }
        }
      }
      
      diff=fabs(oldT-T[Nnodes-1]);
      if(diff>maxdiff) maxdiff=diff;
    }
    
    if(maxdiff <= threshold) Done=TRUE;
    
  }
  
  if (options.TFALLBACK) {
    // HACK to prevent runaway cold nose
    // Handle the case in which the a node was colder than both the nodes above and below
    // in the last time step, and that both differences have increased between the last
    // time step and the current one.
    for(j=1;j<Nnodes-1;j++) {
      if ( ( Tlast[j-1]-Tlast[j] > 0 && Tlast[j+1]-Tlast[j] > 0
             && (T[j-1]-T[j]) - (Tlast[j-1]-Tlast[j]) > 0
             && (T[j+1]-T[j]) - (Tlast[j+1]-Tlast[j]) > 0 )
        || ( Tlast[j-1]-Tlast[j] < 0 && Tlast[j+1]-Tlast[j] < 0
             && (T[j-1]-T[j]) - (Tlast[j-1]-Tlast[j]) < 0
             && (T[j+1]-T[j]) - (Tlast[j+1]-Tlast[j]) < 0 ) ) {
        T[j] = 0.5*(T[j-1]+T[j+1]); // crude fix for now; just average the T's without taking distance, conductivities into account
        Tfbflag[j] = 1;
        Tfbcount[j]++;
      }
    }
  }

  if(!Done && !Error) {
    if (options.TFALLBACK) {
      for(j=0;j<Nnodes;j++) {
        T[j] = T0[j];
        Tfbflag[j] = 1;
        Tfbcount[j]++;
      }
    }
    else {
      fprintf(stderr,"ERROR: Temperature Profile Unable to Converge!!!\n");
      fprintf(stderr,"Dumping Profile Temperatures (last, new).\n");
      for(j=0;j<Nnodes;j++) fprintf(stderr,"%f\t%f\n",T0[j],T[j]);
      fprintf(stderr,"ERROR: Cannot solve temperature profile:\n\tToo Many Iterations in solve_T_profile\n");
      return ( ERROR );
    }
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
  A         = (double) va_arg(ap, double);
  B         = (double) va_arg(ap, double);
  C         = (double) va_arg(ap, double);
  D         = (double) va_arg(ap, double);
  E         = (double) va_arg(ap, double);
  ErrorString = (char *) va_arg(ap, char *);
  
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
  fprintf(stderr,"A\t%f\n",A);
  fprintf(stderr,"B\t%f\n",B);
  fprintf(stderr,"C\t%f\n",C);
  fprintf(stderr,"D\t%f\n",D);
  fprintf(stderr,"E\t%f\n",E);

  fprintf(stderr,"Finished dumping values for solve_T_profile.\nTry increasing SOIL_DT to get model to complete cell.\nThen check output for instabilities.\n");

  return(ERROR);

}

#undef MAXIT



void fda_heat_eqn(double T_2[], double res[], int n, int init, ...)
{
  /**********************************************************************
  Heat Equation for implicit scheme (used to calculate residual of the heat equation)
  passed from solve_T_profile_implicit
  By Ming Pan, mpan@Princeton.EDU
 
  Modifications:
  2006-Aug-08 Integrated with 4.1.0 (from 4.0.3).				JCA
  2006-Aug-08 replaced second term of heat flux with an alternate form.		JCA
  2006-Aug-08 Added NOFLUX option.						JCA
  2006-Aug-08 Added EXP_TRANS option - this allows the heat flux equations
	      to be solved on a transformed grid using an exponential
	      distribution of node depths.					JCA
  2006-Aug-09 Included terms needed for Cs and kappa nodal updating for new
	      ice.								JCA
  2006-Aug-11 Included additional term in the storage term to account for
	      time-varying changes in Cs.					JCA
  2007-Aug-08 Added EXCESS_ICE option.						JCA
  2007-Oct-08 Fixed error in EXP_TRANS formulation.				JCA
  2011-Jun-03 Added options.ORGANIC_FRACT.  Soil properties now take
	      organic fraction into account.					TJB
  2011-Jun-10 Added bulk_dens_min and soil_dens_min to arglist of
	      soil_conductivity() to fix bug in commputation of kappa.		TJB
  2012-Jan-28 Removed restriction of cold nose fix to just top two nodes;
	      now all nodes are checked and corrected if necessary.		TJB
  2013-Jan-08 Excluded bottom node from check in cold nose fix.			TJB
  2013-Dec-26 Removed EXCESS_ICE option.				TJB
  **********************************************************************/
    
  static double  deltat;
  static int     FS_ACTIVE;
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
  static double Dp;
  static double *bulk_dens_min;
  static double *soil_dens_min;
  static double *quartz;
  static double *bulk_density;
  static double *soil_density;
  static double *organic;
  static double *depth;
  static int Nlayers;
  
  // variables used to calculate residual of the heat equation
  // defined here
  static double Ts;
  static double Tb;
  
  // locally used variables
  static double ice_new[MAX_NODES], Cs_new[MAX_NODES], kappa_new[MAX_NODES];
  static double DT[MAX_NODES],DT_down[MAX_NODES],DT_up[MAX_NODES],T_up[MAX_NODES];
  static double Dkappa[MAX_NODES];
  static double Bexp;
  char PAST_BOTTOM;
  double storage_term, flux_term, phase_term, flux_term1, flux_term2;
  double Lsum;
  int i, lidx;
  int focus, left, right;
  
  // argument list handling
  va_list arg_addr;

  // initialize variables if init==1  
  if (init==1) {
    va_start(arg_addr, init);
    deltat     = va_arg(arg_addr, double);
    FS_ACTIVE  = va_arg(arg_addr, int);
    NOFLUX     = va_arg(arg_addr, int);
    EXP_TRANS  = va_arg(arg_addr, int);
    T0         = va_arg(arg_addr, double *);
    moist      = va_arg(arg_addr, double *);
    ice        = va_arg(arg_addr, double *);
    kappa      = va_arg(arg_addr, double *);
    Cs         = va_arg(arg_addr, double *);
    max_moist  = va_arg(arg_addr, double *);
    bubble     = va_arg(arg_addr, double *);
    expt       = va_arg(arg_addr, double *);
    alpha      = va_arg(arg_addr, double *);
    beta       = va_arg(arg_addr, double *);
    gamma      = va_arg(arg_addr, double *);
    Zsum       = va_arg(arg_addr, double *);
    Dp         = va_arg(arg_addr, double);
    bulk_dens_min = va_arg(arg_addr, double *);
    soil_dens_min = va_arg(arg_addr, double *);
    quartz     = va_arg(arg_addr, double *);
    bulk_density = va_arg(arg_addr, double *);
    soil_density = va_arg(arg_addr, double *);
    organic    = va_arg(arg_addr, double *);
    depth      = va_arg(arg_addr, double *);
    Nlayers    = va_arg(arg_addr, int);
    
    if(EXP_TRANS){
      if(!NOFLUX)
	Bexp = logf(Dp+1.)/(double)(n+1); 
      else
	Bexp = logf(Dp+1.)/(double)(n);
    }    

    Ts = T0[0];
    if(!NOFLUX)
      Tb = T0[n+1];
    else
      Tb = T0[n];
    for (i=0; i<n; i++) 
      T_2[i] = T0[i+1];    
  }
  
  // calculate residuals if init==0
  else {
    // get the range of columns to calculate
    va_start(arg_addr, init);
    focus = va_arg(arg_addr, int);
    
    // calculate all entries if focus == -1
    if (focus==-1) {
      
      lidx = 0;
      Lsum = 0.;
      PAST_BOTTOM = FALSE;

      for (i=0; i<n+1; i++) {
	kappa_new[i]=kappa[i];
	if(i>=1) {  //all but surface node
	  // update ice contents
	  if (T_2[i-1]<0) {
	    ice_new[i] = moist[i] - maximum_unfrozen_water(T_2[i-1], 
							   max_moist[i], bubble[i], expt[i]);
	    if (ice_new[i]<0) ice_new[i]=0;
	  }
	  else ice_new[i] = 0;
	  Cs_new[i]=Cs[i];

	  // update other states due to ice content change
	  /***********************************************/
	  if (ice_new[i]!=ice[i]) {
	    kappa_new[i] = soil_conductivity(moist[i], moist[i] - ice_new[i],
					     soil_dens_min[lidx], bulk_dens_min[lidx], quartz[lidx],
					     soil_density[lidx], bulk_density[lidx], organic[lidx]);
	    Cs_new[i] = volumetric_heat_capacity(bulk_density[lidx]/soil_density[lidx], moist[i]-ice_new[i], ice_new[i], organic[lidx]);
	  }
	  /************************************************/	  
	}
	
	if(Zsum[i] > Lsum + depth[lidx] && !PAST_BOTTOM) {
	  Lsum += depth[lidx];
	  lidx++;
	  if( lidx == Nlayers ) {
	    PAST_BOTTOM = TRUE;
	    lidx = Nlayers-1;
	  }
	}
      }
      
      // constants used in fda equation
      for (i=0; i<n; i++) {
	if (i==0) {
	  DT[i]=T_2[i+1]-Ts;
	  DT_up[i]=T_2[i]-Ts;
	  DT_down[i]=T_2[i+1]-T_2[i];
	  T_up[i]=Ts;
	}
	else if (i==n-1) {
	  DT[i]=Tb-T_2[i-1];
	  DT_up[i]=T_2[i]-T_2[i-1];
	  DT_down[i]=Tb-T_2[i];
	  T_up[i]=T_2[i-1];
	}
	else {
	  DT[i]=T_2[i+1]-T_2[i-1];
	  DT_up[i]=T_2[i]-T_2[i-1];
	  DT_down[i]=T_2[i+1]-T_2[i];
	  T_up[i]=T_2[i-1];
	}
	if(i<n-1)
	  Dkappa[i]=kappa_new[i+2]-kappa_new[i];
	else
	  if(!NOFLUX)
	    Dkappa[i]=kappa_new[i+2]-kappa_new[i];
	  else
	    Dkappa[i]=kappa_new[i+1]-kappa_new[i];
      }
      
      for (i=0; i<n; i++) {
	storage_term = Cs_new[i+1]*(T_2[i] - T0[i+1])/deltat + T_2[i]*(Cs_new[i+1]-Cs[i+1])/deltat;
	if(!EXP_TRANS) {
	  flux_term1 = Dkappa[i]/alpha[i]*DT[i]/alpha[i];
	  flux_term2 = kappa_new[i+1]*(DT_down[i]/gamma[i]-DT_up[i]/beta[i])/(0.5*alpha[i]);
	}
	else { //grid transformation
	  flux_term1 = Dkappa[i]/2.*DT[i]/2./(Bexp*(Zsum[i+1]+1.))/(Bexp*(Zsum[i+1]+1.));
	  flux_term2 = kappa_new[i+1]*((DT_down[i]-DT_up[i])/(Bexp*(Zsum[i+1]+1.))/(Bexp*(Zsum[i+1]+1.))  -  DT[i]/2./(Bexp*(Zsum[i+1]+1.)*(Zsum[i+1]+1.)));
	}
	//inelegant fix for "cold nose" problem - when a very cold node skates off to
	//much colder and breaks the second law of thermodynamics (because
	//flux_term1 exceeds flux_term2 in absolute magnitude) - therefore, don't let
	//that node get any colder.  This only seems to happen in the first and
	//second near-surface nodes.
//	if (i<n-1) {
//	  if(fabs(DT[i])>5. && (T_2[i]<T_2[i+1] && T_2[i]<T_up[i])){//cold nose
//	    if((flux_term1<0 && flux_term2>0) && fabs(flux_term1)>fabs(flux_term2)){
//	      flux_term1 = 0;
//#if VERBOSE
//	      fprintf(stderr,"WARNING: resetting thermal flux term in soil heat solution to zero for node %d.\nT[i]=%.2f T[i-1]=%.2f T[i+1]=%.2f flux_term1=%.2f flux_term2=%.2f\n",i+1,T_2[i],T_up[i],T_2[i+1],flux_term1,flux_term2);
//#endif
//	    }
//	  }
//	}
	flux_term = flux_term1+flux_term2;
	phase_term   = ice_density*Lf * (ice_new[i+1] - ice[i+1])/deltat;
        res[i] = flux_term + phase_term - storage_term;
      }
    }
    
    // only calculate entries focus-1, focus, and focus+1 if focus has a value>=0
    else {
      if (focus==0)    left=0;   else  left=focus-1;
      if (focus==n-1) right=n-1; else right=focus+1;

      // update ice content for node focus and its adjacents
      for (i=left; i<=right; i++) {
	if (T_2[i]<0) {
	  ice_new[i+1] = moist[i+1] - maximum_unfrozen_water(T_2[i], 
							     max_moist[i+1], bubble[i+1], expt[i+1]);
	  if (ice_new[i+1]<0) ice_new[i+1]=0;
	}
	else ice_new[i+1]=0;
      }
      
      // update other parameters due to ice content change
      /********************************************************/
      lidx = 0;
      Lsum = 0.;
      PAST_BOTTOM = FALSE;
      for (i=0; i<=right+1; i++) {
	if(i>=left+1) {
	  if (ice_new[i]!=ice[i]) {
	    kappa_new[i] = soil_conductivity(moist[i], moist[i] - ice_new[i],
					     soil_dens_min[lidx], bulk_dens_min[lidx], quartz[lidx],
					     soil_density[lidx], bulk_density[lidx], organic[lidx]);
	    Cs_new[i] = volumetric_heat_capacity(bulk_density[lidx]/soil_density[lidx], moist[i]-ice_new[i], ice_new[i], organic[lidx]);
	  }
	}
	if(Zsum[i] > Lsum + depth[lidx] && !PAST_BOTTOM) {
	  Lsum += depth[lidx];
	  lidx++;
	  if( lidx == Nlayers ) {
	    PAST_BOTTOM = TRUE;
	    lidx = Nlayers-1;
	  }
	}
      }
      /*********************************************************/
      
      // update other states due to ice content change
      for (i=left; i<=right; i++) {
	if (i==0) {
	  DT[i]=T_2[i+1]-Ts;
	  DT_up[i]=T_2[i]-Ts;
	  DT_down[i]=T_2[i+1]-T_2[i];
	  T_up[i]=Ts;
	}
	else if (i==n-1) {
	  DT[i]=Tb-T_2[i-1];
	  DT_up[i]=T_2[i]-T_2[i-1];
	  DT_down[i]=Tb-T_2[i];
	  T_up[i]=T_2[i-1];
	}
	else {
	  DT[i]=T_2[i+1]-T_2[i-1];
	  DT_up[i]=T_2[i]-T_2[i-1];
	  DT_down[i]=T_2[i+1]-T_2[i];
	  T_up[i]=T_2[i-1];
	}
	//update Dkappa due to ice content change
	/*******************************************/
	if(i<n-1)
	  Dkappa[i]=kappa_new[i+2]-kappa_new[i];
	else
	  if(!NOFLUX)
	    Dkappa[i]=kappa_new[i+2]-kappa_new[i];
	  else
	    Dkappa[i]=kappa_new[i+1]-kappa_new[i];
	/********************************************/
      }
      
      for (i=left; i<=right; i++) {
	storage_term = Cs_new[i+1]*(T_2[i] - T0[i+1])/deltat + T_2[i]*(Cs_new[i+1]-Cs[i+1])/deltat;
	if(!EXP_TRANS) {
	  flux_term1 = Dkappa[i]/alpha[i]*DT[i]/alpha[i];
	  flux_term2 = kappa_new[i+1]*(DT_down[i]/gamma[i]-DT_up[i]/beta[i])/(0.5*alpha[i]);
	}
	else { //grid transformation
	  flux_term1 = Dkappa[i]/2.*DT[i]/2./(Bexp*(Zsum[i+1]+1.))/(Bexp*(Zsum[i+1]+1.));
	  flux_term2 = kappa_new[i+1]*((DT_down[i]-DT_up[i])/(Bexp*(Zsum[i+1]+1.))/(Bexp*(Zsum[i+1]+1.))  -  DT[i]/2./(Bexp*(Zsum[i+1]+1.)*(Zsum[i+1]+1.)));
	}
	//inelegant fix for "cold nose" problem - when a very cold node skates off to
	//much colder and breaks the second law of thermodynamics (because
	//flux_term1 exceeds flux_term2 in absolute magnitude) - therefore, don't let
	//that node get any colder.  This only seems to happen in the first and
	//second near-surface nodes.
//	if (i<n-1) {
//	  if(fabs(DT[i])>5. && (T_2[i]<T_2[i+1] && T_2[i]<T_up[i])){//cold nose
//	    if((flux_term1<0 && flux_term2>0) && fabs(flux_term1)>fabs(flux_term2)){
//	      flux_term1 = 0;
//#if VERBOSE
//	      fprintf(stderr,"WARNING: resetting thermal flux term in soil heat solution to zero for node %d.\nT[i]=%.2f T[i-1]=%.2f T[i+1]=%.2f flux_term1=%.2f flux_term2=%.2f\n",i+1,T_2[i],T_up[i],T_2[i+1],flux_term1,flux_term2);
//#endif
//	    }
//	  }
	flux_term = flux_term1+flux_term2;
	phase_term   = ice_density*Lf * (ice_new[i+1] - ice[i+1]) / deltat;
        res[i] = flux_term + phase_term - storage_term;
      }
    } // end of calculation of focus node only
  } // end of non-init
  
}
