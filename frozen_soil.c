#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <vicNl.h>

#define MAXIT 1000

void frozen_soil(soil_con_struct soil_con,
                 layer_data_struct *layer,
		 energy_bal_struct *energy,
		 int rec,
                 int veg,
                 int dt,
                 double surf_temp,
                 int Ulayer,
                 int Llayer) {
/**********************************************************************
	frozen_soil	Keith Cherkauer		May 18, 1996

  This subroutine determines if there is a frozen soil layer, and then
  computes the hydrologic properties of the layer.

  UNITS: m, kg, s, degrees C (unless otherwise noted).

         moisture is returned as (mm) of water in column.

**********************************************************************/
  extern option_struct options;
  extern debug_struct debug;

  char ErrorString[MAXSTRING];
  int i, j, index;
  int Tlayer;
  double *kappa;	/* thermal conductivity of layer (W/m/K) */
  double *T, *Tnew;	/* soil layer temperature: T[0] = surface, T[1] = 0C */
  double *dz;		/* depth between layers */
  double *Cs;		/* volumetric heat capacity of layer (J/m^3/K) */
  double *moist;	/* layer unfrozen water content (mm) */
  double *expt;
  double *max_moist;
  double unfrozen;	/* total amount of water in layer that can exist
			   in an unfrozen condition */
  double frozen;	/* total amount of water in layer that is frozen */
  double fdepth;	/* freezing front depth */
  double tdepth;	/* thawing front depth */
  double *old_fdepth, *old_tdepth;
  double deltat;

  deltat = (double)dt*3600.;	/* convert time step to seconds */
  Tlayer = Ulayer+2;

  dz = (double *)calloc(Tlayer,sizeof(double));
  for(i=0;i<Tlayer;i++) dz[i] = energy->dz[i];

  /** Initialize Soil Layers **/
  T = (double *)calloc(Tlayer,sizeof(double));

  fdepth = energy->fdepth[0];
  tdepth = energy->fdepth[1];

  for(i=0;i<Tlayer;i++) T[i] = energy->T[i];

  kappa = (double *)calloc(Tlayer,sizeof(double));
  Cs = (double *)calloc(Tlayer,sizeof(double));
  moist = (double *)calloc(Tlayer,sizeof(double));
  expt = (double *)calloc(Tlayer,sizeof(double));
  max_moist = (double *)calloc(Tlayer,sizeof(double));
  soil_thermal_calc(soil_con, layer, *energy, kappa,
      Cs, moist, expt, max_moist, options.Nlayer, Tlayer);

  /** Iterate Thermal Solution to Find Freezing Depth **/
  T[0] = surf_temp;
  Tnew = solve_T_profile(T,dz,kappa,Cs,moist,deltat,
         max_moist,soil_con.bubble,expt,energy->ice,Tlayer);
  for(i=0;i<Tlayer;i++) T[i] = Tnew[i];
  free((char *)Tnew);

  if(debug.PRT_BALANCE && debug.DEBUG) {
    printf("Record #%i\nAfter Soil Profile Solution\n",rec);
    write_layer(layer,veg,options.Nlayer,soil_con.depth);
  } 

  /** Calculate New Layer Depths **/
  old_fdepth = (double *)calloc(options.Nlayer,sizeof(double));
  old_tdepth = (double *)calloc(options.Nlayer,sizeof(double));
  for(i=0;i<options.Nlayer;i++) {
    old_fdepth[i] = layer[i].fdepth;
    old_tdepth[i] = layer[i].tdepth;
  }
  find_0_degree_fronts(energy, layer, soil_con.dp, soil_con.depth,
      T, options.Nlayer, Tlayer);

  /***** Find New Layer Temperatures *****/
  find_sublayer_temperatures(layer,T,dz,soil_con.depth,
      energy->fdepth[0],energy->fdepth[1],options.Nlayer,Tlayer);

  /** Store Layer Temperature Values **/
  for(i=0;i<Tlayer;i++) energy->T[i] = T[i];

  if(energy->fdepth[0]>0) energy->frozen = TRUE;
  else energy->frozen = FALSE;

  /** Redistribute Soil Properties for New Frozen Soil Layer Size **/

  if(energy->frozen || fdepth>0.0)
    redistribute_moisture(layer, energy->fdepth,
        soil_con.max_moist, old_fdepth, old_tdepth, soil_con.depth,
        options.Nlayer);
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
	layer[index].unfrozen = unfrozen;

	frozen = layer[index].moist_froz - unfrozen;
	if(frozen < 0.0) {
	  frozen = 0.0;
	  unfrozen = layer[index].moist_froz;
	}
	layer[index].ice = frozen;
	layer[index].moist_froz = unfrozen;
      }
      else if(layer[index].T_froz == 0.) {
	layer[index].unfrozen = soil_con.max_moist[index];
	layer[index].ice = 0.;
      }
      else if(layer[index].T_froz != -999.) {
        sprintf(ErrorString,"ERROR: Frozen Layer Temperature > 0C (%lf)",
            layer[index].T_froz);
	vicerror(ErrorString);
      }
    }
    else layer[index].ice=0.0;
  }

  for(j=0;j<Tlayer;j++) {
    if(T[j]<0.) {
      energy->ice[j] = moist[j] - maximum_unfrozen_water(T[j],max_moist[j],
               soil_con.bubble,expt[j]);
      if(energy->ice[j]<0.) energy->ice[j]=0.;
      if(energy->ice[j]>max_moist[j]) energy->ice[j]=max_moist[j];
    }
    else energy->ice[j]=0.;
  }

  if(debug.PRT_BALANCE && debug.DEBUG) {
    printf("After Refreezing Moisture\n");
    write_layer(layer,veg,options.Nlayer,soil_con.depth);
  } 

  free((char *)moist);
  free((char *)kappa);
  free((char *)Cs);
  free((char *)expt);
  free((char *)max_moist);
  kappa = NULL;
  Cs = NULL;
  moist = NULL;
  expt = NULL;
  max_moist = NULL;

  soil_thermal_calc(soil_con, layer, *energy, kappa,
      Cs, moist, expt, max_moist, options.Nlayer, Tlayer);

  /** Save Thermal Conductivities for Energy Balance **/
  /***** WARNING This will not work if dz or layers are changed *****/
  energy->kappa[0] = layer[0].kappa;
  energy->Cs[0] = layer[0].Cs;
  energy->kappa[1] = layer[1].kappa;
  energy->Cs[1] = layer[1].Cs;
  
  free((char *)T);
  free((char *)dz);

}

double *solve_T_profile(double *T0,
                        double *dz,
                        double *kappa,
                        double *Cs,
                        double *moist,
                        double deltat,
                        double *max_moist,
                        double bubble,
                        double *expt,
                        double *ice,
                        int Tlayer) {
/**********************************************************************
  This subroutine was written to iteratively solve the soil temperature
  profile using a numerical difference equation.  The solution equation
  is second order in space, and first order in time.
**********************************************************************/

  extern debug_struct debug;

  double *T;
  double maxdiff, diff;
  double threshold = 1.e-3;	/* temperature profile iteration threshold */
  double oldT;
  double alpha, beta, gamma, fprime;
  int j, Done, ItCount;

  T = (double *)calloc(Tlayer,sizeof(double));

  for(j=0;j<Tlayer;j++) T[j]=T0[j];

  Done = FALSE;
  ItCount = 0;

  while(!Done && ItCount<MAXIT) {
    ItCount++;
    maxdiff=threshold;
    for(j=1;j<Tlayer-1;j++) {
      oldT=T[j];
  
      /**	2nd order variable kappa equation **/
      alpha = ((dz[j+1] + dz[j]) / 2.0 + (dz[j] + dz[j-1]) / 2.0);
      beta = pow((dz[j+1]+dz[j])/2.0, 2.0) + pow((dz[j]+dz[j-1])/2.0, 2.0);
      gamma = ((dz[j+1] + dz[j]) / 2.0 - (dz[j] + dz[j-1]) / 2.0);
      fprime = (T[j+1]-T[j-1])/alpha;
      alpha = alpha*alpha;

      if(T[j]>=0) {
        T[j] = (beta*deltat*(kappa[j+1]-kappa[j-1])*(T[j+1]-T[j-1])
	     + 2.*alpha*deltat*kappa[j]*(T[j+1]+T[j-1]-gamma*fprime)
             + alpha*beta*Cs[j]*T0[j] 
             + alpha*beta*ice_density*Lf*(0.-ice[j])) 
             / (alpha*beta*Cs[j] + 4.*kappa[j]*alpha*deltat);
      }
      else {
        T[j] = root_brent(T0[j]-5.,T0[j]+5.,soil_thermal_eqn, T[j+1],
             T[j-1], T0[j], kappa[j], kappa[j+1], kappa[j-1], Cs[j], deltat,
             moist[j], max_moist[j], bubble, expt[j], ice[j],
             alpha, beta, gamma, fprime, 1.0);
      }

      diff=fabs(oldT-T[j]);
      if(diff>maxdiff) maxdiff=diff;
    }
    if(maxdiff<=threshold) Done=TRUE;

  }

  if(!Done) {
    fprintf(stderr,"ERROR: Temperature Profile Unable to Converge!!!\n");
    fprintf(stderr,"Dumping Profile Temperatures (last, new).\n");
    for(j=0;j<Tlayer;j++) fprintf(stderr,"%lf\t%lf\n",T0[j],T[j]);
    vicerror("ERROR: Cannot solve temperature profile:\n\tToo Many Iterations in solve_T_profile");
  }

  return(T);

}

#undef MAXIT
