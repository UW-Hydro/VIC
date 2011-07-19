#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

static char vcid[] = "$Id$";

double soil_conductivity(double moist, 
			 double Wu, 
			 double soil_dens_min, 
			 double bulk_dens_min,
			 double quartz,
			 double soil_density, 
			 double bulk_density,
			 double organic) {
/**********************************************************************
  Soil thermal conductivity calculated using Johansen's method.

  Reference: Farouki, O.T., "Thermal Properties of Soils" 1986
	Chapter 7: Methods for Calculating the Thermal Conductivity 
		of Soils

  H.B.H. - refers to the handbook of hydrology.

  porosity = n = porosity
  ratio = Sr = fractionaldegree of saturation
  All K values are conductivity in W/mK
  Wu is the fractional volume of unfrozen water

  UNITS: input in m, kg, s

  Returns K in W/m/K

  double moist         total moisture content (mm/mm)
  double Wu            liquid water content (mm/mm)
  double soil_dens_min mineral soil particle density (kg m-3)
  double bulk_dens_min mineral soil bulk density (kg m-3)
  double quartz        mineral soil quartz content (fraction of mineral soil volume)
  double soil_density  soil particle density (kg m-3)
  double bulk_density  soil bulk density (kg m-3)
  double organic       total soil organic content (fraction of total solid soil volume)
                         i.e., organic fraction of solid soil = organic*(1-porosity)
                               mineral fraction of solid soil = (1-organic)*(1-porosity)

  Modifications:

  2011-Jun-03 Added options.ORGANIC_FRACT.  Soil properties now take
	      organic fraction into account.					TJB
  2011-Jun-10 Added bulk_dens_min and soil_dens_min to arglist of
	      soil_conductivity() to fix bug in commputation of kappa.		TJB
**********************************************************************/
  double Ke;
  double Ki = 2.2;      /* thermal conductivity of ice (W/mK) */
  double Kw = 0.57;     /* thermal conductivity of water (W/mK) */
  double Ksat;
  double Kdry;          /* Dry thermal conductivity of soil (W/mK), including mineral and organic fractions */
  double Kdry_org = 0.05; /* Dry thermal conductivity of organic fraction (W/mK) (Farouki 1981) */
  double Kdry_min;      /* Dry thermal conductivity of mineral fraction (W/mK) */
  double Ks;            /* thermal conductivity of solid (W/mK), including mineral and organic fractions */
  double Ks_org = 0.25; /* thermal conductivity of organic fraction of solid (W/mK) (Farouki 1981) */
  double Ks_min;        /* thermal conductivity of mineral fraction of solid (W/mK) */
  double Sr;            /* fractional degree of saturation */
  double K;
  double porosity;

  /* Calculate dry conductivity as weighted average of mineral and organic fractions. */
  Kdry_min = (0.135*bulk_dens_min+64.7)/(soil_dens_min-0.947*bulk_dens_min);
  Kdry = (1-organic)*Kdry_min + organic*Kdry_org;

  if(moist>0.) {

    porosity = 1.0 - bulk_density / soil_density; //NOTE: if excess_ice present,
                                                  //this is actually effective_porosity

    Sr = moist/porosity;

    // Compute Ks of mineral soil; here "quartz" is the fraction (quartz volume / mineral soil volume)
    if(quartz < .2)
      Ks_min = pow(7.7,quartz) * pow(3.0,1.0-quartz);  // when quartz is less than 0.2
    else
      Ks_min = pow(7.7,quartz) * pow(2.2,1.0-quartz);  // when quartz is greater than 0.2
    Ks = (1-organic)*Ks_min + organic*Ks_org;

    if(Wu==moist) {

      /** Soil unfrozen **/
      Ksat = pow(Ks,1.0-porosity) * pow(Kw,porosity);
      Ke = 0.7 * log10(Sr) + 1.0;

    }
    else {

      /** Soil frozen **/
      Ksat = pow(Ks,1.0-porosity) * pow(Ki,porosity-Wu) * pow(Kw,Wu);
      Ke = Sr;

    }

    K = (Ksat-Kdry)*Ke+Kdry;
    if(K<Kdry) K=Kdry;

  }
  else K=Kdry;

  return (K);
}


double volumetric_heat_capacity(double soil_fract,
                                double water_fract,
                                double ice_fract,
                                double organic_fract) {
/**********************************************************************
  This subroutine calculates the soil volumetric heat capacity based 
  on the fractional volume of its component parts.

  Constant values are volumetric heat capacities in J/m^3/K

  double soil_fract    fraction of soil volume composed of solid soil (fract)
  double organic_fract fraction of solid soil volume composed of organic matter (fract)
  double water_fract   fraction of soil volume composed of liquid water (fract)
  double ice_fract     fraction of soil volume composed of ice (fract)

  Modifications:

  2011-Jun-03 Added options.ORGANIC_FRACT.  Soil properties now take
	      organic fraction into account.					TJB
**********************************************************************/

  double Cs;

  Cs  = 2.0e6 * soil_fract*(1-organic_fract);
  Cs += 2.7e6 * soil_fract*organic_fract;
  Cs += 4.2e6 * water_fract;
  Cs += 1.9e6 * ice_fract;
  Cs += 1.3e3 * (1. - (soil_fract + water_fract + ice_fract)); // air

  return (Cs);

}


void set_node_parameters(double   *dz_node,
			 double   *Zsum_node,
			 double   *max_moist_node,
			 double   *expt_node,
			 double   *bubble_node,
			 double   *alpha,
			 double   *beta,
			 double   *gamma,
			 double   *depth,
			 double   *max_moist,
			 double   *expt,
			 double   *bubble,
			 double   *quartz,
#if QUICK_FS
			 double ***ufwc_table_node,
#endif
#if EXCESS_ICE
			 double    *porosity,
			 double    *effective_porosity,
			 double    *porosity_node,
			 double    *effective_porosity_node,
#endif
			 int       Nnodes,
			 int       Nlayers,
			 char      FS_ACTIVE) {
/**********************************************************************
  This subroutine sets the thermal node soil parameters to constant 
  values based on those defined for the current grid cells soil type.
  Thermal node propertiers for the energy balance solution are also 
  set (these constants are used to reduce the solution time required
  within each iteration).

  double   *dz_node          thermal node thicknes (m)
  double   *Zsum_node        thermal node depth (m)
  double   *max_moist_node   maximum moisture content at thermal node (mm/mm)
  double   *expt_node        exponential at thermal node ()
  double   *bubble_node      bubbling pressure at thermal node (cm)
  double   *alpha            first thermal eqn term ()
  double   *beta             second thermal eqn term ()
  double   *gamma            third thermal eqn term ()
  double   *depth            soil moisture layer thickness (m)
  double   *max_moist        soil moisture layer maximum moisture content (mm)
  double   *bubble           soil moisture layer bubbling pressure (cm)
  double    quartz           soil quartz content (fract)
  double ***ufwc_table_node  table of unfrozen water contents ()
  double   *porosity         soil layer porosity 
  double   *effective_porosity   effective soil layer porosity
  double   *porosity_node       thermal node porosity
  double   *effective_porosity_node  effective thermal node porosity
  int       Nnodes           number of soil thermal nodes
  int       Nlayers          number of soil moisture layers
  char      FS_ACTIVE        TRUE if frozen soils are active in grid cell

  Modifications:

  02-11-00 Modified to remove node zone averages, node parameters
           are now set based on the parameter value of the layer 
           in which they fall.  Averaging of layer properties 
	   only occurs if the node falls directly on a layer
	   boundary.                                        KAC
  11-20-02 Modified to correct differences between calculations
           to determine maximum node moisture and node moisture,
           so that nodes on the boundary between soil layers are
           computed the same way for both.                  KAC
  2007-Apr-24 Passing in Zsum_node rather than recalculating.		JCA
  2007-Apr-24 Rearranged terms in finite-difference heat equation
	      (equation 8 of Cherkauer et al. (1999)).  See note
	      in solve_T_profile.  This affects the equations for
	      beta and gamma.						JCA
  2007-Aug-09 Added features for EXCESS_ICE option.			JCA
  2009-Jul-31 Removed unused layer_node_fract array.			TJB
**********************************************************************/

  extern option_struct options;
#if QUICK_FS
  extern double temps[];
#endif

  char   PAST_BOTTOM;
  int    nidx, lidx;
  double Lsum; /* cumulative depth of moisture layer */
  double Zsum; /* upper boundary of node thermal layer */
  double deltaL[MAX_LAYERS+1];
#if QUICK_FS
  int    ii;
  double Aufwc;
  double Bufwc;
#endif

  PAST_BOTTOM = FALSE;
  lidx = 0;
  Lsum = 0.;

  /* set node parameters */
  for(nidx=0;nidx<Nnodes;nidx++) {
    
    if(Zsum_node[nidx] == Lsum + depth[lidx] && nidx != 0 && lidx != Nlayers-1) {
      /* node on layer boundary */
      max_moist_node[nidx] = (max_moist[lidx] / depth[lidx] 
			      + max_moist[lidx+1] / depth[lidx+1]) / 1000 / 2.;
      expt_node[nidx]      = (expt[lidx] + expt[lidx+1]) / 2.;
      bubble_node[nidx]    = (bubble[lidx] + bubble[lidx+1]) / 2.;
#if EXCESS_ICE
      porosity_node[nidx]    = (porosity[lidx] + porosity[lidx+1]) / 2.;
      effective_porosity_node[nidx]    = (effective_porosity[lidx] + effective_porosity[lidx+1]) / 2.;
#endif
    }
    else { 
      /* node completely in layer */
      max_moist_node[nidx] = max_moist[lidx] / depth[lidx] / 1000;
      expt_node[nidx]      = expt[lidx];
      bubble_node[nidx]    = bubble[lidx];
#if EXCESS_ICE
      porosity_node[nidx]  = porosity[lidx];
      effective_porosity_node[nidx] = effective_porosity[lidx];
#endif
    }
    if(Zsum_node[nidx] > Lsum + depth[lidx] && !PAST_BOTTOM) {
      Lsum += depth[lidx];
      lidx++;
      if( lidx == Nlayers ) {
	PAST_BOTTOM = TRUE;
	lidx = Nlayers-1;
      }
    }

  }

  /* set constant variables for thermal calculations */
  for(nidx=0;nidx<Nnodes-2;nidx++) {
    alpha[nidx] = Zsum_node[nidx+2] - Zsum_node[nidx];
    beta[nidx] = Zsum_node[nidx+1] - Zsum_node[nidx];
    gamma[nidx] = Zsum_node[nidx+2] - Zsum_node[nidx+1];
  }
  if(options.NOFLUX) {
    /* no flux bottom boundary activated */
    alpha[Nnodes-2] = 2. * (Zsum_node[Nnodes-1] - Zsum_node[Nnodes-2]);
    beta[nidx] = Zsum_node[Nnodes-1] - Zsum_node[Nnodes-2];
    gamma[nidx] = Zsum_node[Nnodes-1] - Zsum_node[Nnodes-2];
  }


#if QUICK_FS

  /* If quick frozen soil solution activated, prepare a linearized
     estimation on the maximum unfrozen water content equation */

  if(FS_ACTIVE && options.FROZEN_SOIL) {
    for(nidx=0;nidx<Nnodes;nidx++) { 
      for(ii=0;ii<QUICK_FS_TEMPS;ii++) {
	Aufwc = maximum_unfrozen_water(temps[ii], 1.0, 
				       bubble_node[nidx], 
				       expt_node[nidx]);
	Bufwc = maximum_unfrozen_water(temps[ii+1], 1.0, 
				       bubble_node[nidx], 
				       expt_node[nidx]);
	ufwc_table_node[nidx][ii][0] 
	  = linear_interp(0., temps[ii], temps[ii+1], Aufwc, Bufwc);
	ufwc_table_node[nidx][ii][1] 
	  = (Bufwc - Aufwc) / (temps[ii+1] - temps[ii]);
      }
    }
  }
#endif
}

#define N_INTS 5

int distribute_node_moisture_properties(double *moist_node,
					double *ice_node,
					double *kappa_node,
					double *Cs_node,
					double *Zsum_node,
					double *T_node,
					double *max_moist_node,
#if QUICK_FS
					double ***ufwc_table_node,
#else
					double *expt_node,
					double *bubble_node,
#endif
#if EXCESS_ICE
					double *porosity_node,
					double *effective_porosity_node,
#endif
					double *moist,
					double *depth,
					double *soil_dens_min,
					double *bulk_dens_min,
					double *quartz,
					double *soil_density,
					double *bulk_density,
					double *organic,
					int     Nnodes,
					int     Nlayers,
					char    FS_ACTIVE) {
  /*********************************************************************
  This subroutine determines the moisture and ice contents of each 
  soil thermal node based on the current node temperature and layer
  moisture content.  Thermal conductivity and volumetric heat capacity
  are then estimated for each node based on the division of moisture 
  contents..

  double *moist_node      thermal node moisture content (mm/mm)
  double *ice_node        thermal node ice content (mm/mm)
  double *kappa_node      thermal node thermal conductivity (W m-1 K-1)
  double *Cs_node         thermal node heat capacity (J m-3 K-1)
  double *Zsum_node       thermal node depth (m)
  double *T_node          thermal node temperature (C)
  double *max_moist_node  thermal node maximum moisture content (mm/mm)
  double *expt_node       thermal node exponential
  double *bubble_node     thermal node bubbling pressure (cm)
  double *moist           soil layer moisture (mm)
  double *depth           soil layer thickness (m)
  double *soil_dens_min   mineral soil particle density (kg m-3)
  double *bulk_dens_min   mineral soil bulk density (kg m-3)
  double  quartz          mineral soil quartz content (fract)
  double *soil_density    soil particle density (kg m-3)
  double *bulk_density    soil bulk density (kg m-3)
  double  organic         soil organic content (fract)
  int     Nnodes          number of soil thermal nodes
  int     Nlayers         number of soil moisture layers

  Modifications:

  02-11-00 Modified to remove node zone averages, node parameters
           are now set based on the parameter value of the layer 
           in which they fall.  Averaging of layer properties 
	   only occurs if the node falls directly on a layer
	   boundary.                                        KAC
  11-20-02 Modified to check that node soil moisture is less than
           or equal to maximum node soil moisture, otherwise an 
           error is printed to the screen and the model exits.  KAC
  2005-Mar-24 Removed abs() from check on soil moisture.		TJB
  2007-Apr-04 Modified to handle grid cell errors by returning to
              the main subroutine, rather than ending the simulation.	GCT/KAC
  2007-Apr-24 Passing in Zsum_node rather than recalculating.		JCA
  2007-Aug-09 Added features for EXCESS_ICE.				JCA
  2009-Feb-09 Removed dz_node from call to
	      distribute_node_moisture_properties.			KAC via TJB
  2010-Nov-02 Turned off error reporting for when node moisture exceeds
	      maximum; now node moisture greater than maximum is simply
	      reset to maximum.						TJB
  2011-Jun-03 Added options.ORGANIC_FRACT.  Soil properties now take
	      organic fraction into account.				TJB
  2011-Jun-10 Added bulk_dens_min and soil_dens_min to arglist of
	      soil_conductivity() to fix bug in commputation of kappa.		TJB

*********************************************************************/

  extern option_struct options;

  char PAST_BOTTOM;
  int nidx, lidx;
  double Lsum; /* cumulative depth of moisture layer */
  double soil_fract;

  lidx = 0;
  Lsum = 0.;
  PAST_BOTTOM = FALSE;

  /* node estimates */
  for(nidx=0;nidx<Nnodes;nidx++) {
 
    if(Zsum_node[nidx] == Lsum + depth[lidx] && nidx != 0 && lidx != Nlayers-1) {
      /* node on layer boundary */
      moist_node[nidx] = (moist[lidx] / depth[lidx] 
			      + moist[lidx+1] / depth[lidx+1]) / 1000 / 2.;
      soil_fract = (bulk_density[lidx] / soil_density[lidx] 
		    + bulk_density[lidx+1] / soil_density[lidx+1]) / 2.;
    }
    else { 
      /* node completely in layer */
      moist_node[nidx] = moist[lidx] / depth[lidx] / 1000;
      soil_fract = (bulk_density[lidx] / soil_density[lidx]);
    }      


    // Check that node moisture does not exceed maximum node moisture
//    if (moist_node[nidx]-max_moist_node[nidx] > 0.0001) {
//      fprintf( stderr, "Node soil moisture, %f, exceeds maximum node soil moisture, %f.\n", 
//	       moist_node[nidx], max_moist_node[nidx] );
//      return( ERROR );
//    }
    if (moist_node[nidx]-max_moist_node[nidx] > 0) moist_node[nidx] = max_moist_node[nidx]; // HACK!!!!!!!!!!!

    if(T_node[nidx] < 0 && (FS_ACTIVE && options.FROZEN_SOIL)) {
      /* compute moisture and ice contents */
#if QUICK_FS
      ice_node[nidx] 
	= moist_node[nidx] - maximum_unfrozen_water_quick(T_node[nidx],
						   max_moist_node[nidx], 
						   ufwc_table_node[nidx]);
#else
      ice_node[nidx] 
	= moist_node[nidx] - maximum_unfrozen_water(T_node[nidx],
#if EXCESS_ICE
						    porosity_node[nidx],
						    effective_porosity_node[nidx],
#endif
						    max_moist_node[nidx], 
						    bubble_node[nidx],
						    expt_node[nidx]);
#endif
      if(ice_node[nidx]<0) ice_node[nidx]=0;

      /* compute thermal conductivity */
      kappa_node[nidx] 
	= soil_conductivity(moist_node[nidx], moist_node[nidx] - ice_node[nidx], 
			    soil_dens_min[lidx], bulk_dens_min[lidx], quartz[lidx],
			    soil_density[lidx], bulk_density[lidx], organic[lidx]);

    }
    else {
      /* compute moisture and ice contents */
      ice_node[nidx]   = 0;
      /* compute thermal conductivity */
      kappa_node[nidx] 
	= soil_conductivity(moist_node[nidx], moist_node[nidx], 
			    soil_dens_min[lidx], bulk_dens_min[lidx], quartz[lidx],
			    soil_density[lidx], bulk_density[lidx], organic[lidx]);
    }
    /* compute volumetric heat capacity */
    Cs_node[nidx] = volumetric_heat_capacity(bulk_density[lidx]/soil_density[lidx],
					     moist_node[nidx] - ice_node[nidx], ice_node[nidx], organic[lidx]);

    if(Zsum_node[nidx] > Lsum + depth[lidx] && !PAST_BOTTOM) {
      Lsum += depth[lidx];
      lidx++;
      if( lidx == Nlayers ) {
	PAST_BOTTOM = TRUE;
	lidx = Nlayers-1;
      }
    }
  }
  return (0);

}

#undef N_INTS

int estimate_layer_ice_content(layer_data_struct *layer,
			       double            *Zsum_node,
			       double            *T,
			       double            *max_moist_node,
#if QUICK_FS
			       double          ***ufwc_table_node,
#else
			       double            *expt_node,
			       double            *bubble_node,
#endif // QUICK_FS
			       double            *depth,
			       double            *max_moist,
#if QUICK_FS
			       double          ***ufwc_table_layer,
#else
			       double            *expt,
			       double            *bubble,
#endif // QUICK_FS
#if SPATIAL_FROST
			       double            *frost_fract,
			       double             frost_slope,
#endif // SPATIAL_FROST
#if EXCESS_ICE
			       double            *porosity,
			       double            *effective_porosity,
#endif // EXCESS_ICE
			       int                Nnodes, 
			       int                Nlayers,
			       char               FS_ACTIVE) {
/**************************************************************
  This subroutine estimates the ice content of all soil 
  moisture layers based on the distribution of soil thermal
  node temperatures.

  layer_struct *layer           structure with all soil moisture layer info
  double       *Zsum_node       soil thermal node cumulative thicknesses (m)
  double       *T               soil thermal node temperatures (C)
  double       *max_moist_node  soil thermal node max moisture content (mm/mm)
  double       *expt_node       soil thermal node exponential ()
  double       *bubble_node     soil thermal node bubbling pressure (cm)
  double       *depth           soil moisture layer thickness (m)
  double       *max_moist       soil layer maximum soil moisture (mm)
  double       *expt            soil layer exponential ()
  double       *bubble          soil layer bubling pressure (cm)
  int           Nnodes          number of soil thermal nodes
  int           Nlayer          number of soil moisture layers

  Modifications:
  11-00 Modified to find ice content in spatial frost bands  KAC
  2007-Apr-24 Zsum removed from declaration.			JCA
  2007-Aug-09 Added features for EXCESS_ICE.			JCA
  2009-Mar-15 Fixed missing else condition for SPATIAL_FROST.	KAC via TJB
  2009-Mar-16 Added resid_moist to argument list, so that
	      min_liq (minimum allowable liquid water content)
	      can be computed here for greater efficiency.	TJB
  2009-May-17 Added options.MIN_LIQ.  Now use of min_liq is
	      optional; if MIN_LIQ is FALSE, min_liq = residual
	      residual moisture and model behaves as before
	      the appearance of min_liq.			TJB
  2009-Jul-31 Removed unused layer_node_fract array.		TJB
  2009-Dec-11 Removed min_liq and options.MIN_LIQ.		TJB

**************************************************************/

  extern option_struct options;

  int    nidx, min_nidx, max_nidx;
  int    lidx, frost_area, Nfrost;
  double Lsum[MAX_LAYERS+1];
#if SPATIAL_FROST
  double tmp_ice[MAX_NODES][FROST_SUBAREAS];
  double tmpT[MAX_NODES][FROST_SUBAREAS+1];
#else
  double tmp_ice[MAX_NODES][1];
  double tmpT[MAX_NODES][1+1];
  double frost_fract[1];
#endif
  double tmpZ[MAX_NODES];
  double min_temp, max_temp, tmp_fract;

#if SPATIAL_FROST
  Nfrost = FROST_SUBAREAS;
#else
  Nfrost = 1;
#endif

  // compute cumulative layer depths
  Lsum[0] = 0;
  for ( lidx = 1; lidx <= Nlayers; lidx++ ) Lsum[lidx] = depth[lidx-1] + Lsum[lidx-1];

  // estimate soil layer average variables
  for ( lidx = 0; lidx < Nlayers; lidx++ ) {

    // Initialize layer variables
    layer[lidx].T = 0.;
#if SPATIAL_FROST
    for ( frost_area = 0; frost_area < FROST_SUBAREAS; frost_area++ ) {
      layer[lidx].ice[frost_area] = 0.;
    }
#else
    layer[lidx].ice = 0.;
#endif

    // Bracket current layer between nodes
    min_nidx = Nnodes-2;
    while( Lsum[lidx] < Zsum_node[min_nidx] && min_nidx > 0 ) min_nidx --;
    max_nidx = 1;
    while( Lsum[lidx+1] > Zsum_node[max_nidx] && max_nidx < Nnodes ) max_nidx ++;
    if ( max_nidx >= Nnodes ) {
      fprintf( stderr, "ERROR: Soil thermal nodes do not extend below bottom soil layer, currently unable to handle this condition.\n" );
      return(ERROR);
    }

    // Get soil node temperatures for current layer
    if ( Zsum_node[min_nidx] < Lsum[lidx] )
      tmpT[min_nidx][Nfrost] = linear_interp(Lsum[lidx], Zsum_node[min_nidx], Zsum_node[min_nidx+1], T[min_nidx], T[min_nidx+1]);
    else tmpT[min_nidx][Nfrost] = T[min_nidx];
    tmpZ[min_nidx] = Lsum[lidx];
    for ( nidx = min_nidx+1; nidx < max_nidx; nidx++ ) {
      tmpT[nidx][Nfrost] = T[nidx];
      tmpZ[nidx] = Zsum_node[nidx];
    }
    if ( Zsum_node[max_nidx] > Lsum[lidx+1] )
      tmpT[max_nidx][Nfrost] = linear_interp(Lsum[lidx+1], Zsum_node[max_nidx-1], Zsum_node[max_nidx], T[max_nidx-1], T[max_nidx]);
    else tmpT[max_nidx][Nfrost] = T[max_nidx];
    tmpZ[max_nidx] = Lsum[lidx+1];

    // distribute temperatures for sub-areas
    for ( nidx = min_nidx; nidx <= max_nidx; nidx++ ) {
#if SPATIAL_FROST
      min_temp = tmpT[nidx][Nfrost] - frost_slope / 2.;
      max_temp = min_temp + frost_slope;
#else
      min_temp = max_temp = tmpT[nidx][Nfrost];
#endif
      for ( frost_area = 0; frost_area < Nfrost; frost_area++ ) {
	if ( Nfrost > 1 ) {
	  if ( frost_area == 0 ) tmp_fract = frost_fract[0] / 2.;
	  else tmp_fract += (frost_fract[frost_area-1] / 2. 
			     + frost_fract[frost_area] / 2.);
	  tmpT[nidx][frost_area] = linear_interp(tmp_fract, 0, 1, min_temp, max_temp);
	}
	else tmpT[nidx][frost_area] = tmpT[nidx][Nfrost];
      }
    }

    // Get soil node ice content for current layer
    if (options.FROZEN_SOIL && FS_ACTIVE) {
      for ( nidx = min_nidx; nidx <= max_nidx; nidx++ ) {
        for ( frost_area = 0; frost_area < Nfrost; frost_area++ ) {
	  tmp_ice[nidx][frost_area] = layer[lidx].moist 
#if QUICK_FS
	    - maximum_unfrozen_water_quick(tmpT[nidx][frost_area], max_moist[lidx], 
					   ufwc_table_layer[lidx]);
#else
#if EXCESS_ICE
            - maximum_unfrozen_water(tmpT[nidx][frost_area], porosity[lidx], 
				     effective_porosity[lidx], max_moist[lidx], 
				     bubble[lidx], expt[lidx]);
#else
	    - maximum_unfrozen_water(tmpT[nidx][frost_area], max_moist[lidx], bubble[lidx], 
				     expt[lidx]);
#endif
#endif
	  if ( tmp_ice[nidx][frost_area] < 0 ) tmp_ice[nidx][frost_area] = 0.;
        }
      }
    }
    else {
      for ( nidx = min_nidx; nidx <= max_nidx; nidx++ ) {
        for ( frost_area = 0; frost_area < Nfrost; frost_area++ ) {
	  tmp_ice[nidx][frost_area] = 0; 
        }
      }
    }

    // Compute average soil layer values
    for ( nidx = min_nidx; nidx < max_nidx; nidx++ ) {
#if SPATIAL_FROST
      for ( frost_area = 0; frost_area < Nfrost; frost_area++ ) {
	layer[lidx].ice[frost_area] += (tmpZ[nidx+1]-tmpZ[nidx])*(tmp_ice[nidx+1][frost_area]+tmp_ice[nidx][frost_area])/2.;
      }
#else
      layer[lidx].ice += (tmpZ[nidx+1]-tmpZ[nidx])*(tmp_ice[nidx+1][0]+tmp_ice[nidx][0])/2.;
#endif  // SPATIAL_FROST
      layer[lidx].T += (tmpZ[nidx+1]-tmpZ[nidx])*(tmpT[nidx+1][Nfrost]+tmpT[nidx][Nfrost])/2.;
    }
#if SPATIAL_FROST
    for ( frost_area = 0; frost_area < Nfrost; frost_area++ )
      layer[lidx].ice[frost_area] /= depth[lidx];
#else
    layer[lidx].ice /= depth[lidx];
#endif  // SPATIAL_FROST
    layer[lidx].T /= depth[lidx];


  }
    
  return (0);

}


int estimate_layer_ice_content_quick_flux(layer_data_struct *layer,
			       double            *depth,
			       double             Dp,
			       double             Tsurf,
			       double             T1,
			       double             Tp,
			       double            *max_moist,
#if QUICK_FS
			       double          ***ufwc_table_layer,
#else
			       double            *expt,
			       double            *bubble,
#endif // QUICK_FS
#if SPATIAL_FROST
			       double            *frost_fract,
			       double             frost_slope,
#endif // SPATIAL_FROST
#if EXCESS_ICE
			       double            *porosity,
			       double            *effective_porosity,
#endif // EXCESS_ICE
			       char               FS_ACTIVE) {
/**************************************************************
  This subroutine estimates the temperature and ice content of all soil 
  moisture layers based on the simplified soil T profile described in
  Liang et al. (1999), and used when QUICK_FLUX is TRUE.

  NOTE: these temperature estimates are much less accurate than those
  of the finite element method (Cherkauer et al. (1999); QUICK_FLUX FALSE).
  Since the Liang et al. (1999) approximation does not account for the
  latent heat fluxes associated with freezing and melting of ice, this
  function should not be called when FROZEN_SOIL is TRUE.

  Modifications:
  2011-Jul-19 Fixed bug in how ice content was computed for case of 
	      SPATIAL_FROST = TRUE.					TJB

********************************************************************/

  extern option_struct options;
  int    lidx, frost_area, Nfrost;
  double Lsum[MAX_LAYERS+1];
  double tmpT, tmp_fract, tmp_ice;
  double min_temp, max_temp;

#if SPATIAL_FROST
  Nfrost = FROST_SUBAREAS;
#else
  Nfrost = 1;
#endif

  // compute cumulative layer depths
  Lsum[0] = 0;
  for ( lidx = 1; lidx <= options.Nlayer; lidx++ ) Lsum[lidx] = depth[lidx-1] + Lsum[lidx-1];

  // estimate soil layer average temperatures
  layer[0].T = 0.5*(Tsurf+T1); // linear profile in topmost layer
  for ( lidx = 1; lidx < options.Nlayer; lidx++ ) {
    layer[lidx].T = Tp - Dp/(depth[lidx])*(T1-Tp)*(exp(-(Lsum[lidx+1]-Lsum[1])/Dp)-exp(-(Lsum[lidx]-Lsum[1])/Dp));
  }

  // estimate soil layer ice contents
  for ( lidx = 0; lidx < options.Nlayer; lidx++ ) {

#if SPATIAL_FROST
    for ( frost_area = 0; frost_area < Nfrost; frost_area++ ) layer[lidx].ice[frost_area] = 0;
#else
    layer[lidx].ice = 0;
#endif

    if (options.FROZEN_SOIL && FS_ACTIVE) {

#if SPATIAL_FROST

      min_temp = layer[lidx].T - frost_slope / 2.;
      max_temp = min_temp + frost_slope;
      for ( frost_area = 0; frost_area < Nfrost; frost_area++ ) {
        if ( frost_area == 0 ) tmp_fract = frost_fract[0] / 2.;
        else tmp_fract += (frost_fract[frost_area-1] / 2. + frost_fract[frost_area] / 2.);
        tmpT = linear_interp(tmp_fract, 0, 1, min_temp, max_temp);
        tmp_ice = layer[lidx].moist
#if QUICK_FS
	    - maximum_unfrozen_water_quick(tmpT, max_moist[lidx], ufwc_table_layer[lidx]);
#else
#if EXCESS_ICE
            - maximum_unfrozen_water(tmpT, porosity[lidx], effective_porosity[lidx], max_moist[lidx], bubble[lidx], expt[lidx]);
#else
	    - maximum_unfrozen_water(tmpT, max_moist[lidx], bubble[lidx], expt[lidx]);
#endif
#endif
        layer[lidx].ice[frost_area] = frost_fract[frost_area] * tmp_ice;
        if (layer[lidx].ice[frost_area] < 0) {
          layer[lidx].ice[frost_area] = 0;
        }
        if (layer[lidx].ice[frost_area] > layer[lidx].moist) {
          layer[lidx].ice[frost_area] = layer[lidx].moist;
        }

      }

#else

      layer[lidx].ice = layer[lidx].moist
#if QUICK_FS
	    - maximum_unfrozen_water_quick(layer[lidx].T, max_moist[lidx], ufwc_table_layer[lidx]);
#else
#if EXCESS_ICE
            - maximum_unfrozen_water(layer[lidx].T, porosity[lidx], effective_porosity[lidx], max_moist[lidx], bubble[lidx], expt[lidx]);
#else
	    - maximum_unfrozen_water(layer[lidx].T, max_moist[lidx], bubble[lidx], expt[lidx]);
#endif
#endif

      if (layer[lidx].ice < 0) {
        layer[lidx].ice = 0;
      }
      if (layer[lidx].ice > layer[lidx].moist) {
        layer[lidx].ice = layer[lidx].moist;
      }

#endif

    }

  }

  return (0);

}


void compute_soil_layer_thermal_properties(layer_data_struct *layer,
					   double            *depth,
					   double            *bulk_dens_min,
					   double            *soil_dens_min,
					   double            *quartz,
					   double            *bulk_density,
					   double            *soil_density,
					   double            *organic,
#if SPATIAL_FROST
					   double            *frost_fract,
#endif
					   int                Nlayers) {
/********************************************************************
  This subroutine computes the thermal conductivity and volumetric
  heat capacity of each soil layer based on its current moisture
  and ice contents.  Ice is only present if the frozen soil
  algorithm is activated.

  layer_data_struct *layer          structure with all soil layer variables
  double            *depth          soil layer depths (m)
  double            *bulk_dens_min  mineral soil bulk density (kg/m^3)
  double            *soil_dens_min  mineral soil particle density (kg/m^3)
  double             quartz         mineral soil quartz content (fract)
  double            *bulk_density   soil bulk density (kg/m^3)
  double            *soil_density   soil particle density (kg/m^3)
  double             organic        soil organic content (fract)
  int                Nlayers        number of soil layers

  MODIFICATIONS:
  2007-Aug-10 Added features for EXCESS_ICE option.			JCA
  2011-Jun-03 Added options.ORGANIC_FRACT.  Soil properties now take
	      organic fraction into account.				TJB
  2011-Jun-10 Added bulk_dens_min and soil_dens_min to arglist of
	      soil_conductivity() to fix bug in commputation of kappa.	TJB

********************************************************************/

  int    lidx;
#if SPATIAL_FROST
  int    frost_area;
#endif
  double moist, ice;

  /* compute layer thermal properties */
  for(lidx=0;lidx<Nlayers;lidx++) {
    moist = layer[lidx].moist / depth[lidx] / 1000;
#if SPATIAL_FROST
    ice = 0;
    for ( frost_area = 0; frost_area < FROST_SUBAREAS; frost_area++ )
      ice += layer[lidx].ice[frost_area] / depth[lidx] / 1000 
	* frost_fract[frost_area];
#else
    ice = layer[lidx].ice / depth[lidx] / 1000;
#endif
    layer[lidx].kappa 
      = soil_conductivity(moist, moist - ice, 
			  soil_dens_min[lidx], bulk_dens_min[lidx], quartz[lidx],
			  soil_density[lidx], bulk_density[lidx], organic[lidx]);
    layer[lidx].Cs 
      = volumetric_heat_capacity(bulk_density[lidx]/soil_density[lidx], moist-ice, ice, organic[lidx]);
  }
}

void find_0_degree_fronts(energy_bal_struct *energy,
			  double            *Zsum_node,
			  double            *T,
			  int                Nnodes) {
/***********************************************************************
  This subroutine reads through the soil thermal nodes and determines 
  the depths of all thawing and freezing fronts that are present.

  energy_bal_struct *energy  energy balance variable structure
  double            *Zsum_node      thermal node depth (m)
  double            *T       thermal node temperatures (C)
  int                Nnodes  number of defined thermal nodes

  Modifications:
  2007-Apr-24 Passing in Zsum_node rather than recalculating.		JCA
  2007-Oct-10 Fixed reference to Zsum_node.				JCA
  2009-Feb-09 Removed dz_node from calculation, now uses Zsum_node for
              consistency.                                              KAC
***********************************************************************/

  int    nidx, fidx; 
  int    Nthaw; /* number of thawing fronts found */
  int    Nfrost; /* number of frost fronts found */
  double tdepth[MAX_FRONTS]; /* thawing frost depths */
  double fdepth[MAX_FRONTS]; /* freezing front depths */
  double deltaz;
  
  /* Initialize parameters */
  Nthaw = Nfrost = 0;
  for(fidx=0;fidx<MAX_FRONTS;fidx++) {
    fdepth[fidx] = MISSING;
    tdepth[fidx] = MISSING;
  }

  /* find 0 degree fronts */
  for(nidx=Nnodes-2;nidx>=0;nidx--) {
    if(T[nidx] > 0 && T[nidx+1] <= 0 && Nthaw<MAX_FRONTS) {
      tdepth[Nthaw] = linear_interp(0,T[nidx],T[nidx+1],Zsum_node[nidx],Zsum_node[nidx+1]);
      Nthaw++;
    }
    else if(T[nidx] < 0 && T[nidx+1] >= 0 && Nfrost<MAX_FRONTS) {
      fdepth[Nfrost] = linear_interp(0,T[nidx],T[nidx+1],Zsum_node[nidx],Zsum_node[nidx+1]);
      Nfrost++;
    }
  }

  /* store thaw depths */
  for(fidx=0;fidx<MAX_FRONTS;fidx++) energy->tdepth[fidx] = tdepth[fidx];
  /* store frost depths */
  for(fidx=0;fidx<MAX_FRONTS;fidx++) energy->fdepth[fidx] = fdepth[fidx];
  energy->Nthaw = Nthaw;
  energy->Nfrost = Nfrost;

}

double maximum_unfrozen_water(double T,
#if EXCESS_ICE
			      double porosity,
			      double effective_porosity,
#endif //EXCESS_ICE
                              double max_moist,
                              double bubble,
                              double expt) {
/**********************************************************************
  This subroutine computes the maximum amount of unfrozen water that
  can exist at the current temperature.

  Modifications:

  2007-Apr-24 Removed T from denominator in equation -
	      according to the proof in Zhang et al. (2007),
	      "Development and Testing of a Frozen Soil Model for
	      the Cold Region Climate Study".				JCA
  2007-Aug-09 Added features for EXCESS_ICE option.			JCA
  2009-Feb-10 Modified to return max_moist if T > 0C.			KAC via TJB
**********************************************************************/

  double unfrozen;

  if ( T <= 0 ) {
    unfrozen = max_moist * pow((-Lf * T) / 273.16 / (9.81 * bubble / 100.), -(2.0 / (expt - 3.0)));
    //INCORRECT: unfrozen = max_moist * pow((-Lf * T) / (T + 273.16) / (9.81 *
    //bubble / 100.), -(2.0 / (expt - 3.0)));
    
    if(unfrozen > max_moist) unfrozen = max_moist;
    if(unfrozen < 0) unfrozen = 0;
    
  }
  else unfrozen = max_moist;

  return (unfrozen);
  
}

#if QUICK_FS
double maximum_unfrozen_water_quick(double   T,
				    double   max_moist,
				    double **table) {
/**********************************************************************
  This subroutine computes the maximum amount of unfrozen water that
  can exist at the current temperature.
**********************************************************************/

  extern double temps[];

  int i;
  double unfrozen;

  i = 1;
  while(T < temps[i] && i < QUICK_FS_TEMPS) i++;
  unfrozen = max_moist * (table[i-1][0] + table[i-1][1] * T);
  if(unfrozen > max_moist) unfrozen = max_moist;
  if(unfrozen < 0) unfrozen = 0;

  return (unfrozen);
}
#endif

layer_data_struct find_average_layer(layer_data_struct *wet,
				     layer_data_struct *dry,
				     double             depth,
				     double             mu) {
/*************************************************************
  This subroutine computes the average soil layer moistures
  between the wet and dry fraction for use in computing 
  energy balance parameters.  Other layer variables are copied 
  from the wet fraction structure since they are they same for 
  wet and dry fractions.

  Modifications:

  2011-Jun-07 Added condition that wet and dry portions are
	      only averaged together if DIST_PRCP is TRUE.	TJB
**************************************************************/

  extern option_struct options;
  layer_data_struct layer;
#if SPATIAL_FROST
  int frost_area;
#endif

  layer = *wet;

  if (options.DIST_PRCP) {

#if SPATIAL_FROST
    for ( frost_area = 0; frost_area < FROST_SUBAREAS; frost_area++ )
      layer.ice[frost_area] = ((wet->ice[frost_area] * mu) 
			       + (dry->ice[frost_area] * (1. - mu)));
#else
    layer.ice = ((wet->ice * mu) + (dry->ice * (1. - mu)));
#endif
    layer.moist = ((wet->moist * mu) + (dry->moist * (1. - mu)));

  }

  return(layer);

}

