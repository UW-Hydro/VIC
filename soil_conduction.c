#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

static char vcid[] = "$Id: soil_conduction.c,v 4.1.2.3 2004/05/11 19:46:33 tbohn Exp $";

double soil_conductivity(double moist, 
			 double Wu, 
			 double soil_density, 
			 double bulk_density,
			 double quartz) {
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
  double soil_density  soil density (kg m-3)
  double bulk_density  soil bulk density (kg m-3)
  double quartz        soil quartz content (fraction)

**********************************************************************/
  double Ke;
  double Ki = 2.2;	/* thermal conductivity of ice (W/mK) */
  double Kw = 0.57;	/* thermal conductivity of water (W/mK) */
  double Ksat;
  double Ks;		/* thermal conductivity of solid (W/mK)
			   function of quartz content */
  double Kdry;
  double Sr;		/* fractional degree of saturation */
  double K;
  double porosity;

  Kdry = (0.135*bulk_density+64.7)/(soil_density-0.947*bulk_density);

  if(moist>0.) {

    porosity = 1.0 - bulk_density / soil_density;

    Sr = moist/porosity;

    Ks = pow(7.7,quartz) * pow(2.2,1.0-quartz);

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


#define organic_fract 0.00

double volumetric_heat_capacity(double soil_fract,
                                double water_fract,
                                double ice_fract) {
/**********************************************************************
  This subroutine calculates the soil volumetric heat capacity based 
  on the fractional volume of its component parts.

  Constant values are volumetric heat capacities in J/m^3/K
	Soil value is for clay or quartz - assumed for all other types

  double soil_fract   fraction of soil volume composed of actual soil (fract)
  double water_fract  fraction of soil volume composed of liquid water (fract)
  double ice_fract    fraction of soil volume composed of ice (fract)

**********************************************************************/

  double Cs;

  Cs = 2.0e6 * (soil_fract - organic_fract);
  Cs += 4.2e6 * water_fract;
  Cs += 1.9e6 * ice_fract;
  Cs += 2.7e6 * organic_fract;
  Cs += 1.3e3 * (1. - (soil_fract + water_fract + ice_fract + organic_fract));

  return (Cs);

}

#undef organic_fract

void set_node_parameters(double   *dz_node,
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
			 float   **layer_node_fract,
#if QUICK_FS
			 double ***ufwc_table_node,
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
  int       Nnodes           number of soil thermal nodes
  int       Nlayers          number of soil moisture layers
  char      FS_ACTIVE        TRUE if frozen soils are active in grid cell

  Modifications:

  02-11-00 Modified to remove node zone averages, node parameters
           are now set based on the parameter value of the layer 
           in which they fall.  Averaging of layer properties 
	   only occurs if the node falls directly on a layer
	   boundary.                                        KAC
  11-May-04 (Port from 4.1.0) Modified to correct differences between
	    calculations to determine maximum node moisture and node
	    moisture, so that nodes on the boundary between soil
	    layers are computed the same way for both.		TJB

**********************************************************************/

  extern option_struct options;
#if QUICK_FS
  extern double temps[];
#endif

  char   PAST_BOTTOM;
  int    nidx, lidx;
  int    tmplidx;
  double Lsum; /* cumulative depth of moisture layer */
  double Zsum; /* cumulative depth of thermal node */
  double deltaL[MAX_LAYERS+1];
  double fract;
  double tmpfract;
#if QUICK_FS
  int    ii;
  double Aufwc;
  double Bufwc;
#endif

  PAST_BOTTOM = FALSE;
  lidx = 0;
  Lsum = 0.;
  Zsum = 0.;

  /* set node parameters */
  for(nidx=0;nidx<Nnodes;nidx++) {

    if(Zsum == Lsum + depth[lidx] && nidx != 0 && lidx != Nlayers-1) {
      /* node on layer boundary */
      max_moist_node[nidx] = (max_moist[lidx] / depth[lidx]
                              + max_moist[lidx+1] / depth[lidx+1]) / 1000 / 2.;
      expt_node[nidx]      = (expt[lidx] + expt[lidx+1]) / 2.;
      bubble_node[nidx]    = (bubble[lidx] + bubble[lidx+1]) / 2.;
    }
    else { 
      /* node completely in layer */
      max_moist_node[nidx] = max_moist[lidx] / depth[lidx] / 1000;
      expt_node[nidx]      = expt[lidx];
      bubble_node[nidx]    = bubble[lidx];
    }      
    Zsum += (dz_node[nidx] + dz_node[nidx+1]) / 2.;
    if(Zsum > Lsum + depth[lidx] && !PAST_BOTTOM) {
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
    alpha[nidx] = ((dz_node[nidx+2] + dz_node[nidx+1]) / 2.0 
		   + (dz_node[nidx+1] + dz_node[nidx]) / 2.0);
    beta[nidx] = ((dz_node[nidx+2] + dz_node[nidx+1]) 
		  * (dz_node[nidx+2] + dz_node[nidx+1])) / 4.0 
      + ((dz_node[nidx+1]+dz_node[nidx]) 
	 * (dz_node[nidx+1]+dz_node[nidx])) / 4.0;
    gamma[nidx] = ((dz_node[nidx+2] + dz_node[nidx+1]) / 2.0 
		   - (dz_node[nidx+1] + dz_node[nidx]) / 2.0);
  }
  if(options.NOFLUX) {
    /* no flux bottom boundary activated */
    alpha[Nnodes-2] = ((dz_node[Nnodes-1] + dz_node[Nnodes-1]) / 2.0 
		       + (dz_node[Nnodes-1] + dz_node[Nnodes-2]) / 2.0);
    beta[Nnodes-2] = ((dz_node[Nnodes-1] + dz_node[Nnodes-1]) 
		      * (dz_node[Nnodes-1] + dz_node[Nnodes-1]))/4.0 
      + ((dz_node[Nnodes-1]+dz_node[Nnodes-2])
	 * (dz_node[Nnodes-1]+dz_node[Nnodes-2])) / 4.0;
    gamma[Nnodes-2] = ((dz_node[Nnodes-1] + dz_node[Nnodes-1]) / 2.0 
		       - (dz_node[Nnodes-1] + dz_node[Nnodes-2]) / 2.0);
  }

  /* set fraction of soil thermal node in each soil layer */
  Lsum = 0;
  deltaL[Nlayers] = 0;
  for(lidx=0;lidx<Nlayers;lidx++) {
    deltaL[lidx] = depth[lidx];
    deltaL[Nlayers] -= depth[lidx];
  }
  for(nidx=0;nidx<Nnodes-1;nidx++) 
    deltaL[Nlayers] += (dz_node[nidx] + dz_node[nidx+1]) / 2.;
  for(lidx=0;lidx<=Nlayers;lidx++) {
    Zsum = -dz_node[0] / 2.;
    for(nidx=0;nidx<Nnodes;nidx++) {
      if( Zsum < Lsum && Zsum + dz_node[nidx] >= Lsum ) {
	layer_node_fract[lidx][nidx] = 1. 
	  - (float)linear_interp(Lsum, Zsum, Zsum + dz_node[nidx], 0, 1);
	if(Lsum + deltaL[lidx] < Zsum + dz_node[nidx])
	  layer_node_fract[lidx][nidx] -= 
	    (float)linear_interp(Lsum + deltaL[lidx], Zsum, Zsum + dz_node[nidx], 1, 0);
      }
      else if( Zsum < Lsum + deltaL[lidx] && 
	       Zsum + dz_node[nidx] >= Lsum + deltaL[lidx] )
	layer_node_fract[lidx][nidx] = 
	  (float)linear_interp(Lsum + deltaL[lidx], Zsum, 
			       Zsum + dz_node[nidx], 0, 1);
      else if( Zsum >= Lsum && Zsum + dz_node[nidx] <= Lsum + deltaL[lidx] )
	layer_node_fract[lidx][nidx] = 1;
      else layer_node_fract[lidx][nidx] = 0;
      Zsum += dz_node[nidx];
    }
    Lsum += deltaL[lidx];
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

void distribute_node_moisture_properties(double *moist_node,
					 double *ice_node,
					 double *kappa_node,
					 double *Cs_node,
					 double *dz_node,
					 double *T_node,
					 double *max_moist_node,
#if QUICK_FS
					 double ***ufwc_table_node,
#else
					 double *expt_node,
					 double *bubble_node,
#endif
					 double *moist,
					 double *depth,
					 double *soil_density,
					 double *bulk_density,
					 double *quartz,
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
  double *dz_node         thermal node thickness (m)
  double *T_node          thermal node temperature (C)
  double *max_moist_node  thermal node maximum moisture content (mm/mm)
  double *expt_node       thermal node exponential
  double *bubble_node     thermal node bubbling pressure (cm)
  double *moist           soil layer moisture (mm)
  double *depth           soil layer thickness (m)
  double  soil_density    soil density (kg m-3)
  double *bulk_density    soil layer bulk density (kg m-3)
  double  quartz          soil quartz content (fract)
  int     Nnodes          number of soil thermal nodes
  int     Nlayers         number of soil moisture layers

  Modifications:

  02-11-00 Modified to remove node zone averages, node parameters
           are now set based on the parameter value of the layer 
           in which they fall.  Averaging of layer properties 
	   only occurs if the node falls directly on a layer
	   boundary.                                        KAC
  11-May-04 (Port from 4.1.0) Modified to check that node soil
	    moisture is less than or equal to maximum node soil
	    moisture, otherwise an error is printed to the screen
	    and the model exits.				TJB

*********************************************************************/

  extern option_struct options;

  char ErrStr[MAXSTRING];
  char PAST_BOTTOM;
  int nidx, lidx, i;
  int tmplidx;
  double Lsum; /* cumulative depth of moisture layer */
  double Zsum; /* cumulative depth of thermal node */
  double tmp_moist;
  double tmp_T;
  double fract;
  double tmpfract;
  double ice_sum;
  double dz, dz_int;
  double T_upper;
  double T_mid;
  double T_lower;
  double T_int;
  double factor;
  double soil_fract;

  lidx = 0;
  Lsum = 0.;
  Zsum = 0.;
  PAST_BOTTOM = FALSE;

  /* node estimates */
  for(nidx=0;nidx<Nnodes;nidx++) {
 
    if(Zsum == Lsum + depth[lidx] && nidx != 0 && lidx != Nlayers-1) {
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
    if (abs(moist_node[nidx]-max_moist_node[nidx]) > SMALL) {
      sprintf( ErrStr, "Node soil moisture, %f, exceeds maximum node soil moisuttre, %f.", moist_node[nidx], max_moist_node[nidx] );
      vicerror(ErrStr);
    }

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
					     max_moist_node[nidx], 
					     bubble_node[nidx],
					     expt_node[nidx]);
#endif
      if(ice_node[nidx]<0) ice_node[nidx]=0;

      /* compute thermal conductivity */
      kappa_node[nidx] 
	= soil_conductivity(moist_node[nidx], moist_node[nidx] 
			    - ice_node[nidx], soil_density[lidx],
			    bulk_density[lidx], quartz[lidx]);

    }
    else {
      /* compute moisture and ice contents */
      ice_node[nidx]   = 0;
      /* compute thermal conductivity */
      kappa_node[nidx] 
	= soil_conductivity(moist_node[nidx], moist_node[nidx], 
			    soil_density[lidx], bulk_density[lidx], 
			    quartz[lidx]);
    }
    /* compute volumetric heat capacity */
    Cs_node[nidx] = volumetric_heat_capacity(bulk_density[lidx] 
					     / soil_density[lidx],
					     moist_node[nidx] - ice_node[nidx],
					     ice_node[nidx]);

    Zsum += (dz_node[nidx] + dz_node[nidx+1]) / 2.;
    if(Zsum > Lsum + depth[lidx] && !PAST_BOTTOM) {
      Lsum += depth[lidx];
      lidx++;
      if( lidx == Nlayers ) {
	PAST_BOTTOM = TRUE;
	lidx = Nlayers-1;
      }
    }
  }

}

#undef N_INTS

void estimate_layer_ice_content(layer_data_struct *layer,
				double            *dz,
				double            *T,
				double            *max_moist_node,
#if QUICK_FS
				double          ***ufwc_table_node,
#else
				double            *expt_node,
				double            *bubble_node,
#endif
				double            *depth,
				double            *max_moist,
#if QUICK_FS
				double          ***ufwc_table_layer,
#else
				double            *expt,
				double            *bubble,
#endif
				double            *bulk_density,
				double            *soil_density,
				double            *quartz,
				float            **layer_node_fract,
				int                Nnodes, 
				int                Nlayers,
				char               FS_ACTIVE) {
/**************************************************************
  This subroutine estimates the ice content of all soil 
  moisture layers based on the distribution of soil thermal
  node temperatures.

  layer_struct *layer           structure with all soil moisture layer info
  double       *dz              soil thermal node thicknesses (m)
  double       *T               soil thermal node temperatures (C)
  double       *max_moist_node  soil thermal node max moisture content (mm/mm)
  double       *expt_node       soil thermal node exponential ()
  double       *bubble_node     soil thermal node bubbling pressure (cm)
  double       *depth           soil moisture layer thickness (m)
  double       *max_moist       soil layer maximum soil moisture (mm)
  double       *expt            soil layer exponential ()
  double       *bubble          soil layer bubling pressure (cm)
  double       *bulk_density    soil layer bulk density (kg m-3)
  double        soil_density    soil layer soil density (kg m-3)
  double        quartz          soil layer quartz content (fract)
  int           Nnodes          number of soil thermal nodes
  int           Nlayer          number of soil moisture layers

**************************************************************/

  extern option_struct options;

  int    nidx;
  int    lidx;
  double Lsum; /* cumulative depth of moisture layer */
  double Zsum; /* cumulative depth of thermal node */
  double deltaz;
  double fract; /* fract of internodal region in layer */
  double boundT; /* soil temperature between layers */
  double tmp_ice;
  double ice_content[2]; /* stores estimated nodal ice content */
  double kappa_layer[2]; /* stores node thermal conductivity */
  double Cs_layer[2]; /* stores node volumetric heat capacity */
  double T_sum; /* summation of nodal temperatures within the layer */
  double ice_sum; /* summation of nodal ice contents within the layer */
  double dz_sum; /* summation of depth used in computing statistics */

  for(lidx=0;lidx<Nlayers;lidx++) {
    layer[lidx].T = 0.;
    layer[lidx].ice = 0.;
    for(nidx=0;nidx<Nnodes;nidx++) {
      if(layer_node_fract[lidx][nidx] > 0) {
	layer[lidx].T += T[nidx] * layer_node_fract[lidx][nidx] * dz[nidx];
	if(T[nidx] < 0 && options.FROZEN_SOIL && FS_ACTIVE) {
	  tmp_ice = layer[lidx].moist 
#if QUICK_FS
		     - maximum_unfrozen_water_quick(T[nidx], 
						    max_moist[lidx], 
						    ufwc_table_layer[lidx]);
#else
	- maximum_unfrozen_water(T[nidx], max_moist[lidx], bubble[lidx], 
				 expt[lidx]);
#endif
	  if(tmp_ice < 0) tmp_ice = 0;
	}
	else tmp_ice = 0;
	layer[lidx].ice += tmp_ice * layer_node_fract[lidx][nidx] * dz[nidx];
      }
    }
    layer[lidx].T /= depth[lidx];
    layer[lidx].ice /= depth[lidx];
  }

}

void compute_soil_layer_thermal_properties(layer_data_struct *layer,
					   double            *depth,
					   double            *bulk_density,
					   double            *soil_density,
					   double            *quartz,
					   int                Nlayers) {
/********************************************************************
  This subroutine computes the thermal conductivity and volumetric
  heat capacity of each soil layer based on its current moisture
  and ice contents.  Ice is only present if the frozen soil
  algorithm is activated.

  layer_data_struct *layer          structure with all soil layer variables
  double            *depth          soil layer depths (m)
  double            *bulk_density   soil layer bulk density (kg/m^3)
  double             soil_density   soil layer soil density (kg/m^3)
  double             quartz         soil layer quartz content (fract)
  int                Nlayers        number of soil layers

********************************************************************/

  int lidx;
  double moist, ice;

  /* compute layer thermal properties */
  for(lidx=0;lidx<Nlayers;lidx++) {
    moist = layer[lidx].moist / depth[lidx] / 1000;
    ice = layer[lidx].ice / depth[lidx] / 1000;
    layer[lidx].kappa 
      = soil_conductivity(moist, moist - ice, soil_density[lidx], 
			  bulk_density[lidx], quartz[lidx]);
    layer[lidx].Cs 
      = volumetric_heat_capacity(bulk_density[lidx] / soil_density[lidx], 
				 moist - ice, ice);
  }
}

void find_0_degree_fronts(energy_bal_struct *energy,
			  double            *dz,
			  double            *T,
			  int                Nnodes) {
/***********************************************************************
  This subroutine reads through the soil thermal nodes and determines 
  the depths of all thawing and freezing fronts that are present.

  energy_bal_struct *energy  energy balance variable structure
  double            *dz      thermal node thicknesses (m)
  double            *T       thermal node temperatures (C)
  int                Nnodes  number of defined thermal nodes

***********************************************************************/

  int    nidx, fidx; 
  int    Nthaw; /* number of thawing fronts found */
  int    Nfrost; /* number of frost fronts found */
  double tdepth[MAX_FRONTS]; /* thawing frost depths */
  double fdepth[MAX_FRONTS]; /* freezing front depths */
  double Zsum;
  double deltaz;
  
  /* Initialize parameters */
  Zsum = 0;
  for(nidx=0;nidx<Nnodes-1;nidx++)
    Zsum += (dz[nidx] + dz[nidx+1]) / 2.;
  Nthaw = Nfrost = 0;
  for(fidx=0;fidx<MAX_FRONTS;fidx++) {
    fdepth[fidx] = MISSING;
    tdepth[fidx] = MISSING;
  }

  /* find 0 degree fronts */
  for(nidx=Nnodes-2;nidx>=0;nidx--) {
    deltaz = (dz[nidx] + dz[nidx+1]) / 2.;
    if(T[nidx] > 0 && T[nidx+1] <= 0 && Nthaw<MAX_FRONTS) {
      tdepth[Nthaw] = linear_interp(0,T[nidx],T[nidx+1],Zsum-deltaz,Zsum);
      Nthaw++;
    }
    else if(T[nidx] < 0 && T[nidx+1] >= 0 && Nfrost<MAX_FRONTS) {
      fdepth[Nfrost] = linear_interp(0,T[nidx],T[nidx+1],Zsum-deltaz,Zsum);
      Nfrost++;
    }
    Zsum -= deltaz;
  }

  /* store thaw depths */
  for(fidx=0;fidx<MAX_FRONTS;fidx++) energy->tdepth[fidx] = tdepth[fidx];
  /* store frost depths */
  for(fidx=0;fidx<MAX_FRONTS;fidx++) energy->fdepth[fidx] = fdepth[fidx];
  energy->Nthaw = Nthaw;
  energy->Nfrost = Nfrost;

}

double maximum_unfrozen_water(double T,
                              double max_moist,
                              double bubble,
                              double expt) {
/**********************************************************************
  This subroutine computes the maximum amount of unfrozen water that
  can exist at the current temperature.
**********************************************************************/

  double unfrozen;

  unfrozen = max_moist * pow((-Lf * T) / (T
      + 273.16) / (9.18 * bubble / 100.), -(2.0 / (expt - 3.0)));
  if(unfrozen > max_moist) unfrozen = max_moist;
  if(unfrozen < 0) unfrozen = 0;

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
**************************************************************/

  layer_data_struct layer;

  layer = *wet;

  layer.ice = ((wet->ice * mu) + (dry->ice * (1. - mu)));
  layer.moist = ((wet->moist * mu) + (dry->moist * (1. - mu)));

  return(layer);

}

