#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

#define MAX_CHANGE 0.9

static char vcid[] = "$Id$";

double snow_density(snow_data_struct *snow, double new_snow, double sswq, double Tgrnd, double Tair, double dt)
{
/**********************************************************************
  snow_density		Keith Cherkauer		May 28, 1997

  This subroutine computes the snow density based on swe and snow metamorphism.

  One of two algorithms may be used, depending on the value of options.SNOW_DENSITY:

    DENS_SNTHRM = Algorithm is taken from SNTHERM89, adjusted for an essentially
		  single-layer model.

      UNITS:
	new_snow	mm	new snow
	sswq		m	initial snow pack snow water equivalent (before adding new snowfall)
	Tgrnd		C	ground temperature
	Tair		C	air temperature
	dt		h	time step length
	density		kg/m^3	snow density



    DENS_BRAS   = Original algorithm, originally based on a plot of seasonal
		  variation of typical snow densities found in Bras (figure 6.10,
		  p 258).  Because this equation was developed by regression
		  against data from Southern Manitoba, it therefore is limited
		  in applicability to other regions.

      UNITS:
	new_snow	         mm	new snow
	air_temp	         C	current air temperature
	swq		m	snow water equivalent	
	depth		m	snow pack depth
	density            kg/m^3   snow density



  Modified:
  08-19-99 Added check to make sure that the change in snowpack depth
           due to new snow does not exceed the actual depth of the 
	   pack.							Bart
  06-30-03 Added check to keep compression from aging from exceeding
           the actual depth of the snowpack.				KAC
  08-Oct-03 Modified the checks on delta_depth (mentioned above)
	    so that the condition is delta_depth > MAX_CHANGE*depth.	TJB
  08-Oct-03 Modified compression due to aging to only be calculated
	    if depth > 0.						TJB
  2008-Feb-17 Moved parameters related to snow densification to
	      snow.h.							TJB
  2008-Feb-17 Replaced previous algorithm with one based on SNTHERM89
	      adjusted for an essentially single-layer model.		KMA via TJB
  2008-Apr-21 Added option to use previous algorithm for backwards
	      compatibility.  DENS_BRAS = original algorithm; DENS_SNTHRM
	      = SNTHERM89.						TJB

**********************************************************************/

  extern option_struct   options;
  double density_new;
  double density;
  double depth;
  double swq;
  double CR; /* compaction rate */
  double dexpf;
  double ddz1; /* rate of settling of snowpack due to destructive metamorphism */
  double ddz2; /* rate of compaction of snowpack due to overburden */
  double Tavg; /* average snowpack temperature */
  double c3,c4; /* snow densification factors */
  double dm; /* upper snow density limit for the settling process */
  double Ps; /* snow load pressure (N/m^2) */
  double f; /* effective compaction coefficient */
  double overburden;
  double viscosity;
  double delta_depth;
  double depth_new;

  density = 0.0;

  /* Estimate density of new snow based on air temperature */
  if (new_snow > 0.)
    density_new = new_snow_density(Tair);
  else
    density_new = 0.0;

  /* Estimate average snowpack temperature */
//  Tavg = (snow->surf_temp+Tgrnd)/2.0+KELVIN;
  Tavg = snow->surf_temp+KELVIN;

  if (options.SNOW_DENSITY == DENS_SNTHRM) {

    if (new_snow > 0.) {
      if (snow->depth>0.0)
        density = snow->density;
      else
        density = density_new;
    }
    else {
      density = snow->density;
    }

    dexpf = exp(-SNDENS_C1*(KELVIN-Tavg));

    /* Settling due to destructive metamorphism */
    if (new_snow>0.0 && density_new>0.0)
      dm = (SNDENS_DMLIMIT > 1.15*density_new) ? SNDENS_DMLIMIT : 1.15*density_new;
    else
      dm = SNDENS_DMLIMIT;

    if (density <= dm) {
      c3 = 1.0;
      c4 = 1.0;
    }
    else {
      c3 = exp(-0.046*(density-dm));
      c4 = 1.0;
    }
    if ((snow->surf_water+snow->pack_water)/snow->depth > 0.01)
      c4 = 2.0; /* presence of wet snow */

    ddz1 = -SNDENS_C2 * c3* c4 * dexpf;

    /* Compaction due to overburden */
    // swq in this context is the amount of snow whose weight contributes to compaction
//    f = sswq/80. * exp(-sswq/100.);
    f = SNDENS_F;
    swq = new_snow/1000. + f*sswq; /* Currently VIC essentially has only one layer of snow, so compaction due to overburden will come mostly from new snowfall. */

    if (new_snow > 0.0) {
      Ps = 0.5*G*RHO_W*swq;
      ddz2 = -Ps / SNDENS_ETA0 * exp(-(-SNDENS_C5*(Tavg-KELVIN) + SNDENS_C6*density));
    }
    else {
      ddz2 = 0.0;
    }

    /* Calculate compaction rate and new snow density */
    CR = -ddz1-ddz2;
    density = density*(1+CR*dt*SECPHOUR);

  }

  else if (options.SNOW_DENSITY == DENS_BRAS) {

    depth = snow->depth;
    swq = sswq;

    /** Compaction of snow pack by new snow fall **/
    /** Bras pg. 257 **/

    if ( new_snow > 0 ) {

      if ( depth > 0. ) {

        /* Compact current snowpack by weight of new snowfall */
        delta_depth = ( ((new_snow / 25.4) * (depth / 0.0254)) / (swq / 0.0254)
                      * pow( (depth / 0.0254) / 10., 0.35) ) * 0.0254;
        if (delta_depth > MAX_CHANGE * depth) delta_depth = MAX_CHANGE * depth;
        depth_new = new_snow / density_new;
        depth = depth - delta_depth + depth_new;
        swq += new_snow / 1000.;
        density = 1000. * swq / depth;

      }
      else {

        /* no snowpack present, so snow density equals that of new snow */
        density  = density_new;
        swq     += new_snow / 1000.;
        depth    = 1000. * swq / density;

      }

    }
    else density = 1000. * swq / snow->depth;

    /** Densification of the snow pack due to aging **/
    /** based on SNTHRM89 R. Jordan 1991 - used in Bart's DHSVM code **/
    if ( depth > 0. ) {

      overburden  = 0.5 * G * RHO_W * swq;
      viscosity   = SNDENS_ETA0 * exp(-SNDENS_C5*(Tavg-KELVIN) + SNDENS_C6*density);
      delta_depth = overburden / viscosity * depth * dt * SECPHOUR;
      if (delta_depth > MAX_CHANGE * depth) delta_depth = MAX_CHANGE * depth;
      depth      -= delta_depth;
      density     = 1000. * swq / depth;

    }

  }

  return (density);

}

double new_snow_density(double air_temp) {
/**********************************************************************
  new_snow_density		Keith Cherkauer		May 28, 1997

  This routine estimates the density of new snow.

  Modified:
  2008-Feb-17 Modified to use the algorithm of Lundberg and
	      Pomeroy (1998).						KMA via TJB
  2008-Apr-21 Added option to use previous algorithm for backwards
	      compatibility.  DENS_BRAS = original algorithm; DENS_SNTHRM
	      = Lundberg and Pomeroy (1998).				TJB
**********************************************************************/
  extern option_struct   options;
  double density_new;

  density_new = 0.0;

  if (options.SNOW_DENSITY == DENS_SNTHRM) {
    density_new = 67.9 + 51.3 * exp(air_temp/2.6);
  }
  else if (options.SNOW_DENSITY == DENS_BRAS) {
    air_temp = air_temp * 9. / 5. + 32.;
    if(air_temp > 0) density_new = (double)NEW_SNOW_DENSITY + 1000.
		     * (air_temp / 100.) * (air_temp / 100.);
    else density_new = (double)NEW_SNOW_DENSITY;
  }

  return (density_new);

}


double snow_albedo(double new_snow,
                   double swq,
                   double depth,
                   double albedo,
                   double cold_content,
                   double dt,
                   int    last_snow,
		   char   MELTING) {
/**********************************************************************
  snow_albedo		Keith Cherkauer		June 10, 1997

  This subroutine computes the snow pack surface albedo.  Original version
  computed albedo as a function of snow age and season, based on the algorithm
  of the US Army Corps of Engineers.  More recently, added the option to
  use the algorithm of Sun et al., 1999, which depends only on snow age and
  cold content (independent of time of year).

  Modified:
  06-15-02 Added MELTING flag which tells the algorithm whether or not
           the pack was melting previously.  This locks the albedo
           onto the lower ablation albedo curve until a sufficiently
           large new snow event resets the surface albedo.       KAC
  2008-Apr-21 Added flag to include new albedo calculation based on
	      albedo algorithm described in Sun, S., J. Jin, and Y. Xue,
	      A simple snow-atmosphere-soil transfer model, J. Geophys. Res.,
	      104 (D16), 19,587-19,597, 1999.				KAC via TJB
**********************************************************************/

  extern option_struct   options;

  /** New Snow **/
  if(new_snow > TraceSnow  && cold_content < 0.0 ) albedo = NEW_SNOW_ALB;

  /** Aged Snow **/
  else if(swq > 0.0) {

    /* Accumulation season */
    if ( cold_content < 0.0 && !MELTING )
      albedo = NEW_SNOW_ALB*pow(SNOW_ALB_ACCUM_A, 
			  pow((double)last_snow * dt / 24.,
			      SNOW_ALB_ACCUM_B));

    /* Melt Season */
    else
      albedo = NEW_SNOW_ALB*pow(SNOW_ALB_THAW_A, 
			  pow((double)last_snow * dt / 24.,
			      SNOW_ALB_THAW_B));
  }

  else
    /* No snow falling or present */
    albedo = 0;

  return(albedo);

}
