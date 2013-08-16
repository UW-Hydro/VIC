#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

static char vcid[] = "$Id: snow_utility.c,v 4.2.2.4 2004/05/06 19:57:35 tbohn Exp $";

#define ETA0       (3.6e6)  /* viscosity of snow at T = 0C and density = 0
			       used in calculation of true viscosity (Ns/m2) */

#define G          9.81     /* gravitational accelleration (m/(s^2)) */
#define C5         0.08     /* constant used in snow viscosity calculation,
			       taken from SNTHRM.89 (/C)*/
#define C6         0.021    /* constant used in snow viscosity calculation,
			       taken from SNTHRM.89 (kg/m3) */
#define MAX_CHANGE 0.9      /* maximum fraction of snowpack depth change
			       caused by new snow */

double snow_density(int date,
                    double new_snow,
                    double air_temp,
                    double swq,
                    double depth, 
                    double coldcontent,
                    double dt,
		    double Tsurf) {
/**********************************************************************
  snow_density		Keith Cherkauer		May 28, 1997

  This subroutine computes the snow density based on the day of the 
  year.  Density information comes from a plot of seasonal variation
  of typical snow densities found in Bras (Figure 6.10, p 258).  The
  equation was developed by regressing against the curve for Southern 
  Manitoba, so this routine should be modified if used outside the 
  plains of south central Canada, and the north central US.

  UNITS:
	new_snow	         mm	new snow
	air_temp	         C	current air temperature
	swq		m	snow water equivalent	
	depth		m	snow pack depth
	density            kg/m^3   snow density

  Modified:
  08-19-99 Added check to make sure that the change in snowpack depth
           due to new snow does not exceed the actual depth of the 
	   pack.                                               Bart
  06-30-03 Added check to keep compression from aging from exceeding
           the actual depth of the snowpack.                   KAC
  08-Oct-03 Modified the checks on delta_depth (mentioned above)
	    so that the condition is delta_depth > MAX_CHANGE*depth. TJB
  08-Oct-03 Modified compression due to aging to only be calculated
	    if depth > 0.					TJB

**********************************************************************/

  double density;
  double delta_depth;
  double density_new;
  double depth_new;
  double overburden;
  double viscosity;
  double deltadepth;

  /** Compaction of snow pack by new snow fall **/
  /** Bras pg. 257 **/

  if( new_snow > 0 ) {

    /* Estimate density of new snow based on air temperature */

    density_new = new_snow_density(air_temp);

    if( depth>0. ) {

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

      density = density_new;

      swq     += new_snow / 1000.;

      depth    = 1000. * swq / density;

    }

  }
  else density = 1000. * swq / depth;

  /** Densification of the snow pack due to aging **/
  /** based on SNTHRM89 R. Jordan 1991 - used in Bart's DHSVM code **/

  if ( depth > 0. ) {

    overburden  = 0.5 * G * RHO_W * swq;

    viscosity   = ETA0 * exp(-C5 * Tsurf + C6 * density);

    delta_depth = overburden / viscosity * depth * dt * SECPHOUR;

    if (delta_depth > MAX_CHANGE * depth) delta_depth = MAX_CHANGE * depth;
      
    depth      -= delta_depth;

    density     = 1000. * swq / depth;

  }

  return (density);

}

double new_snow_density(double air_temp) {
  /**************************************************
    This routine estimates the density of new snow.
  **************************************************/
  double density_new;

  air_temp = air_temp * 9. / 5. + 32.;
  if(air_temp > 0) density_new = (double)NEW_SNOW_DENSITY + 1000.
		     * (air_temp / 100.) * (air_temp / 100.);
  else density_new = (double)NEW_SNOW_DENSITY;
  return (density_new);
}

double snow_albedo(double new_snow,
                   double swq,
                   double cold_content,
                   double dt,
                   int last_snow,
		   char MELTING) {
/**********************************************************************
  snow_albedo		Keith Cherkauer		June 10, 1997

  This subroutine computes the snow pack surface albedo based on snow
  age and season, using the tables generated in snow_table_albedo.

  Modified:
  06-15-02 Added MELTING flag which tells the algorithm whether or not
           the pack was melting previously.  This locks the albedo
           onto the lower ablation albedo curve until a sufficiently
           large new snow event resets the surface albedo.       KAC
**********************************************************************/

  double albedo;

  /** New Snow **/
  if(new_snow > 0.0) albedo = NEW_SNOW_ALB;

  /** Aged Snow **/
  else if(swq > 0.0) {

    /* Accumulation season */
    if(cold_content < 0.0 && !MELTING )
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

#undef ETA0
#undef G
#undef C5     
#undef C6
