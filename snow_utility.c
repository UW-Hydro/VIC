#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

static char vcid[] = "$Id$";

#define ETA0           (3.6e6)  /* viscosity of snow at T = 0C and density = 0
                                   used in calculation of true viscosity 
                                   (Ns/m2) */

#define G              9.81     /* gravitational accelleration (m/(s^2)) */
#define C5         .08          /* constant used in snow viscosity calculation,
                                   taken from SNTHRM.89 (/C)*/
#define C6         .021         /* constant used in snow viscosity calculation,
                                   taken from SNTHRM.89 (kg/m3) */

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

  if(new_snow > 0) {

    air_temp = air_temp * 9. / 5. + 32.;
    if(air_temp > 0) density_new = (double)NEW_SNOW_DENSITY + 1000.
                                 * (air_temp / 100.) * (air_temp / 100.);
    else density_new = (double)NEW_SNOW_DENSITY;

    if(depth>0.) {
      delta_depth = ( ((new_snow / 25.4) * (depth / 0.0254)) / (swq / 0.0254)
                    * pow( (depth / 0.0254) / 10., 0.35) ) * 0.0254;
  
      depth_new = new_snow / density_new;
  
      depth = depth - delta_depth + depth_new;

      swq += new_snow / 1000.;

      density = 1000. * swq / depth;

    }
    else {

      density = density_new;

      swq += new_snow / 1000.;

    }

  }
  else if(coldcontent < 0) {

    density = 1000. * swq / depth;

    density += 1.2 * dt / 24.;

  }
  else density = 1000. * swq / depth;

  /** Densification of the snow pack due to aging **/
  /** based on SNTHRM89 R. Jordan 1991 - used in Bart's DHSVM code **/

  depth = 1000. * swq / density;
  overburden = 0.5 * G * RHO_W * swq;
  viscosity = ETA0 * exp(-C5 * Tsurf + C6 * density);
  deltadepth = -overburden/viscosity*depth*dt*SECPHOUR;
  depth += deltadepth;
  density = 1000. * swq / depth;

  return (density);

}

double snow_albedo(double new_snow,
                   double swq,
                   double cold_content,
                   double dt,
                   int last_snow) {
/**********************************************************************
  snow_albedo		Keith Cherkauer		June 10, 1997

  This subroutine computes the snow pack surface albedo based on snow
  age and season, using the tables generated in snow_table_albedo.
**********************************************************************/

  double albedo;

  /** New Snow **/
  if(new_snow > 0.0) albedo = NEW_SNOW_ALB;

  /** Aged Snow **/
  else if(swq>0.0) {

    /* Accumulation season */
    if(cold_content < 0.0)
      albedo = NEW_SNOW_ALB*pow(SNOW_ALB_ACCUM_A, pow((double)last_snow*dt/24.,
                   SNOW_ALB_ACCUM_B));

    /* Melt Season */
    else
      albedo = NEW_SNOW_ALB*pow(SNOW_ALB_THAW_A, pow((double)last_snow*dt/24.,
                   SNOW_ALB_THAW_B));

  }

  return(albedo);
}

#undef ETA0
#undef G
#undef C5     
#undef C6
