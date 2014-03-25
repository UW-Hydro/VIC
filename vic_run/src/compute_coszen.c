#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
#include <math.h>
#include <mtclim_constants_vic.h>
#include <mtclim_parameters_vic.h>

static char vcid[] = "$Id: compute_coszen.c,v 5.7 2004/07/07 01:46:14 tbohn Exp $";

double compute_coszen(double lat,
                      double lng,
                      double time_zone_lng,
                      dmy_struct dmy)
/**********************************************************************
	compute_coszen.c	Ted Bohn		Feb 23, 2007

  This subroutine computes the cosine of the solar zenith angle, given
  the current location and date.

**********************************************************************/
{  
  extern option_struct options;
  double coslat;
  double sinlat;
  double decl;
  double cosdecl;
  double sindecl;
  double cosegeom;
  double sinegeom;
  double coshss;
  double hss;
  double hour_offset;
  double cosh;
  double coszen;

  /* calculate cos and sin of latitude */
  coslat = cos(lat*PI/180);
  sinlat = sin(lat*PI/180);

  /* calculate cos and sin of declination */
  decl = MINDECL * cos(((double)dmy.day_in_year + DAYSOFF) * RADPERDAY);
  cosdecl = cos(decl);
  sindecl = sin(decl);

  /* calculate daylength as a function of lat and decl */
  cosegeom = coslat * cosdecl;
  sinegeom = sinlat * sindecl;
  coshss = -(sinegeom) / cosegeom;
  if (coshss < -1.0)
    coshss = -1.0;  /* 24-hr daylight */
  if (coshss > 1.0)
    coshss = 1.0;    /* 0-hr daylight */
  hss = acos(coshss);                /* hour angle at sunset (radians) */
//  /* daylength (seconds) */
//  daylength[i] = 2.0 * hss * SECPERRAD;

  /* calculate cos of hour angle */
  hour_offset = (time_zone_lng - lng) * 24/360;
  cosh = cos(((double)dmy.hour + hour_offset - 12)*PI/12);

  /* calculate cosine of solar zenith angle */
  coszen = cosegeom * cosh + sinegeom;

  return coszen;

}
