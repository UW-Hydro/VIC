#include <vic_def.h>
#include <vic_run.h>
#include <mtclim_constants_vic.h>

double
compute_coszen(double     lat,
               double     lng,
               double     time_zone_lng,
               dmy_struct dmy)

/**********************************************************************
        compute_coszen.c	Ted Bohn		Feb 23, 2007

   This subroutine computes the cosine of the solar zenith angle, given
   the current location and date.

**********************************************************************/
{
    double               coslat;
    double               sinlat;
    double               decl;
    double               cosdecl;
    double               sindecl;
    double               cosegeom;
    double               sinegeom;
    double               coshss;
    double               hour_offset;
    double               cosh;
    double               coszen;

    /* calculate cos and sin of latitude */
    coslat = cos(lat * PI / 180);
    sinlat = sin(lat * PI / 180);

    /* calculate cos and sin of declination */
    decl = MINDECL * cos(((double)dmy.day_in_year + DAYSOFF) * RADPERDAY);
    cosdecl = cos(decl);
    sindecl = sin(decl);

    /* calculate daylength as a function of lat and decl */
    cosegeom = coslat * cosdecl;
    sinegeom = sinlat * sindecl;
    coshss = -(sinegeom) / cosegeom;
    if (coshss < -1.0) {
        coshss = -1.0; /* 24-hr daylight */
    }
    if (coshss > 1.0) {
        coshss = 1.0; /* 0-hr daylight */
    }

    /* calculate cos of hour angle */
    hour_offset = (time_zone_lng - lng) * 24 / 360;
    cosh = cos(((double)dmy.hour + hour_offset - 12) * PI / 12);

    /* calculate cosine of solar zenith angle */
    coszen = cosegeom * cosh + sinegeom;

    return coszen;
}
