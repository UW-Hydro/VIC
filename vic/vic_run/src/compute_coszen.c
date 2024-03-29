/******************************************************************************
 * @section DESCRIPTION
 *
 * This subroutine computes the cosine of the solar zenith angle, given the
 * current location and date.
 *****************************************************************************/

#include <vic_run.h>

/******************************************************************************
 * @brief    This subroutine computes the cosine of the solar zenith angle.
 *****************************************************************************/
double
compute_coszen(double         lat,
               double         lng,
               double         time_zone_lng,
               unsigned short day_in_year,
               unsigned       second)
{
    double coslat;
    double sinlat;
    double decl;
    double cosdecl;
    double sindecl;
    double cosegeom;
    double sinegeom;
    double coshss;
    double hour_offset;
    double cosh;
    double coszen;
    double hour;

    hour = second / SEC_PER_HOUR;

    /* calculate cos and sin of latitude */
    coslat = cos(lat * CONST_PI / 180);
    sinlat = sin(lat * CONST_PI / 180);

    /* calculate cos and sin of declination */
    decl = CONST_MINDECL * cos(((double) day_in_year + CONST_DAYSOFF) *
                               CONST_RADPERDAY);
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
    hour_offset = (time_zone_lng - lng) * HOURS_PER_DAY / 360;
    cosh = cos((hour + hour_offset - 12) * CONST_PI / 12);

    /* calculate cosine of solar zenith angle */
    coszen = cosegeom * cosh + sinegeom;

    return coszen;
}
