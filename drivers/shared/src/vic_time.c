/******************************************************************************
 * @section DESCRIPTION
 *
 * VIC time and calendar module
 *
 * Ported from Unidata and Jeff Whitaker's netcdftime.py module:
 * https://github.com/Unidata/netcdf4-python/blob/master/netcdftime/netcdftime.py
 *
 * @section LICENSE
 *
 * The Variable Infiltration Capacity (VIC) macroscale hydrological model
 * Copyright (C) 2014 The Land Surface Hydrology Group, Department of Civil
 * and Environmental Engineering, University of Washington.
 *
 * The VIC model is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 *****************************************************************************/

#include <math.h>

#include <vic_def.h>
#include <vic_run.h>
#include <vic_driver_shared.h>

/******************************************************************************
 * @brief   Get fractional day of month from dmy structure.
 * @return  Fractional Day.
 *****************************************************************************/
double
fractional_day_from_dmy(dmy_struct *dmy)
{
    double day;

    day = (double) dmy->day + (double) dmy->dayseconds / (double) (SEC_PER_DAY);

    return day;
}

/******************************************************************************
 * @brief   Creates a Julian Day from a dmy_struct.
 * @return  Fractional Julian Day.
 * @note    if calendar="standard" or "gregorian" (default), Julian day follows
 *          Julian Calendar on and before 1582-10-5, Gregorian calendar after
 *          1582-10-15.
 *
 *          if calendar="proleptic_gregorian", Julian Day follows gregorian
 *          calendar.
 *
 *          if calendar="julian", Julian Day follows julian calendar.
 *
 *          Algorithm:
 *          Meeus, Jean (1998) Astronomical Algorithms (2nd Edition).
 *          Willmann-Bell, Virginia. p. 63
 *
 *          based on redate.py by David Finlayson.
 *****************************************************************************/
double
julian_day_from_dmy(dmy_struct    *dmy,
                    unsigned short calendar)
{
    double jd, day;
    int    year, month;
    int    A, B;

    year = dmy->year;
    month = dmy->month;

    // Convert time to fractions of a day
    day = fractional_day_from_dmy(dmy);

    // Start Meeus algorithm (variables are in his notation)
    if (month < 3) {
        month += MONTHS_PER_YEAR;
        year -= 1;
    }

    A = year / 100;

    jd = DAYS_PER_YEAR * year +
         (double)(int)(0.25 * (double)year + 2000.) +
         (double)(int)(30.6001 * (double)(month + 1)) +
         day + 1718994.5;

    // optionally adjust the jd for the switch from
    // the Julian to Gregorian Calendar
    // here assumed to have occurred the day after 1582 October 4
    if (calendar == CALENDAR_STANDARD || calendar == CALENDAR_GREGORIAN) {
        if (jd >= 2299170.5) {
            // 1582 October 15 (Gregorian Calendar)
            B = 2 - A + (int)(A / 4);
        }
        else if (jd < 2299160.5) {
            // 1582 October 5 (Julian Calendar)
            B = 0;
        }
        else {
            log_err("impossible date (falls in gap between end of Julian "
                    "calendar and beginning of Gregorian calendar");
        }
    }
    else if (calendar == CALENDAR_PROLEPTIC_GREGORIAN) {
        B = 2 - A + (int)(A / 4);
    }
    else if (calendar == CALENDAR_JULIAN) {
        B = 0;
    }
    else {
        log_err("unknown calendar, must be one of julian,standard,gregorian,"
                "proleptic_gregoria");
    }

    // adjust for Julian calendar if necessary
    jd += B;

    return jd;
}

/******************************************************************************
 * @brief   creates a Julian Day for a calendar with no leap years from a
 *          dmy structure.
 * @return  Fractional Julian Day.
 *****************************************************************************/
double
no_leap_day_from_dmy(dmy_struct *dmy)
{
    double         jd, day;
    unsigned short year, month;

    year = dmy->year;
    month = dmy->month;

    // Convert time to fractions of a day
    day = fractional_day_from_dmy(dmy);

    // Start Meeus algorithm (variables are in his notation)
    if (month < 3) {
        month += MONTHS_PER_YEAR;
        year -= 1;
    }

    jd = (double)((unsigned short)(DAYS_PER_YEAR * (year + 4716))) +
         (double)((unsigned short)(30.6001 * (month + 1))) +
         day - 1524.5;

    return jd;
}

/******************************************************************************
 * @brief   creates a Julian Day for a calendar where all years have 366 days
 *          from a dmy structure.
 * @return  Fractional Julian Day.
 *****************************************************************************/
double
all_leap_from_dmy(dmy_struct *dmy)
{
    double         jd, day;
    unsigned short year, month;

    year = dmy->year;
    month = dmy->month;

    // Convert time to fractions of a day
    day = fractional_day_from_dmy(dmy);

    // Start Meeus algorithm (variables are in his notation)
    if (month < 3) {
        month += MONTHS_PER_YEAR;
        year -= 1;
    }

    jd = (double)(unsigned short)(DAYS_PER_LYEAR * (year + 4716)) +
         (double)(unsigned short)(30.6001 * (month + 1)) +
         day - 1524.5;

    return jd;
}

/******************************************************************************
 * @brief   creates a Julian Day for a calendar where all months have 30 days
 *          from a dmy structure.
 * @return  Fractional Julian Day.
 *****************************************************************************/
double
all_30_day_from_dmy(dmy_struct *dmy)
{
    double         jd, day;
    unsigned short year, month;

    year = dmy->year;
    month = dmy->month;

    // Convert time to fractions of a day
    day = fractional_day_from_dmy(dmy);

    jd = (double)((unsigned short)(360. * (year + 4716))) +
         (double)((unsigned short)(30. * (month - 1))) +
         day;

    return jd;
}

/******************************************************************************
 * @brief   Calculate the day, month, year and seconds given a Fractional
 *          Julian Day.
 * @return  dmy stucture.
 * @note    if calendar="standard" or "gregorian" (default), Julian day follows
 *          Julian Calendar on and before 1582-10-5, Gregorian calendar after
 *          1582-10-15.
 *
 *          if calendar="proleptic_gregorian", Julian Day follows gregorian
 *          calendar.
 *
 *          if calendar="julian", Julian Day follows julian calendar.
 *
 *          Algorithm:
 *          Meeus, Jean (1998) Astronomical Algorithms (2nd Edition).
 *          Willmann-Bell, Virginia. p. 63
 *
 *          based on redate.py by David Finlayson.
 *****************************************************************************/
void
dmy_julian_day(double         julian,
               unsigned short calendar,
               dmy_struct    *dmy)
{
    double day, F, eps;
    double hour, minute, second;
    int    Z, alpha;
    int    A, B, C, D, E;
    int    nday, dayofyr;
    int    year, month;
    bool   is_leap;

    if (julian < 0) {
        log_err("Julian Day must be positive");
    }

    // get the day (Z) and the fraction of the day (F)
    // add 0.000005 which is 452 ms in case of jd being after
    // second 23:59:59 of a day we want to round to the next day
    Z = (int)round(julian + 0.00005);
    F = (double)(julian + 0.5 - Z);
    if (calendar == CALENDAR_STANDARD || calendar == CALENDAR_GREGORIAN) {
        alpha = (int)(((Z - 1867216.0) - 0.25) / 36524.25);
        A = Z + 1 + alpha - (int)(0.25 * alpha);
        // check if dates before oct 5th 1582 are in the array
        if (julian < 2299160.5) {
            A = Z;
        }
    }
    else if (calendar == CALENDAR_PROLEPTIC_GREGORIAN) {
        alpha = (int)(((Z - 1867216.0) - 0.25) / 36524.25);
        A = Z + 1 + alpha - (int)(0.25 * alpha);
    }
    else if (calendar == CALENDAR_JULIAN) {
        A = Z;
    }
    else {
        log_err("unknown calendar, must be one of julian,standard,gregorian,"
                "proleptic_gregorian");
    }

    B = A + 1524;
    C = (int)(6680. + ((B - 2439870.) - 122.1) / DAYS_PER_JYEAR);
    D = (int)(DAYS_PER_YEAR * C + (int)(0.25 * C));
    E = (int)((B - D) / 30.6001);

    // Convert to date
    day = B - D - (int)(30.6001 * E) + F;
    if (day < 1) {
        day = 1;
    }
    nday = B - D - 123;
    dayofyr = nday - 305;
    if (nday <= 305) {
        dayofyr = nday + 60;
    }

    month = E - 1;
    if (month > MONTHS_PER_YEAR) {
        month -= MONTHS_PER_YEAR;
    }
    year = C - 4715;
    if (month > 2) {
        year -= 1;
    }
    if (year <= 0) {
        year -= 1;
    }

    // a leap year?
    is_leap = leap_year(year, calendar);

    if (is_leap && (month > 2)) {
        dayofyr += 1;
    }

    eps = 1e-12 * abs(Z);
    if (eps < 1e-12) {
        eps = 1e-12;
    }

    hour = F * HOURS_PER_DAY + eps;
    if (hour < 0) {
        hour = 0;
    }
    if (hour > 23) {
        hour = 23;
    }
    F -= hour / HOURS_PER_DAY;

    minute = F * MIN_PER_DAY + eps;
    if (minute < 0) {
        minute = 0;
    }
    if (minute > 59) {
        minute = 59;
    }

    second = (F - minute / MIN_PER_DAY) * CONST_CDAY;
    if (second < 0) {
        second = 0;
    }

    dmy->day = day;
    dmy->day_in_year = dayofyr;
    dmy->month = month;
    dmy->year = year;
    dmy->dayseconds = hour * SEC_PER_HOUR + minute * SEC_PER_MIN + second;

    return;
}

/******************************************************************************
 * @brief   Calculate the day, month, year and seconds given a Fractional
 *          Julian Day for the NoLeap calendar.
 * @return  dmy structure.
 *****************************************************************************/
void
dmy_no_leap_day(double      julian,
                dmy_struct *dmy)
{
    double         F;
    double         day, dfrac;
    unsigned       A, B, C, D, E, year;
    double         I, days, seconds;
    unsigned short month, nday, dayofyr;

    if (julian < 0) {
        log_err("Julian Day must be positive");
    }

    F = modf(julian + 0.5, &I);
    A = (unsigned)I;
    B = A + 1524;
    C = (unsigned)((B - 122.1) / DAYS_PER_YEAR);
    D = (unsigned)(DAYS_PER_YEAR * C);
    E = (unsigned)((B - D) / 30.6001);

    // Convert to date
    day = B - D - (unsigned)(30.6001 * E) + (unsigned)F;
    nday = B - D - 123;
    if (nday <= 305) {
        dayofyr = nday + 60;
    }
    else {
        dayofyr = nday - 305;
    }

    if (E < 14) {
        month = E - 1;
    }
    else {
        month = E - 13;
    }

    if (month > 2) {
        year = C - 4716;
    }
    else {
        year = C - 4715;
    }

    // Convert fractions of a day to time
    dfrac = modf(day, &days);
    seconds = dfrac * CONST_CDAY;

    dmy->year = year;
    dmy->month = month;
    dmy->day = (unsigned short)days;
    dmy->day_in_year = dayofyr;
    dmy->dayseconds = seconds;

    return;
}

/******************************************************************************
 * @brief   Calculate the day, month, year and seconds given a Fractional
 *          Julian Day for the AllLeap calendar.
 * @return  dmy structure.
 * @note    based on redate.py by David Finlayson.
 *****************************************************************************/
void
dmy_all_leap(double      julian,
             dmy_struct *dmy)
{
    double         F, day, days, seconds;
    double         dfrac, I;
    unsigned       A, B, C, D, E, year, nday;
    unsigned short dayofyr, month;

    if (julian < 0) {
        log_err("Julian Day must be positive");
    }

    F = modf(julian + 0.5, &I);
    A = (unsigned)I;
    B = A + 1524;
    C = (unsigned)((B - 122.1) / DAYS_PER_LYEAR);
    D = (unsigned)(DAYS_PER_LYEAR * C);
    E = (unsigned)((B - D) / 30.6001);

    // Convert to date
    day = B - D - (unsigned)(30.6001 * E) + (unsigned)F;
    nday = B - D - 123;
    if (nday <= 305) {
        dayofyr = nday + 60;
    }
    else {
        dayofyr = nday - 305;
    }
    if (E < 14) {
        month = E - 1;
    }
    else {
        month = E - 13;
    }
    if (month > 2) {
        dayofyr += 1;
    }

    if (month > 2) {
        year = C - 4716;
    }
    else {
        year = C - 4715;
    }

    // Convert fractions of a day to time
    dfrac = modf(day, &days);
    seconds = dfrac * CONST_CDAY;

    dmy->year = year;
    dmy->month = month;
    dmy->day = (unsigned short)days;
    dmy->day_in_year = dayofyr;
    dmy->dayseconds = seconds;

    return;
}

/******************************************************************************
 * @brief   Calculate the day, month, year and seconds given a Fractional
 *          Julian Day for the 360day calendar.
 * @return  dmy structure.
 *****************************************************************************/
void
dmy_all_30_day(double      julian,
               dmy_struct *dmy)
{
    double         F, Z, dfrac, day, days, seconds;
    unsigned short dayofyr, month;
    unsigned       year;

    if (julian < 0) {
        log_err("Julian Day must be positive");
    }

    F = modf(julian, &Z);
    year = (int)((Z - 0.5) / 360.) - 4716;
    dayofyr = Z - (year + 4716) * 360;
    month = (int)((dayofyr - 0.5) / 30) + 1;
    day = dayofyr - (month - 1) * 30 + F;

    // Convert fractions of a day to time
    dfrac = modf(day, &days);
    seconds = dfrac * CONST_CDAY;

    dmy->year = year;
    dmy->month = month;
    dmy->day = (unsigned short)days;
    dmy->day_in_year = dayofyr;
    dmy->dayseconds = seconds;

    return;
}

/******************************************************************************
 * @brief   Calculate the numeric date given a dmy_structure.
 * @return  time value in units described by "time_units", using the
 *          specified "calendar"
 * @note    If calendar = "standard" or "gregorian" (indicating that the
 *          mixed Julian/Gregorian calendar is to be used), an error will be
 *          raised if the dmy structure describes a date between
 *          1582-10-5 and 1582-10-15.
 *****************************************************************************/
double
date2num(double         origin,
         dmy_struct    *dmy,
         double         tzoffset,
         unsigned short calendar,
         unsigned short time_units)
{
    double jdelta;

    if (calendar == CALENDAR_JULIAN ||
        calendar == CALENDAR_STANDARD ||
        calendar == CALENDAR_GREGORIAN ||
        calendar == CALENDAR_PROLEPTIC_GREGORIAN) {
        jdelta = julian_day_from_dmy(dmy, calendar) - origin;
    }
    else if (calendar == CALENDAR_NOLEAP || calendar == CALENDAR_365_DAY) {
        if ((dmy->month == 2) && (dmy->day == 29)) {
            log_err("there is no leap day in the noleap calendar");
        }
        jdelta = no_leap_day_from_dmy(dmy) - origin;
    }
    else if ((calendar == CALENDAR_ALL_LEAP) ||
             (calendar == CALENDAR_366_DAY)) {
        jdelta = all_leap_from_dmy(dmy) - origin;
    }
    else if (calendar == CALENDAR_360_DAY) {
        if (dmy->day > 30) {
            log_err("there are only 30 days in every month with the 360_day "
                    "calendar");
        }
        jdelta = all_30_day_from_dmy(dmy) - origin;
    }
    else {
        log_err("Unknown Calendar Flag: %hu", calendar);
    }


    // convert to desired units, add time zone offset.
    if (time_units == TIME_UNITS_SECONDS) {
        jdelta = jdelta * CONST_CDAY + tzoffset * SEC_PER_HOUR;
    }
    else if (time_units == TIME_UNITS_MINUTES) {
        jdelta = jdelta * MIN_PER_DAY + tzoffset * MIN_PER_HOUR;
    }
    else if (time_units == TIME_UNITS_HOURS) {
        jdelta = jdelta * HOURS_PER_DAY + tzoffset;
    }
    else if (time_units == TIME_UNITS_DAYS) {
        jdelta = jdelta + tzoffset / HOURS_PER_DAY;
    }
    else {
        log_err("Unknown Time Units Flag: %hu", time_units);
    }
    return jdelta;
}

/******************************************************************************
 * @brief   Create a dmy structure given numeric date
 * @return  Return a dmy structure given a numeric "time_value" in "time_units"
            using "calendar".
 *****************************************************************************/
void
num2date(double         origin,
         double         time_value,
         double         tzoffset,
         unsigned short calendar,
         unsigned short time_units,
         dmy_struct    *dmy)
{
    double jd, jdelta;

    if (time_units == TIME_UNITS_SECONDS) {
        jdelta = time_value / CONST_CDAY - tzoffset / HOURS_PER_DAY;
    }
    else if (time_units == TIME_UNITS_MINUTES) {
        jdelta = time_value / MIN_PER_DAY - tzoffset / HOURS_PER_DAY;
    }
    else if (time_units == TIME_UNITS_HOURS) {
        jdelta = time_value / HOURS_PER_DAY - tzoffset / HOURS_PER_DAY;
    }
    else if (time_units == TIME_UNITS_DAYS) {
        jdelta = time_value - tzoffset / HOURS_PER_DAY;
    }
    else {
        log_err("Unknown Time Units Flag: %hu", time_units);
    }

    jd = jdelta + origin;

    if (calendar == CALENDAR_JULIAN ||
        calendar == CALENDAR_STANDARD ||
        calendar == CALENDAR_GREGORIAN ||
        calendar == CALENDAR_PROLEPTIC_GREGORIAN) {
        dmy_julian_day(jd, calendar, dmy);
    }
    else if (calendar == CALENDAR_NOLEAP || calendar == CALENDAR_365_DAY) {
        dmy_no_leap_day(jd, dmy);
    }
    else if (calendar == CALENDAR_ALL_LEAP || calendar == CALENDAR_366_DAY) {
        dmy_all_leap(jd, dmy);
    }
    else if (calendar == CALENDAR_360_DAY) {
        dmy_all_30_day(jd, dmy);
    }
    else {
        log_err("Unknown Calendar Flag: %hu", calendar);
    }

    return;
}

/******************************************************************************
 * @brief   Make array of last day of months
 * @return  Return pointer to last day of month array
 *****************************************************************************/
void
make_lastday(unsigned short calendar,
             unsigned short year,
             unsigned short lastday[])
{
    size_t         i;
    unsigned short temp[MONTHS_PER_YEAR] = {
        31, 28, 31, 30, 31, 30,
        31, 31, 30, 31, 30, 31
    };

    for (i = 0; i < MONTHS_PER_YEAR; i++) {
        if (calendar == CALENDAR_360_DAY) {
            lastday[i] = 30;
        }
        else {
            lastday[i] = temp[i];
        }
    }
    // leap year?
    if (calendar == CALENDAR_JULIAN ||
        calendar == CALENDAR_STANDARD ||
        calendar == CALENDAR_GREGORIAN ||
        calendar == CALENDAR_PROLEPTIC_GREGORIAN) {
        if (leap_year(year, calendar)) {
            lastday[1] = 29;
        }
    }
    else if (calendar == CALENDAR_366_DAY ||
             calendar == CALENDAR_ALL_LEAP) {
        lastday[1] = 29;
    }
    return;
}

/******************************************************************************
 * @brief   Make array of last day of months
 * @return  Return pointer to last day of month array
 *****************************************************************************/
bool
leap_year(unsigned short year,
          unsigned short calendar)
{
    bool leap = false;

    if ((calendar == CALENDAR_JULIAN ||
         calendar == CALENDAR_STANDARD ||
         calendar == CALENDAR_GREGORIAN ||
         calendar == CALENDAR_PROLEPTIC_GREGORIAN) &&
        (year % 4 == 0)) {
        leap = true;
        if ((calendar == CALENDAR_PROLEPTIC_GREGORIAN) &&
            (year % 100 == 0) && (year % 400 != 0)) {
            leap = false;
        }
        else if ((calendar == CALENDAR_STANDARD ||
                  calendar == CALENDAR_GREGORIAN) &&
                 (year % 100 == 0) && (year % 400 != 0) && (year < 1583)) {
            leap = false;
        }
    }
    return leap;
}

/******************************************************************************
 * @brief   Initialize global time origin
 *****************************************************************************/
void
initialize_time()
{
    extern global_param_struct global_param;

    dmy_struct                 dmy;

    // Origin is 0001-01-01
    dmy.year = 1;
    dmy.month = 1;
    dmy.day = 1;
    dmy.day_in_year = 1;
    dmy.dayseconds = 0;

    // Set origin using date2num with numeric origin of 0.
    // This calculates the julian day (in units = global_param.time_units)
    // for 0001-01-01 for the given calendar
    // This is later used as the reference (origin for advancing numeric dates).
    // See make_dmy.c for more details on how global_param.time_origin_num is used.
    global_param.time_origin_num = date2num(0., &dmy, 0., global_param.calendar,
                                            global_param.time_units);
    return;
}

/******************************************************************************
 * @brief Validate dmy structure
 * @return 0 if ok, else return > 0
 *****************************************************************************/

int
valid_date(unsigned short calendar,
           dmy_struct    *dmy)
{
    unsigned short lastday[MONTHS_PER_YEAR];
    unsigned short days_in_year;
    size_t         i;

    // Get array of last days of month
    make_lastday(calendar, dmy->year, lastday);

    // Calculate the sum of lastday (days in year)
    days_in_year = 0;
    for (i = 0; i < MONTHS_PER_YEAR; i++) {
        days_in_year += lastday[i];
    }

    if (dmy->dayseconds > SEC_PER_DAY) {
        return 1;
    }
    else if (dmy->month > MONTHS_PER_YEAR) {
        return 2;
    }
    else if (dmy->month < 1) {
        return 3;
    }
    else if (dmy->day > lastday[dmy->month - 1]) {
        return 4;
    }
    else if (dmy->day < 1) {
        return 5;
    }
    else if (dmy->day_in_year > days_in_year) {
        return 6;
    }
    else if (dmy->day_in_year < 1) {
        return 7;
    }
    else {
        return 0;
    }
}

void
dt_seconds_to_time_units(unsigned short int time_units,
                         double             dt_seconds,
                         double            *dt_time_units)
{
    if (time_units == TIME_UNITS_SECONDS) {
        *dt_time_units = dt_seconds;
    }
    else if (time_units == TIME_UNITS_MINUTES) {
        *dt_time_units = dt_seconds / SEC_PER_MIN;
    }
    else if (time_units == TIME_UNITS_HOURS) {
        *dt_time_units = dt_seconds / SEC_PER_HOUR;
    }
    else if (time_units == TIME_UNITS_DAYS) {
        *dt_time_units = dt_seconds / SEC_PER_DAY;
    }
    else {
        log_err("Unknown Time Units Flag: %hu", time_units);
    }
}
