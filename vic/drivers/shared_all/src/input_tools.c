/******************************************************************************
 * @section DESCRIPTION
 *
 * This file includes routines to help parse and process VIC input files.
 *
 * @section LICENSE
 *
 * The Variable Infiltration Capacity (VIC) macroscale hydrological model
 * Copyright (C) 2016 The Computational Hydrology Group, Department of Civil
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

#include <vic_driver_shared_all.h>

/******************************************************************************
 * @brief   String to bool conversion
 *****************************************************************************/
bool
str_to_bool(char str[])
{
    if (strcasecmp("TRUE", str) == 0) {
        return true;
    }
    else if (strcasecmp("FALSE", str) == 0) {
        return false;
    }
    else {
        log_err("%s is neither TRUE nor FALSE", str);
    }
    return false; // To avoid warnings.
}

/******************************************************************************
 * @brief    This routine determines the counts the number of output variables
             in each output file specified in the global parameter file.
 *****************************************************************************/
void
count_nstreams_nvars(FILE   *gp,
                     size_t *nstreams,
                     size_t  nvars[])
{
    unsigned long start_position;
    char          cmdstr[MAXSTRING];
    char          optstr[MAXSTRING];
    size_t        i;

    // Figure out where we are in the input file
    fflush(gp);
    start_position = ftell(gp);

    // Move the position to the begining of the file
    rewind(gp);

    // read the first line
    fgets(cmdstr, MAXSTRING, gp);

    // initialize nstreams and nvars
    *nstreams = 0;
    for (i = 0; i < MAX_OUTPUT_STREAMS; i++) {
        nvars[i] = 0;
    }

    // Loop through the lines
    while (!feof(gp)) {
        if (cmdstr[0] != '#' && cmdstr[0] != '\n' && cmdstr[0] != '\0') {
            // line is not blank or a comment
            sscanf(cmdstr, "%s", optstr);

            // if the line starts with OUTFILE, increment nstreams
            if (strcasecmp("OUTFILE", optstr) == 0) {
                (*nstreams)++;
            }

            // if the line starts with OUTVAR, add another variable to nvars
            if (strcasecmp("OUTVAR", optstr) == 0) {
                nvars[*nstreams - 1]++;
            }
        }
        fgets(cmdstr, MAXSTRING, gp);
    }

    if (*nstreams > MAX_OUTPUT_STREAMS) {
        log_err("Too many output streams specified.");
    }

    // put the position in the file back to where we started
    fseek(gp, start_position, SEEK_SET);
}

/******************************************************************************
 * @brief    Convert string version of AGG_TYPE_* to enum value
 *****************************************************************************/
unsigned short int
str_to_agg_type(char aggstr[])
{
    if ((strcasecmp("", aggstr) == 0) || (strcasecmp("*", aggstr) == 0)) {
        return AGG_TYPE_DEFAULT;
    }
    else {
        if (strcasecmp("AGG_TYPE_AVG", aggstr) == 0) {
            return AGG_TYPE_AVG;
        }
        else if (strcasecmp("AGG_TYPE_BEG", aggstr) == 0) {
            return AGG_TYPE_BEG;
        }
        else if (strcasecmp("AGG_TYPE_END", aggstr) == 0) {
            return AGG_TYPE_END;
        }
        else if (strcasecmp("AGG_TYPE_MAX", aggstr) == 0) {
            return AGG_TYPE_MAX;
        }
        else if (strcasecmp("AGG_TYPE_MIN", aggstr) == 0) {
            return AGG_TYPE_MIN;
        }
        else if (strcasecmp("AGG_TYPE_SUM", aggstr) == 0) {
            return AGG_TYPE_SUM;
        }
        else {
            log_err("Unknown aggregation type found: %s", aggstr);
        }
    }
    return 0; // To avoid warnings.
}

/******************************************************************************
 * @brief    Convert string version of OUT_TYPE* to enum value
 *****************************************************************************/
unsigned short int
str_to_out_type(char typestr[])
{
    if ((strcasecmp("", typestr) == 0) || (strcasecmp("*", typestr) == 0)) {
        return OUT_TYPE_DEFAULT;
    }
    else {
        if (strcasecmp("OUT_TYPE_USINT", typestr) == 0) {
            return OUT_TYPE_USINT;
        }
        else if (strcasecmp("OUT_TYPE_SINT", typestr) == 0) {
            return OUT_TYPE_SINT;
        }
        else if (strcasecmp("OUT_TYPE_INT", typestr) == 0) {
            return OUT_TYPE_INT;
        }
        else if (strcasecmp("OUT_TYPE_CHAR", typestr) == 0) {
            return OUT_TYPE_CHAR;
        }
        else if (strcasecmp("OUT_TYPE_FLOAT", typestr) == 0) {
            return OUT_TYPE_FLOAT;
        }
        else if (strcasecmp("OUT_TYPE_DOUBLE", typestr) == 0) {
            return OUT_TYPE_DOUBLE;
        }
        else {
            log_err("Unknown out type found: %s", typestr);
        }
    }
    return 0; // To avoid warnings.
}

/******************************************************************************
 * @brief    Convert string version of mult to double
 *****************************************************************************/
double
str_to_out_mult(char multstr[])
{
    if ((strcasecmp("", multstr) == 0) || (strcasecmp("*", multstr) == 0)) {
        return OUT_MULT_DEFAULT;
    }
    else {
        return (double) atof(multstr);
    }
    return 0.; // To avoid warnings.
}

/******************************************************************************
 * @brief    Convert string version of frequency flags to enum value
 *****************************************************************************/
unsigned short int
str_to_freq_flag(char freq[])
{
    if (strcasecmp("NEVER", freq) == 0) {
        return FREQ_NEVER;
    }
    else if (strcasecmp("NSTEPS", freq) == 0) {
        return FREQ_NSTEPS;
    }
    else if (strcasecmp("NSECONDS", freq) == 0) {
        return FREQ_NSECONDS;
    }
    else if (strcasecmp("NMINUTES", freq) == 0) {
        return FREQ_NMINUTES;
    }
    else if (strcasecmp("NHOURS", freq) == 0) {
        return FREQ_NHOURS;
    }
    else if (strcasecmp("NDAYS", freq) == 0) {
        return FREQ_NDAYS;
    }
    else if (strcasecmp("NMONTHS", freq) == 0) {
        return FREQ_NMONTHS;
    }
    else if (strcasecmp("NYEARS", freq) == 0) {
        return FREQ_NYEARS;
    }
    else if (strcasecmp("DATE", freq) == 0) {
        return FREQ_DATE;
    }
    else if (strcasecmp("END", freq) == 0) {
        return FREQ_END;
    }
    else {
        log_err("Unknown frequency flag found: %s", freq);
    }
    return 0; // To avoid warnings.
}

/******************************************************************************
 * @brief    Convert string version of frequency flags to enum value
 *****************************************************************************/
void
str_to_ascii_format(char *format)
{
    if ((strcasecmp("", format) == 0) || (strcasecmp("*", format) == 0)) {
        strcpy(format, OUT_ASCII_FORMAT_DEFAULT);
    }
    // else do nothing
}

/******************************************************************************
 * @brief  Parse chars of calendar and return calendar integer
 * @return enum integer representing calendar
 *****************************************************************************/
unsigned short int
str_to_calendar(char *cal_chars)
{
    if (strcasecmp("STANDARD", cal_chars) == 0) {
        return CALENDAR_STANDARD;
    }
    else if (strcasecmp("GREGORIAN", cal_chars) == 0) {
        return CALENDAR_GREGORIAN;
    }
    else if (strcasecmp("PROLEPTIC_GREGORIAN", cal_chars) == 0) {
        return CALENDAR_PROLEPTIC_GREGORIAN;
    }
    else if ((strcasecmp("NOLEAP", cal_chars) == 0) ||
             (strcasecmp("NO_LEAP", cal_chars) == 0)) {
        return CALENDAR_NOLEAP;
    }
    else if (strcasecmp("365_DAY", cal_chars) == 0) {
        return CALENDAR_365_DAY;
    }
    else if (strcasecmp("360_DAY", cal_chars) == 0) {
        return CALENDAR_360_DAY;
    }
    else if (strcasecmp("JULIAN", cal_chars) == 0) {
        return CALENDAR_JULIAN;
    }
    else if (strcasecmp("ALL_LEAP", cal_chars) == 0) {
        return CALENDAR_ALL_LEAP;
    }
    else if (strcasecmp("366_DAY", cal_chars) == 0) {
        return CALENDAR_366_DAY;
    }
    else {
        log_err("Unknown calendar specified: %s", cal_chars);
    }
    return 0; // To avoid warnings.
}

/******************************************************************************
 * @brief  Parse chars of time units and return time units integer
 * @return enum integer representing time units
 *****************************************************************************/
unsigned short int
str_to_timeunits(char units_chars[])
{
    if (strcasecmp("SECONDS", units_chars) == 0) {
        return TIME_UNITS_SECONDS;
    }
    else if (strcasecmp("MINUTES", units_chars) == 0) {
        return TIME_UNITS_MINUTES;
    }
    else if (strcasecmp("HOURS", units_chars) == 0) {
        return TIME_UNITS_HOURS;
    }
    else if (strcasecmp("DAYS", units_chars) == 0) {
        return TIME_UNITS_DAYS;
    }
    else {
        log_err("Unknown time units specified: %s", units_chars);
    }
    return 0; // To avoid warnings.
}

/******************************************************************************
 * @brief  Convert enum time units to string
 *****************************************************************************/
void
str_from_time_units(unsigned short int time_units,
                    char              *unit_str)
{
    if (time_units == TIME_UNITS_SECONDS) {
        sprintf(unit_str, "seconds");
    }
    else if (time_units == TIME_UNITS_MINUTES) {
        sprintf(unit_str, "minutes");
    }
    else if (time_units == TIME_UNITS_HOURS) {
        sprintf(unit_str, "hours");
    }
    else if (time_units == TIME_UNITS_DAYS) {
        sprintf(unit_str, "days");
    }
    else {
        log_err("Invalid value, or no value for OUT_TIME_UNITS (%d).",
                time_units);
    }
}

/******************************************************************************
 * @brief  get CF calendar string based on enum calendar value
 *****************************************************************************/
void
str_from_calendar(unsigned short int calendar,
                  char              *calendar_str)
{
    if (calendar == CALENDAR_STANDARD) {
        sprintf(calendar_str, "standard");
    }
    else if (calendar == CALENDAR_GREGORIAN) {
        sprintf(calendar_str, "gregorian");
    }
    else if (calendar == CALENDAR_PROLEPTIC_GREGORIAN) {
        sprintf(calendar_str, "proleptic_gregorian");
    }
    else if (calendar == CALENDAR_NOLEAP) {
        sprintf(calendar_str, "noleap");
    }
    else if (calendar == CALENDAR_365_DAY) {
        sprintf(calendar_str, "365_day");
    }
    else if (calendar == CALENDAR_360_DAY) {
        sprintf(calendar_str, "360_day");
    }
    else if (calendar == CALENDAR_JULIAN) {
        sprintf(calendar_str, "julian");
    }
    else if (calendar == CALENDAR_ALL_LEAP) {
        sprintf(calendar_str, "all_leap");
    }
    else if (calendar == CALENDAR_366_DAY) {
        sprintf(calendar_str, "366_day");
    }
    else {
        log_err("Invalid, or no calendar specified");
    }
}

/******************************************************************************
 * @brief  determine if the aggtype is a cell_method, if so fill the cell_method
           argument with the "cell_method" CF string attriubute
 * @return bool true if aggtype is a windowed aggregation, otherwise false
 *****************************************************************************/
bool
cell_method_from_agg_type(unsigned short int aggtype,
                          char               cell_method[])
{
    if (aggtype == AGG_TYPE_AVG) {
        strcpy(cell_method, "time: mean");
        return true;
    }
    else if (aggtype == AGG_TYPE_MAX) {
        strcpy(cell_method, "time: maximum");
        return true;
    }
    else if (aggtype == AGG_TYPE_MIN) {
        strcpy(cell_method, "time: minimum");
        return true;
    }
    else if (aggtype == AGG_TYPE_SUM) {
        strcpy(cell_method, "time: sum");
        return true;
    }
    else if (aggtype == AGG_TYPE_END) {
        strcpy(cell_method, "time: end");
        return true;
    }
    else if (aggtype == AGG_TYPE_BEG) {
        strcpy(cell_method, "time: beg");
        return true;
    }
    else {
        return false;
    }
    return false; // To avoid warnings.
}
