/******************************************************************************
 * @section DESCRIPTION
 *
 * calculates nitrogen scaling factors for all canopy layers, following eqns
 * 106 and 107 in Knorr 1997.
 *
 * Note: this should only be applied to veg types that have a canopy, e.g.
 * trees and shrubs, but not grass or tundra vegetation.
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

#include <time.h>
#include <vic_def.h>
#include <vic_run.h>
#include <vic_driver_shared.h>

/******************************************************************************
 * @brief    Finalize logging - called after all logging is completed
 *****************************************************************************/
void
finalize_logging(void)
{
    extern FILE *LOG_DEST;

    if (!(LOG_DEST == stdout || LOG_DEST == stderr)) {
        fclose(LOG_DEST);
        LOG_DEST = stderr;
    }
}

/******************************************************************************
 * @brief    Get string of current date and time.  Format YYMMDD-SSSSS.
 *****************************************************************************/
void
get_current_datetime(char *cdt)
{
    char         ymd[8];
    struct tm    timeinfo;
    unsigned int seconds_since_midnight;
    time_t       curr_date_time;

    curr_date_time = time(NULL);
    if (curr_date_time == -1) {
        return;
    }

    localtime_r(&curr_date_time, &timeinfo);

    seconds_since_midnight = (unsigned int) curr_date_time % CONST_CDAY;

    if (strftime(ymd, 7, "%y%m%d", &timeinfo) == 0) {
        return;
    }

    sprintf(cdt, "%s-%05d", ymd, seconds_since_midnight);
}

/******************************************************************************
 * @brief    Make logfile name string.
 *****************************************************************************/
void
get_logname(const char *path,
            int         id,
            char       *filename)
{
    char  timestamp[MAXSTRING];
    char *ext = ".txt";
    char *prefix = "vic.log.";

    memset(timestamp, 0, MAXSTRING);
    get_current_datetime(timestamp);

    memset(filename, 0, MAXSTRING);
    if (id != MISSING) {
        snprintf(filename, MAXSTRING - 1, "%s%s%s.%06d%s", path, prefix,
                 timestamp, id, ext);
    }
    else {
        snprintf(filename, MAXSTRING - 1, "%s%s%s%s", path, prefix,
                 timestamp, ext);
    }
}

/******************************************************************************
 * @brief    Set global log destination
 *****************************************************************************/
void
initialize_log(void)
{
    extern FILE *LOG_DEST;

    LOG_DEST = stderr;
}

/******************************************************************************
 * @brief    Set global log destination
 *****************************************************************************/
void
setup_logging(int id)
{
    extern filenames_struct filenames;
    extern filep_struct     filep;
    extern FILE            *LOG_DEST;
    char                    logfilename[MAXSTRING];

    if (strcmp(filenames.log_path, "MISSING") != 0) {
        // Create logfile name
        get_logname(filenames.log_path, id, logfilename);

        // Open Logfile
        filep.logfile = open_file(logfilename, "w");

        // Print log file name to stderr
        log_info("Initialized Log File: %s", logfilename);

        // Set Log Destination
        LOG_DEST = filep.logfile;

        // Write first line of log file
        log_info("Initialized Log File: %s", logfilename);
    }
    else {
        log_info("Logging to stderr");
    }
}
