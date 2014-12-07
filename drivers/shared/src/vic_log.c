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
 * @brief    Get string of current date and time.  Format YYMMDD-SSSSS.
 *****************************************************************************/
char*
get_current_datetime()
{
    const int size = 20;
    char     *cdt = (char*)malloc(sizeof(char) * size);
    char      ymd[8];
    unsigned  seconds_since_midnight;

    if (cdt == NULL) {
        return NULL;
    }

    memset(cdt, 0, size);

    time_t currDateTime;
    currDateTime = time(NULL);

    if (currDateTime == -1) {
        return NULL;
    }

    struct tm *timeinfo = localtime(&currDateTime);

    seconds_since_midnight = (unsigned)currDateTime % 86400;

    if (strftime(ymd, 7, "%y%m%d", timeinfo) == 0) {
        return NULL;
    }

    sprintf(cdt, "%s-%05d", ymd, seconds_since_midnight);

    return cdt;
}

/******************************************************************************
 * @brief    Make logfile name string.
 *****************************************************************************/
char*
get_logname(const char *path)
{
    char *timestamp = get_current_datetime();
    char *prefix = "vic.log.";
    char *ext = ".txt";
    int   size = (strlen(path) + strlen(prefix) + strlen(ext) +
                  strlen(timestamp) + 1);
    char *filename = (char*)malloc(sizeof(char) * size);

    if (filename == NULL) {
        return NULL;
    }

    memset(filename, 0, size);
    strcpy(filename, path);
    strcpy(filename + strlen(path), prefix);
    strcpy(filename + strlen(path) + strlen(prefix), timestamp);
    strcpy(filename + strlen(path) + strlen(prefix) + strlen(timestamp),
           ext);

    return filename;
}

/******************************************************************************
 * @brief    Set global log destination
 *****************************************************************************/
void
initialize_log()
{
    extern FILE *LOG_DEST;

    LOG_DEST = stderr;
}

/******************************************************************************
 * @brief    Set global log destination
 *****************************************************************************/
void
setup_logging()
{
    extern filenames_struct filenames;
    extern filep_struct     filep;
    extern FILE            *LOG_DEST;

    char                   *logfilename;

    if (strcmp(filenames.log_path, "MISSING") != 0) {
        // Create logfile name
        logfilename = get_logname(filenames.log_path);

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
