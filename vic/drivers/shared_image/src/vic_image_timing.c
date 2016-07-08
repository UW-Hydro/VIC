/******************************************************************************
 * @section DESCRIPTION
 *
 * Write vic timing table for image-like drivers
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

#include <vic_driver_shared_image.h>

/******************************************************************************
 * @brief    VIC timing file
 *****************************************************************************/
void
write_vic_timing_table(timer_struct *timers,
                       char         *driver)
{
    extern FILE               *LOG_DEST;
    extern filenames_struct    filenames;
    extern global_param_struct global_param;
    extern int                 mpi_size;

    char                       machine[MAXSTRING];
    char                       user[MAXSTRING];
    time_t                     curr_date_time;
    struct tm                 *timeinfo;
    uid_t                      uid;
    struct passwd             *pw;
    double                     ndays;
    double                     nyears;

    // datestr
    curr_date_time = time(NULL);
    if (curr_date_time == -1) {
        log_err("Something went wrong getting the current time!");
    }
    timeinfo = localtime(&curr_date_time);

    // hostname
    if (gethostname(machine, MAXSTRING) != 0) {
        strcpy(machine, "unknown");
    }

    // username
    uid = geteuid();
    pw = getpwuid(uid);

    if (pw) {
        strcpy(user, pw->pw_name);
    }
    else {
        strcpy(user, "unknown");
    }

    // calculate run length
    ndays = global_param.dt * global_param.nrecs / SEC_PER_DAY;
    nyears = ndays / DAYS_PER_YEAR;

    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST,
            "------------------------------"
            " VIC TIMING PROFILE "
            "------------------------------\n\n");
    fprintf(LOG_DEST, "  Date                      : %s", asctime(timeinfo));
    fprintf(LOG_DEST, "  Compiler                  : %s (%s)\n", COMPILER,
            COMPILER_VERSION);
    fprintf(LOG_DEST, "  Machine                   : %s\n", machine);
    fprintf(LOG_DEST, "  VIC User                  : %s\n", user);
    fprintf(LOG_DEST, "  VIC Version               : %s\n", GIT_VERSION);
    fprintf(LOG_DEST, "  VIC GIT Version           : %s\n", VERSION);
    fprintf(LOG_DEST, "  VIC_DRIVER                : %s\n", driver);
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "  Global Param File         : %s\n", filenames.global);
    fprintf(LOG_DEST, "  Domain File               : %s\n", filenames.domain);
    fprintf(LOG_DEST, "  Start Date                : %04hu-%02hu-%02hu-%05u\n",
            global_param.startyear, global_param.startmonth,
            global_param.startday, global_param.startsec);
    fprintf(LOG_DEST, "  Stop Date                 : %04hu-%02hu-%02hu\n",
            global_param.endyear, global_param.endmonth, global_param.endday);
    fprintf(LOG_DEST, "  Nrecs                     : %zu\n",
            global_param.nrecs);
    fprintf(LOG_DEST, "  Model Timestep (seconds)  : %lf\n", global_param.dt);
    fprintf(LOG_DEST, "  Snow Timestep (seconds)   : %lf\n",
            global_param.snow_dt);
    fprintf(LOG_DEST, "  Runoff Timestep (seconds) : %lf\n",
            global_param.runoff_dt);
    fprintf(LOG_DEST, "  Atmos Timestep (seconds)  : %lf\n",
            global_param.atmos_dt);
    fprintf(LOG_DEST, "\n");

    fprintf(LOG_DEST, "  Total pes active          : %d\n", mpi_size);
    fprintf(LOG_DEST, "  pes per node              : %ld\n",
            sysconf(_SC_NPROCESSORS_ONLN));

    fprintf(LOG_DEST, "\n");

    fprintf(LOG_DEST, "  Overall Metrics\n");
    fprintf(LOG_DEST, "  ---------------\n");
    fprintf(LOG_DEST, "    Model Cost       : %lf pe-hrs/simulated_year\n",
            mpi_size * timers[TIMER_VIC_ALL].delta_wall / SEC_PER_HOUR /
            nyears);
    fprintf(LOG_DEST, "    Model Throughput : %lf simulated_years/day\n",
            nyears / (timers[TIMER_VIC_ALL].delta_wall / SEC_PER_DAY));
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "    Total Time       : %lf seconds \n",
            timers[TIMER_VIC_ALL].delta_wall);
    fprintf(LOG_DEST, "    Init Time        : %lf seconds \n",
            timers[TIMER_VIC_INIT].delta_wall);
    fprintf(LOG_DEST, "    Run Time         : %lf seconds  %lf seconds/day \n",
            timers[TIMER_VIC_RUN].delta_wall,
            timers[TIMER_VIC_RUN].delta_wall / ndays);
    fprintf(LOG_DEST, "    Final time       : %lf seconds \n",
            timers[TIMER_VIC_FINAL].delta_wall);
    fprintf(LOG_DEST, "\n");

    // fprintf(LOG_DEST, "Timing table for individual timers\n");
    // fprintf(LOG_DEST, "----------------------------------\n");
    // TODO...

    fprintf(LOG_DEST,
            "\n------------------------------"
            " END VIC TIMING PROFILE "
            "------------------------------\n\n");

    free(timeinfo);
}
