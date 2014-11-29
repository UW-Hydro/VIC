/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine checks the command line for valid program options.
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

#include <unistd.h>
#include <vic_def.h>
#include <vic_run.h>
#include <vic_driver_shared.h>

char *optstring = "g:vo";

/**********************************************************************
   cmd_proc                  Keith Cherkauer                1997

   This routine checks the command line for valid program options.  If
   no options are found, or an invalid combination of them appear, the
   routine calls usage() to print the model usage to the screen, before
   exiting execution.
**********************************************************************/
void
cmd_proc(int    argc,
         char **argv,
         char  *globalfilename)
{
    char GLOBAL_SET;
    int  optchar;

    if (argc == 1) {
        usage(argv[0]);
        exit(1);
    }

    GLOBAL_SET = FALSE;

    while ((optchar = getopt(argc, argv, optstring)) != EOF) {
        switch ((char)optchar) {
        case 'v':
            /** Version information **/
            display_current_settings(DISP_VERSION, NULL, NULL);
            exit(0);
            break;
        case 'o':
            /** Compile-time options information **/
            display_current_settings(DISP_COMPILE_TIME, NULL, NULL);
            exit(0);
            break;
        case 'g':
            /** Global Parameters File **/
            strncpy(globalfilename, optarg, MAXSTRING);
            GLOBAL_SET = TRUE;
            break;
        default:
            /** Print Usage if Invalid Command Line Arguments **/
            usage(argv[0]);
            exit(1);
            break;
        }
    }

    if (!GLOBAL_SET) {
        fprintf(stderr,
                "ERROR: Must set global control file using the '-g' flag\n");
        usage(argv[0]);
        exit(1);
    }
}

/******************************************************************************
 * @brief    This routine prints out usage details.
 *****************************************************************************/
void
usage(char *executable)
{
    fprintf(stderr,
            "Usage: %s [-v | -o | -g<global_parameter_file>]\n", executable);
    fprintf(stderr, "  v: display version information\n");
    fprintf(stderr,
            "  o: display compile-time options settings"
            " (set in vicNl_def.h)\n");
    fprintf(stderr,
            "  g: read model parameters from <global_parameter_file>.\n");
    fprintf(stderr,
            "       <global_parameter_file> is a file that contains all"
            " needed model\n");
    fprintf(stderr,
            "       parameters as well as model option flags, and the names"
            " and\n");
    fprintf(stderr, "       locations of all other files.\n");
}
