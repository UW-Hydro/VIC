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
#include <getopt.h>
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
        print_usage(argv[0]);
        exit(1);
    }

    GLOBAL_SET = false;

    while ((optchar = getopt(argc, argv, optstring)) != EOF) {
        switch ((char)optchar) {
        case 'v':
            /** Version information **/
            display_current_settings(DISP_VERSION);
            exit(0);
            break;
        case 'o':
            /** Compile-time options information **/
            display_current_settings(DISP_COMPILE_TIME);
            exit(0);
            break;
        case 'g':
            /** Global Parameters File **/
            strncpy(globalfilename, optarg, MAXSTRING);
            GLOBAL_SET = true;
            break;
        default:
            /** Print Usage if Invalid Command Line Arguments **/
            print_usage(argv[0]);
            exit(1);
            break;
        }
    }

    if (!GLOBAL_SET) {
        fprintf(stderr,
                "ERROR: Must set global control file using the '-g' flag\n");
        print_usage(argv[0]);
        exit(1);
    }
}

/******************************************************************************
 * @brief    This routine prints out usage details.
 *****************************************************************************/
void
print_usage(char *executable)
{
    fprintf(stdout,
            "Usage: %s [-v | -o | -g <global_parameter_file>]\n", executable);
    fprintf(stdout, "  v: display version information\n");
    fprintf(stdout,
            "  o: display compile-time options settings"
            " (set in .h files)\n");
    fprintf(stdout,
            "  g: read model parameters from <global_parameter_file>.\n");
    fprintf(stdout,
            "       <global_parameter_file> is a file that contains all"
            " needed model\n");
    fprintf(stdout,
            "       parameters as well as model option flags, and the names"
            " and\n");
    fprintf(stdout, "       locations of all other files.\n");
}

/******************************************************************************
 * @brief    This routine prints out the model Version
 *****************************************************************************/
void
print_version(char *driver)
{
    fprintf(stdout, "  VIC Version : %s\n", SHORT_VERSION);
    fprintf(stdout, "  VIC Driver  : %s\n", driver);

    print_license();
}

/******************************************************************************
 * @brief    This routine prints out license information
 *****************************************************************************/
void
print_license()
{
    fprintf(stdout,
            "\n  Variable Infiltration Capacity (VIC) macroscale hydrologic\n");
    fprintf(stdout,
            "  model version %s, Copyright (C) 2014 Land Surface\n",
            SHORT_VERSION);
    fprintf(stdout,
            "  HydrologyGroup, Dept. of Civil and Environmental Engineering,\n");
    fprintf(stdout,
            "  University of Washington.  VIC comes with ABSOLUTELY NO\n");
    fprintf(stdout,
            "  WARRANTY. This is free software, you may redistribute it\n");
    fprintf(stdout,
            "  under certain conditions; see LICENSE.txt for details.\n\n");

    fprintf(stdout,
            "  Report Bugs and Issues to : https://github.com/UW-Hydro/VIC/issues\n");
    fprintf(stdout,
            "  VIC Users Email Listserve : vic_users@u.washington.edu \n\n");
}
