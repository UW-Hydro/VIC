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
    fprintf(stderr,
            "Usage: %s [-v | -o | -g <global_parameter_file>]\n", executable);
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

/******************************************************************************
 * @brief    This routine prints out the model Version
 *****************************************************************************/
void print_version(char *driver)
{
    fprintf(stderr, "  VIC Version : %s\n", SHORT_VERSION);
    fprintf(stderr, "  VIC Driver  : %s\n", driver);

    print_license();

}

/******************************************************************************
 * @brief    This routine prints out license information
 *****************************************************************************/
void print_license()
{

 fprintf(stderr,
         "\n  Variable Infiltration Capacity (VIC) macroscale hydrologic\n");
 fprintf(stderr,
         "  model version %s, Copyright (C) 2014 Land Surface\n",
         SHORT_VERSION);
 fprintf(stderr,
         "  HydrologyGroup, Dept. of Civil and Environmental Engineering,\n");
 fprintf(stderr,
         "  University of Washington.  VIC comes with ABSOLUTELY NO\n");
 fprintf(stderr,
         "  WARRANTY. This is free software, you may redistribute it\n");
 fprintf(stderr,
         "  under certain conditions; see LICENSE.txt for details.\n\n");

 fprintf(stderr,
         "  Report Bugs to: https://github.com/UW-Hydro/VIC/issues\n");
 fprintf(stderr,
         "  VIC Users Email List:  vicadmin@hydro.washington.edu\n\n");
}
