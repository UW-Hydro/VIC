/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine checks the command line for valid program options.
 *****************************************************************************/

#include <unistd.h>
#include <getopt.h>

#include <vic_driver_shared_all.h>

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
            exit(EXIT_SUCCESS);
            break;
        case 'o':
            /** Compile-time options information **/
            display_current_settings(DISP_COMPILE_TIME);
            exit(EXIT_SUCCESS);
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
        exit(EXIT_FAILURE);
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
    fprintf(stdout, "VIC Driver  : %s\n", driver);
    fprintf(stdout, "VIC Version : %s\n", VERSION);
    fprintf(stdout, "VIC Git Tag : %s\n", GIT_VERSION);
    fprintf(stdout, "Compiled    : by %s on %s (%s) %s %s\n",
            USERNAME, HOSTNAME, PLATFORM, BUILD_DATE, BUILD_TIME);
    fprintf(stdout, "Compiler    : %s\n", COMPILER);
    fprintf(stdout, " version    : %s\n", COMPILER_VERSION);

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
            "  model, version %s\n\n", SHORT_VERSION);
    fprintf(stdout,
            "  MIT LICENSE\n\n");
    fprintf(stdout,
            "  Copyright (C) 2021 VIC Developers\n\n");
    fprintf(stdout,
            "  See LICENSE.txt for details.\n\n");
    fprintf(stdout,
            "  Report Bugs and Issues to : https://github.com/UW-Hydro/VIC/issues\n");
    fprintf(stdout,
            "  VIC Users Email Listserve : vic_users@u.washington.edu \n\n");
}
