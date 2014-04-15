#include <unistd.h>
#include <vic_def.h>
#include <vic_run.h>
#include <vic_driver_image.h>

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

/**********************************************************************
        usage   Keith Cherkauer   May 27, 1996

   This routine prints out usage details.
**********************************************************************/
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
