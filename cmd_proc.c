#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vicNl.h>

static char vcid[] = "$Id$";

filenames_struct cmd_proc(int argc, char *argv[]) 
/**********************************************************************
  cmd_proc                  Keith Cherkauer                1997

  This routine checks the command line for valid program options.  If
  no options are found, or an invalid combination of them appear, the
  routine calls usage() to print the model usage to the screen, before
  exiting execution.

  Modifications:
  11-18-98  Added comment block to cmd_proc() and fixed routine so
            that it will exit if global command file is not defined
            using the "-g" flag.                                KAC
  2003-Oct-03 Added -v option to display version information.		TJB
  2012-Jan-16 Removed LINK_DEBUG code					BN
**********************************************************************/
{
  extern option_struct options;
  extern int getopt();
  extern char *optarg;
  extern char *optstring;

  filenames_struct names;
  int              optchar;
  char             GLOBAL_SET;
  
  if(argc==1) {
    usage(argv[0]);
    exit(1);
  }
  
  GLOBAL_SET = FALSE;

  while((optchar = getopt(argc, argv, optstring)) != EOF) {
    switch((char)optchar) {
    case 'v':
      /** Version information **/
      display_current_settings(DISP_VERSION,(filenames_struct*)NULL,(global_param_struct*)NULL);
      exit(0);
      break;
    case 'o':
      /** Compile-time options information **/
      display_current_settings(DISP_COMPILE_TIME,(filenames_struct*)NULL,(global_param_struct*)NULL);
      exit(0);
      break;
    case 'g':
      /** Global Parameters File **/
      strcpy(names.global, optarg);
      GLOBAL_SET = TRUE;
      break;
    default:
      /** Print Usage if Invalid Command Line Arguments **/
      usage(argv[0]);
      exit(1);
      break;
    }
  }

  if(!GLOBAL_SET) {
    fprintf(stderr,"ERROR: Must set global control file using the '-g' flag\n");
    usage(argv[0]);
    exit(1);
  }

  return names;
}


void usage(char *temp)
/**********************************************************************
	usage		Keith Cherkauer		May 27, 1996

  This routine prints out usage details.

**********************************************************************/
{
  fprintf(stderr,"Usage: %s [-v | -o | -g<global_parameter_file>]\n",temp);
  fprintf(stderr,"  v: display version information\n");
  fprintf(stderr,"  o: display compile-time options settings (set in user_def.h)\n");
  fprintf(stderr,"  g: read model parameters from <global_parameter_file>.\n");
  fprintf(stderr,"       <global_parameter_file> is a file that contains all needed model\n");
  fprintf(stderr,"       parameters as well as model option flags, and the names and\n");
  fprintf(stderr,"       locations of all other files.\n");
}
