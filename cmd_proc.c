#include <stdio.h>
#include <stdlib.h>
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

**********************************************************************/
{
  extern option_struct options;
#if LINK_DEBUG
  extern debug_struct debug;
#endif
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
  fprintf(stderr,"Usage: %s -g<model_control_file>\n",temp);
  fprintf(stderr,"\t<model_control_file> is a file that contains all needed model\n\t\tparameters as well as model option flags, and the names and\n\t\tlocations of all other files,\n");
}
