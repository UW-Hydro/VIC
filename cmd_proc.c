#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

filenames_struct cmd_proc(int argc, char *argv[]) 
{
  extern option_struct options;
  extern debug_struct debug;

  filenames_struct names;
  extern int getopt();
  extern char *optarg;
  extern char *optstring;
  int optchar;
  
  if(argc==1) {
    usage(argv[0]);
    exit(1);
  }
  
  while((optchar = getopt(argc, argv, optstring)) != EOF) {
    switch((char)optchar) {
    case 'g':
      /** Global Parameters File **/
      strcpy(names.global, optarg);
      break;
    case 'I':
      /** Soil Initialization File **/
      strcpy(names.init_soil, optarg);
      options.INIT_SOIL=TRUE;
      break;
    case 'S':
      /** Snow Initialization File **/
      strcpy(names.init_snow, optarg);
      options.INIT_SNOW=TRUE;
      break;
    default:
      /** Print Usage if Invalid Command Line Arguments **/
      usage(argv[0]);
      exit(1);
      break;
    }
  }
  return names;
}


void usage(char *temp)
/**********************************************************************
	usage		Keith Cherkauer		May 27, 1996

  This routine prints out usage details.

**********************************************************************/
{
  fprintf(stderr,"Usage: %s -g<model_control_file> ",temp);
  fprintf(stderr,"[-I<soil_init_file>] [-S<snow_init_file>]\n");
  fprintf(stderr,"\t<model_control_file> is a file that contains all needed model\n\t\tparameters as well as model option flags, and the names and\n\t\tlocations of all other files,\n");
  fprintf(stderr,"\t<soil_init_file> is the initialization file for soil thermal and\n\t\tmoisture profiles,\n");
  fprintf(stderr,"\t<soil_init_file> is the initialization file for the snowpack\n\t\t[no snow assumed],\n");
}
