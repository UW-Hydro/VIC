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
  fprintf(stderr,"Usage: %s -i<prec_temp_path> -s<soil_param>",temp);
  fprintf(stderr," -v<veg_param>\n\t\t-L<veg_lib>");
  fprintf(stderr," -g<global_param> -r<result_path> [-o<snowout_path>]\n");
  fprintf(stderr,"\t\t[-a<sawd_path>] [-I<soil_init_file>] [-S<snow_init_file>]\n");
  fprintf(stderr,"\t\t[-H]\n");
  fprintf(stderr,"where:\tthe <prec_temp_path> is the location of gridded precipitation station data,\n");
  fprintf(stderr,"\t<soil_param> is a file that contains soil parameters for each grid cell in the model,\n");
  fprintf(stderr,"\t<veg_param> is a file that contains vegetation coverage frzctions for each grid cell in the model,\n");
  fprintf(stderr,"\t<veg_lib> is a file that contains parameters for all vegetation types in the model,\n");
  fprintf(stderr,"\t<global_param> is a file that contains global parameters for the model run,\n");
  fprintf(stderr,"\t<result_path> is the location where gridded output files will be written.\n");
  fprintf(stderr,"\t<snowout_path> is the location of gridded output files from the NWS snow melt model (required for water balance model),\n");
  fprintf(stderr,"\t\tif set to INTERNAL, VIC will use its own snow model instead of reading from files.\n");
  fprintf(stderr,"\t<soil_init_file> is the initialization file for soil thermal and moisture profiles,\n");
  fprintf(stderr,"\t<soil_init_file> is the initialization file for the snowpack [no snow assumed],\n");
  fprintf(stderr,"\t-H turns on hourly precipitation in surface airays data file\n");
  fprintf(stderr,"\t-a <sawd_path> is the location of gridded hourly surface airways data files (required for full energy balance model),\n");
  fprintf(stderr,"\n\tNOTE: in the case of all paths, that path also contains the start of the file name, grid cell information is tacked on to the end of the file name by the model code.\n");
}
