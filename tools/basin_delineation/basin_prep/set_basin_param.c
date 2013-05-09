#include <stdio.h>
#include <math.h>

main(int argc, char *argv[]) {
/****************************************************************
  set_basin_param.c        Keith Cherkauer       August 14, 1998

  This program was written to read a runfile mask and replaces
  parameters in the activated grid cells with the given value.
****************************************************************/

  FILE *frun, *foldparam, *fnewparam;
  int   i, j, ncol, nrow;
  int   runflag;
  char  tmpstr[512];
  float value;
  float param;

  if(argc!=5) {
    fprintf(stderr,"Usage: %s <run mask> <old param file> <new param file> <value>\n",argv[0]);
    fprintf(stderr,"\tThis program sets the activated grid cells in the old parameter file to the given value, and writes the results to the new parameter file.  Activated cells have values other than 0 in the run mask file.\n");
    exit(0);
  }

  if((frun=fopen(argv[1],"r"))==NULL) {
    fprintf(stderr,"ERROR: Unable to open %s\n",argv[1]);
    exit(0);
  }

  if((foldparam=fopen(argv[2],"r"))==NULL) {
    fprintf(stderr,"ERROR: Unable to open %s\n",argv[2]);
    exit(0);
  }

  if((fnewparam=fopen(argv[3],"w"))==NULL) {
    fprintf(stderr,"ERROR: Unable to open %s\n",argv[3]);
    exit(0);
  }

  value = atof(argv[4]);

  /** Skip and copy header **/
  for(i=0;i<6;i++) {
    fgets(tmpstr,512,frun);
    fgets(tmpstr,512,foldparam);
    if(i==0) sscanf(tmpstr,"%*s %i",&ncol);
    if(i==1) sscanf(tmpstr,"%*s %i",&nrow);
    fprintf(fnewparam,"%s",tmpstr);
  }

  for(i=0;i<nrow;i++) {
    for(j=0;j<ncol;j++) {
      fscanf(frun,"%i",&runflag);
      fscanf(foldparam,"%f",&param);
      if(runflag!=0) {
	fprintf(fnewparam,"%f",value);
      }
      else {
	fprintf(fnewparam,"%f",param);
      }
      if(j==ncol-1) fprintf(fnewparam,"\n");
      else fprintf(fnewparam," ");
    }
  }
}

