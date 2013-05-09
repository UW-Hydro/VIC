#include <stdio.h>
#include <math.h>

main(int argc, char *argv[]) {
/****************************************************************
  sum_basin_param.c        Keith Cherkauer       August 14, 1998

  This program was written to read a runfile mask and computes
  the sum and average parameter values for the activated grid 
  cells.

****************************************************************/

  FILE *frun, *fparam;
  int   i, j, ncol, nrow;
  int   runflag;
  int   N;
  char  tmpstr[512];
  float factor;
  float param;
  float sum;

  if(argc!=4) {
    fprintf(stderr,"Usage: %s <run mask> <param file> <factor>\n",argv[0]);
    fprintf(stderr,"\tThis program sums the parameters divided by factor in the activated grid cells, and prints the results to stdout.  Activated cells have values other than 0 in the run mask file.\n");
    exit(0);
  }

  if((frun=fopen(argv[1],"r"))==NULL) {
    fprintf(stderr,"ERROR: Unable to open %s\n",argv[1]);
    exit(0);
  }

  if((fparam=fopen(argv[2],"r"))==NULL) {
    fprintf(stderr,"ERROR: Unable to open %s\n",argv[2]);
    exit(0);
  }

  factor = atof(argv[3]);

  /** Skip and copy header **/
  for(i=0;i<6;i++) {
    fgets(tmpstr,512,frun);
    fgets(tmpstr,512,fparam);
    if(i==0) sscanf(tmpstr,"%*s %i",&ncol);
    if(i==1) sscanf(tmpstr,"%*s %i",&nrow);
  }

  sum = 0.;
  N   = 0;
  for(i=0;i<nrow;i++) {
    for(j=0;j<ncol;j++) {
      fscanf(frun,"%i",&runflag);
      fscanf(fparam,"%f",&param);
      if(runflag!=0) {
	sum += param / factor;
	N ++;
      }
    }
  }

  fprintf(stdout,"Parameter Sum = %f\nParameter Average = %f\n",
    sum,sum/(float)N);

}

