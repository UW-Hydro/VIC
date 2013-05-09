#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/** Define the units conversion factor between the simulated 
    discharge (cfs), and the observed dischare (???) **/
#define sim_to_obs 1

/** Define the conversion factor from observed simulation rate to total
  volume **/


main(int argc, char *argv[]) {
/************************************************************************
  compute_R2.c              Keith Cherkauer            January 10, 1999

  This program was written to comute the R^2 value for the current 
  simulated daily discharge.

************************************************************************/

  FILE *fvic, *fobs, *fout;
  int startrec, endrec;
  int Nrecs, i;
  float basinfactor;
  float meanobs;
  float meanvic;
  float Rsqr;
  float tmpdata;
  float *vicdata;
  float *obsdata;
  float vicflow, obsflow;
  float numsum, densum;

  if(argc!=7) {
    fprintf(stderr,"Usage: %s <VIC discharge> <observed discharge> <start rec> <end rec> <basin size factor> <R2 outfile>\n",argv[0]);
    exit(0);
  }

  if((fvic=fopen(argv[1],"r"))==NULL) {
    fprintf(stderr,"ERROR: Unable to open VIC simulation file %s.\n",argv[1]);
    exit(1);
  }

  if((fobs=fopen(argv[2],"r"))==NULL) {
    fprintf(stderr,"ERROR: Unable to open observed flow file %s.\n",argv[2]);
    exit(2);
  }

  if((fout=fopen(argv[6],"w"))==NULL) {
    fprintf(stderr,"ERROR: Unable to open output file %s.\n","R2.out");
    exit(1);
  }

  startrec = atoi(argv[3]);
  endrec = atoi(argv[4]);
  basinfactor = atof(argv[5]);

  /** Process Data **/
  Nrecs   = endrec - startrec + 1;
  vicdata = (float*)calloc(Nrecs,sizeof(float));
  obsdata = (float*)calloc(Nrecs,sizeof(float));

  for(i=0;i<startrec;i++) {
    fscanf(fvic,"%*s %f\n",&tmpdata);
    fscanf(fobs,"%*s %f\n",&tmpdata);
  }

  meanvic = meanobs = 0;
  vicflow = obsflow = 0;
  for(i=0;i<Nrecs;i++) {
    fscanf(fvic,"%*s %f\n",&vicdata[i]);
    fscanf(fobs,"%*s %f\n",&obsdata[i]);

    /* Convert Simulated Flow to Match Observed */
    vicdata[i] *= sim_to_obs;

    vicdata[i] *= basinfactor;
    meanvic += vicdata[i];
    meanobs += obsdata[i];
    /* Convert discharge from cubic feet per sec to cubic meters 
       (per time step) to look at total flow volumes */
    vicflow += vicdata[i] * 0.02831685 * 24. * 60. * 60.;
    obsflow += obsdata[i] * 0.02831685 * 24. * 60. * 60.;
  }
  meanvic /= (float)Nrecs;
  meanobs /= (float)Nrecs;
  vicflow /= 1.e9;
  obsflow /= 1.e9;

  numsum = densum = 0;
  for(i=0;i<Nrecs;i++) {
    numsum += pow( obsdata[i] - vicdata[i], 2. );
    densum += pow( obsdata[i] - meanobs, 2. );
  }
  Rsqr = -1. * ( 1. - ( numsum / densum) );

  fprintf(fout,"Total observed flow  = %f km^3.\n",obsflow);
  fprintf(fout,"Total simulated flow = %f km^3.\n",vicflow);
  fprintf(fout,"Mean observed flow   = %f cfs.\n",meanobs);
  fprintf(fout,"Mean simulated flow  = %f cfs.\n",meanvic);
  fprintf(fout,"R^2 model error      = %f.\n",Rsqr);

  fprintf(stdout,"%f\n", Rsqr);

  free((char*)vicdata);
  free((char*)obsdata);

  fclose(fvic);
  fclose(fobs);
  fclose(fout);

}
