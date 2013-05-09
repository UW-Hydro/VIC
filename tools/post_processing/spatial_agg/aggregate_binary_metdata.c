/*
 * SUMMARY:      Program used to aggregate metdata to a coarser resolution. 
 * USAGE:        Needs the lat,lon,value from fraction file of the 2 resolutions. 
 *               It also needs the inputdirectory of the old metdata.    
 *
 * AUTHOR:       Bernt Viggo Matheussen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       tempbvm@hydro.washington.edu
 * ORIG-DATE:     3-May-99 at 15:14:18
 * LAST-MOD: Fri Jul  2 14:14:37 1999 by Keith Aric Cherkauer <cherkaue@u.washington.edu>
 * DESCRIPTION:  
 * DESCRIP-END.
 * FUNCTIONS:    
 * COMMENTS:     
 */

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <float.h>
#include <math.h>
#include <string.h>

#define MAX 5000
#define STEPS 17897
#define MAX_PARAM 6
#define OLD_DECIMAL_PLACES 4   /* Number of decimal places used by 
				  old gridded met station files */
#define NEW_DECIMAL_PLACES 4   /* Number of decimal places used by 
				  new gridded met station files */

/* Counts number of lines in file */
int line_count(char file[400]);

/* Reads in data */
void read_data(char file[400], float FR4[MAX][3], int cells);

/* Check to see if file exists */
float exist_file(float FR8[MAX][3],float lat, float lon, int cells);

/* Reads in metdata to array */
void read_metdata(char indir[400], float AFR[STEPS][MAX_PARAM], float lat, float lon, int Nparam);

void main(int argc, char *argv[]) {
  FILE *fp;
  char fr4[400];
  char fr8[400];
  char indir[400]; 
  char outdir[400];
  char outfile[400];
  char LATLON[50];
  char fmtstr[512];

  int i,j;
  int cells8;
  int cells4;
  int param, Nparam;
 
  float FR4[MAX][3];
  float FR8[MAX][3];

  float AFR[STEPS][MAX_PARAM];
  float BFR[STEPS][MAX_PARAM];
  float CFR[STEPS][MAX_PARAM];
  float DFR[STEPS][MAX_PARAM];
  float NEW[STEPS][MAX_PARAM];

  float res8;
  float res4;
  float lat,lon;
  float Alat,Blat,Clat,Dlat; 
  float Alon,Blon,Clon,Dlon; 

  float Afr,Bfr,Cfr,Dfr;
  int Aok,Bok,Cok,Dok;
  float sum_fraction;  
  int valid_cells;

  short int tmpval[6];
  
  if(argc!= 7) {
    fprintf(stderr,"Usage: %s <old fraction xyz> <old fraction data dir> <new fraction xyz> <new fraction data dir> <old resolution> <number of parameters>\n",argv[0]);
    fprintf(stderr,"\t<old fraction xyz> is the 3 coloumn fraction file for the old resolution.\n");
    fprintf(stderr,"\t<old fraction data dir> is the directory for the forcing data for the old resolution.\n");
    fprintf(stderr,"\t<new fraction xyz> is the 3 coloumn fraction file for the new resolution.\n");
    fprintf(stderr,"\t<new fraction data dir> is the directory for the forcing data for the new resolution.\n");
    fprintf(stderr,"\t<old resolution> is the old resoluation in fractions of degrees.\n");
    fprintf(stderr,"\t<number of parameters> is the number of parameters (data types / columns) in the forcing file.\n");
    exit(0);
  }

  strcpy(fr8,argv[1]);
  strcpy(indir,argv[2]); 
  strcpy(fr4,argv[3]);
  strcpy(outdir,argv[4]);

  res8 = atof(argv[5]); /* Remember to change this one */

  res4 = 2*res8;

  Nparam = atoi(argv[6]);

  cells4 = line_count(fr4);
  cells8 = line_count(fr8);

  printf("File %s has %d lines\n",fr4,cells4);
  printf("File %s has %d lines\n",fr8,cells8);
  
  read_data(fr4,FR4,cells4);
  read_data(fr8,FR8,cells8);
  
  for (i=0;i<cells4;i++) {  
    for(j=0;j<STEPS;j++) { /* Setting initial values */
      for(param=0;param<Nparam;param++) {
	AFR[j][param] = 0.0;
	BFR[j][param] = 0.0;
	CFR[j][param] = 0.0;
	DFR[j][param] = 0.0;
      }
    }
    
    lat      = FR4[i][0];
    lon      = FR4[i][1];
    
    Alat = lat + res8/2;
    Alon = lon - res8/2;
    Blat = lat + res8/2;
    Blon = lon + res8/2;
    Clat = lat - res8/2;
    Clon = lon - res8/2;
    Dlat = lat - res8/2;
    Dlon = lon + res8/2;
    
    Afr = exist_file(FR8,Alat,Alon,cells8);
    Bfr = exist_file(FR8,Blat,Blon,cells8);
    Cfr = exist_file(FR8,Clat,Clon,cells8);
    Dfr = exist_file(FR8,Dlat,Dlon,cells8);
    
    Aok = Bok = Cok = Dok = 0;
    
    if(Afr > 0.0) Aok = 1;
    if(Bfr > 0.0) Bok = 1;
    if(Cfr > 0.0) Cok = 1;
    if(Dfr > 0.0) Dok = 1;
    
    sum_fraction = Afr + Bfr + Cfr + Dfr;
    valid_cells = Aok + Bok + Cok + Dok;
    
    sprintf(fmtstr,"lat  %%.%if lon  %%.%if valid_cells %%d sum_fraction %%.3f\n",OLD_DECIMAL_PLACES,OLD_DECIMAL_PLACES);
    printf(fmtstr,lat,lon,valid_cells,sum_fraction); 
    printf("Afr %.3f Bfr %.3f Cfr %.3f Dfr %.3f\n",Afr,Bfr,Cfr,Dfr);
    
    /* Reading files */
    if(Afr > 0) read_metdata(indir,AFR,Alat,Alon,Nparam);
    if(Bfr > 0) read_metdata(indir,BFR,Blat,Blon,Nparam);
    if(Cfr > 0) read_metdata(indir,CFR,Clat,Clon,Nparam);
    if(Dfr > 0) read_metdata(indir,DFR,Dlat,Dlon,Nparam);
    
    
    for(j=0;j<STEPS;j++) { /* Calculating new values */
      for(param=0;param<Nparam;param++) {
	NEW[j][param] = (AFR[j][param]*Afr + BFR[j][param]*Bfr 
			 + CFR[j][param]*Cfr + DFR[j][param]*Dfr)
	  / sum_fraction;
      }
    }
    
    for(param=0;param<Nparam;param++)
      printf("%.4f",NEW[0][param]);
    printf("\n");
    
    strcpy(outfile,outdir);
    sprintf(fmtstr,"data_%%.%if_%%.%if",NEW_DECIMAL_PLACES,NEW_DECIMAL_PLACES);
    sprintf(LATLON,fmtstr,lat,lon);
    strcat(outfile,LATLON);
    printf("outfile %s\n",outfile);
    
    if((fp = fopen(outfile,"w"))==NULL) {
      printf("Cannot open file %s \n",outfile);exit(0);
    }   
    
    for(j=0;j<STEPS;j++) {
      for(param=0;param<Nparam;param++) {
	tmpval[param] = (short int)rint(NEW[j][param]);
      }
      fwrite(tmpval,6,sizeof(short int),fp);
    }
    
    fclose(fp); 
  }

}/* END MAIN ******************************************************************/


/* Function returns number of lines in file */
int line_count(char file[400])
{
  FILE *fp;
  int c, lines;
  if((fp = fopen(file,"r"))==NULL){
    printf("Cannot open file %s \n",file);exit(0);}
  lines = 0;
  while((c = fgetc(fp)) !=EOF)  if (c=='\n') lines++;
  fclose(fp);
  return lines;
}/* END function int line_count(char file[200])   */



void read_data(char file[400], float FR4[MAX][3], int cells) {
  FILE *fp;
  int i;

  if((fp = fopen(file,"r"))==NULL){
    printf("Cannot open file %s \n",file);exit(0);
  }
  for(i=0;i<cells;i++) fscanf(fp,"%f %f %f",&FR4[i][0],&FR4[i][1],&FR4[i][2]);
  fclose(fp);
}
/* END function void read_data() */


float exist_file(float FR8[MAX][3],float lat, float lon, int cells) {
  int i, factor;
  float return_fraction=0.0;

  factor = pow(10,OLD_DECIMAL_PLACES-1);

  for(i=0;i<cells;i++) {
    if(  ((long)(factor*lat) == (long)(factor*FR8[i][0])) &&  
         ((long)(factor*lon) == (long)(factor*FR8[i][1])) ) {
      return_fraction = FR8[i][2];
    }
  }
  return return_fraction;
}/* END function */ 




void read_metdata(char indir[400], 
		  float AFR[STEPS][MAX_PARAM], 
		  float lat, 
		  float lon, 
		  int Nparam)
{
  FILE *fp;
  int i, param;
  char LATLON[50];
  char file[400];
  char fmtstr[15];
  short int tmpval[6];

  /* Make filename */
  sprintf(fmtstr,"data_%%.%if_%%.%if",OLD_DECIMAL_PLACES,OLD_DECIMAL_PLACES);
  strcpy(file,indir);
  sprintf(LATLON,fmtstr,lat,lon);
  strcat(file,LATLON);
  printf("%s\n",file);
  if((fp = fopen(file,"r"))==NULL) { 
    printf("Cannot open file %s \n",file);exit(0);
  } 
  for(i=0;i<STEPS;i++) {
    fread(tmpval,6,sizeof(short int),fp);
    for(param=0;param<Nparam;param++)
      AFR[i][param] = tmpval[param];
  }
  fclose(fp); 
}/* ENDE function */

