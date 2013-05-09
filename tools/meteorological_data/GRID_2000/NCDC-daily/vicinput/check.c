/* File:             check.c                                                       */
/* Programmer :      Bernt Viggo Matheussen                                         */
/* Date:             02.24.98                                                       */
/* Version:          1.0                                                            */

/* This program filters the final regridded data checking for                       */
/* maximum values, minimum values, monthly means, yearly mean                       */
/* Input is one of the datafiles data_42.88_-119.88                                 */
/* Output is a file containing max, min, SD, etc                                    */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#define COLS (3)       /* Number of coloumns in inputfile */

int nr_timesteps(char pp[400]);  /* Function reurns number of timesteps in inputfile */
void one_mnth(FILE *fp1, int yy, int mn, float month[30]);
void main(int argc, char **argv)     /* commandline arguments */
{
  FILE *fpin, *fpout;       /* Filepointers to in and output files */
  char infile[400], startyear[6],outfile[400]; /* Strings to commandline arguments */
  int syear, run_years;
  int steps, month;
  float mnth_prcp[600], mnth_tmax[600], mnth_tmin[600];
  float mnth[30];

  /* mnth[0] = monthly sum prcp     mnth[1]= monthly avg  tmax            */
  /* mnth[2] = monthly avg tmin     mnth[3]= monthly min  prcp            */
  /* mnth[4] = monthly min tmax     mnth[5]= monthly min tmin             */
  /* mnth[6] = monthly max prcp     mnth[7]= monthly max tmax             */
  /* mnth[8] = monthly max tmin     mnth[9]= monthly SD prcp              */
  /* mnth[10] = monthly SD tmax     mnth[11]= monthly SD tmin             */
  int i,j;

  if (argc<4) {           /* Must be exactly 3 arguments behind the program name */
  printf("Not correct number of commandline arguments \n");
  printf("filter \"inputfile\" \"startyear\" \"outputfile\"  \n"); 
  exit(EXIT_FAILURE); }
  strcpy(infile,argv[1]);    printf("%s \n",infile);
  strcpy(startyear,argv[2]); printf("%s \n",startyear); syear = atoi(startyear);
  strcpy(outfile,argv[3]);   printf("%s \n",outfile);
  steps = nr_timesteps(infile);  printf("steps = %d\n",steps);
  if((fpin = fopen(infile,"r"))==NULL){ 
    printf("Cannot open file %s \n",infile);exit(0);}
  if((fpout = fopen(outfile,"w"))==NULL){ 
    printf("Cannot open file %s \n",outfile);exit(0);}

  run_years = (int)(steps/365);
  printf("run_years %d\n",run_years);

  fprintf(fpout,"    sump   avgtmax   avgtmin minprcp   mintmax   mintmin   maxtprcp   maxtmax   maxtmin SDprcp   SDtmax  SDtmin \n");
  for (j=syear;j<(syear+run_years);j++)
    {
    for (i=0;i<12;i++)
      {
      one_mnth(fpin,j,i,mnth);
      fprintf(fpout,"%8.2f  %8.2f  %8.2f",mnth[0],mnth[1],mnth[2]);
      fprintf(fpout,"%8.2f  %8.2f  %8.2f",mnth[3],mnth[4],mnth[5]);
      fprintf(fpout,"%8.2f  %8.2f  %8.2f",mnth[6],mnth[7],mnth[8]);
      fprintf(fpout,"%8.2f  %8.2f  %8.2f\n",mnth[9],mnth[10],mnth[11]);
      }
    }
  fclose(fpin);
  fclose(fpout);
}/*** END main ************************************************************************/
/**************************************************************************************/

int nr_timesteps(char pp[400]) /* Returns number of timesteps in inputfile */
{
  FILE *fp1;  char str1[20];  long i;  int r;
  printf("%s\n", pp);
  if((fp1 = fopen(pp,"r"))==NULL){       /* opens file  */
    printf("Cannot open file %s \n",pp);  exit(0);}  i =0;
  while (fscanf(fp1,"%s",str1)!=EOF)i++;  r =(int)(i/COLS);  return r;
  fclose(fp1);
}/*** END function int nr_timesteps(char pp[400]) ********/
/*********************************************************/


void one_mnth(FILE *fp1, int yy, int mn, float  month[30])
{
  int steps, i;
  char *regular[] = {"31","28","31","30","31","30","31","31","30","31","30","31"};
  char *leap[]    = {"31","29","31","30","31","30","31","31","30","31","30","31"};
  float cc, min_prcp, min_tmax, min_tmin;
  float max_prcp, max_tmax, max_tmin;
  float prcp[31], tmax[31], tmin[31];
  double mean,SDprcp, SDtmax, SDtmin;
  double VARprcp, VARtmax, VARtmin;

  if (yy%4==0) steps = atoi(leap[mn%12]);  
  if (yy%4!=0) steps = atoi(regular[mn%12]);
  printf("steps %d\n",steps);
  month[0] = month[1]= month[2] =0;
  min_prcp = min_tmax = min_tmin = 1000;
  max_prcp = max_tmax = max_tmin = -1000;
  for (i=0;i<steps;i++)
    {
    fscanf(fp1,"%f",&cc); month[0]=month[0]+cc;  /* Read precipitation */
    if (cc<min_prcp) min_prcp = cc;              /* Checking min value */
    if (cc>max_prcp) max_prcp = cc;              /* Checking max value */ 
    prcp[i]=cc;
    fscanf(fp1,"%f",&cc); month[1]=month[1]+cc;  /* Read Tmax */
    if (cc<min_tmax) min_tmax = cc;              /* Checking min value */
    if (cc>max_tmax) max_tmax = cc;              /* Checking max value */ 
    tmax[i]=cc;
    fscanf(fp1,"%f",&cc); month[2]=month[2]+cc;  /* Read Tmin */
    if (cc<min_tmin) min_tmin = cc;              /* Checking min value */
    if (cc>max_tmin) max_tmin = cc;              /* Checking max value */ 
    tmin[i]=cc;
    }
  month[1]=month[1]/steps;
  month[2]=month[2]/steps;
  month[3]=min_prcp;
  month[4]=min_tmax;
  month[5]=min_tmin;
  month[6]=max_prcp;
  month[7]=max_tmax;
  month[8]=max_tmin;
  SDprcp = SDtmax = SDtmin = VARprcp = VARtmax = VARtmin = 0.00;

  for (i=0;i<steps;i++)      /*  Finding the SD for prcp,tmin, tmax */
    {
    SDprcp= SDprcp + prcp[i];    
    SDtmax= SDtmax + tmax[i];
    SDtmin= SDtmin + tmin[i];
    }
  SDprcp = SDprcp/steps;   /* Mean values */
  SDtmax = SDtmax/steps;
  SDtmin = SDtmin/steps;

  for (i=0;i<steps;i++)   /* VARIANS  */
    {
    VARprcp = VARprcp + ((prcp[i]-SDprcp)*(prcp[i]-SDprcp));    
    VARtmax = VARtmax + ((tmax[i]-SDtmax)*(tmax[i]-SDtmax));    
    VARtmin = VARtmin + ((tmin[i]-SDtmin)*(tmin[i]-SDtmin));    
    }
  VARprcp = VARprcp/(steps-1);
  VARtmax = VARtmax/(steps-1);
  VARtmin = VARtmin/(steps-1);
  SDprcp = sqrt(VARprcp); month[9]  = SDprcp; /* Standard deviation */
  SDtmax = sqrt(VARtmax); month[10] = SDtmax; /* Standard deviation */
  SDtmin = sqrt(VARtmin); month[11] = SDtmin; /* Standard deviation */

}/**** END void one_mnth(FILE *fp1, int yy, int mn) **********************************/
/*************************************************************************************/
