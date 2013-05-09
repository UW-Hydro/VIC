/* File:             rescale_airtemp.c                                              */
/* Programmer :      Bernt Viggo Matheussen                                         */
/* Date:             02.27.99                                                       */
/* Version:          1.0                                                            */

/* This program rescales the gridded airtemperature data from the regrid program    */
/* with Tmax = Tmax + TmaxPrism -TmaxMonth                                          */
/* This is done to get the long term mean airtemperatures to be the same as PRISM   */

/* Inputfiles are the gridded data file, monthly values from mk_monthly_airtemp.c   */
/* and get_prism_airtemp.aml                                                        */
/* Output is a new rescaled gridded file that can be used by vicinput.c             */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

int timesteps(int yy, int mt);

void main(int argc, char **argv) 
{ 
  FILE *fpmask, *fpmonthly, *fpprism,*fpin_grd, *fpout_grd;        
  char maskfile[400];                    /* Total length of directory and filename < 400     */
  char monthly[400],old_monthly[400];    /* Using files named monthly.jan, monthly.feb, .... */
  char prism[400],old_prism[400];        /* Using files named prism.jan, prism.feb......     */
  char in_grd[400];       
  char out_grd[400];
  char start_year[400];    /* The first month in the start year has to be january            */
  char *monthstr[]={"jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec"};
  char void_nr[20];
  char str1[20],mm[20],pp[20],grd[20];
  char end_year[400];
  float value, valpp;
  int year,endyy;    /* month number 0 is january and 1 is february */
  int i, j, k, l;      /* Used only as counters in loops */


  if (argc<8) {           /* Must be exactly 7 arguments behind the program name */
  printf("Not correct number of commandline arguments \n");
  exit(EXIT_FAILURE); }
  strcpy(maskfile,argv[1]);      printf("%s \n",maskfile);
  strcpy(monthly,argv[2]);       printf("%s \n",monthly);
  strcpy(prism,argv[3]);         printf("%s \n",prism);
  strcpy(in_grd,argv[4]);        printf("%s \n",in_grd);
  strcpy(out_grd,argv[5]);       printf("%s \n",out_grd);
  strcpy(start_year,argv[6]);    printf("%d \n",atoi(start_year));   year = atoi(start_year);
  strcpy(end_year,argv[7]);      printf("%d \n",atoi(end_year));     endyy = atoi(end_year);
  strcpy(old_monthly,monthly); 
  strcpy(old_prism,prism);  


  if((fpmask = fopen(maskfile,"r"))==NULL){                  /* opens files that will be open all time*/
      printf("Cannot open file %s \n",maskfile);exit(0);}
  if((fpin_grd = fopen(in_grd,"r"))==NULL){
      printf("Cannot open file %s \n",in_grd);exit(0);}
  if((fpout_grd = fopen(out_grd,"w"))==NULL){
      printf("Cannot open file %s \n",out_grd);exit(0);}
  endyy++;
  for(l=year;l<(endyy);l++)     /* Controling how many years to run  */
    {  
    for (k=0;k<12;k++)   /* month number 0 is january and 1 is february */
      {
      strcpy(monthly,old_monthly);strcat(monthly,monthstr[k]); /* open monthly.mnth */
      strcpy(prism,old_prism);strcat(prism,monthstr[k]);       /* open prism.mnth   */

      if((fpmonthly = fopen(monthly,"r"))==NULL){
        printf("Cannot open file %s \n",monthly);exit(0);}
      if((fpprism = fopen(prism,"r"))==NULL){
        printf("Cannot open file %s \n",prism);exit(0);}

      for(j=0;j<timesteps(l,k);j++)     /* Number of timesteps to run */
        { 
        /***** ONE TIMESTEP STARTS HERE  ****************************************/
        for(i=0;i<12;i++)fscanf(fpmonthly,"%s",str1);    /* Moving down to data in monthly.jan, .... */
        for(i=0;i<12;i++){fscanf(fpmask,"%s",str1);strcpy(void_nr,str1);}/*Moving down in maskfile*/
        for(i=0;i<12;i++)fscanf(fpprism,"%s",str1);    /* Moving down to data in prism.jan, .... */
        while(fscanf(fpmask,"%s",str1)!=EOF){ /* Using the mask to decide which gridcells to read */
          
          fscanf(fpmonthly,"%s",mm); 
          fscanf(fpprism,"%s",pp);
          valpp = atof(pp);
          if (atof(pp)<0.01) valpp = atof(mm); /* Checking if value exists. */

          if(atof(str1)>0.001)   /* Asuming all stations have elevation above 0.001 meter  */
            {
            fscanf(fpin_grd,"%s",grd); /* Reading the value from the inputfile (test.grd)  */
            value = atof(grd) + valpp - atof(mm);  /* reascales the airtemperature */
            fprintf(fpout_grd,"%4.2f ",value);
            }
          }
        fprintf(fpout_grd,"\n");
        rewind(fpmask);        /* Must rewind the filepointers for every timestep */
        rewind(fpmonthly);
        rewind(fpprism);
        /***** ONE TIMESTEP ENDS HERE   ******************************************/
        /*printf("Timestep %d\n",j);*/
        }
        fclose(fpmonthly);
        fclose(fpprism);          /* close temporary files */
      }
    printf("year %d \n",l);
    }
 
  fclose(fpout_grd);
  fclose(fpin_grd);
  fclose(fpmask);  
} /** END MAIN  ***/



int timesteps(int yy, int mt)
{
  int val;
  char *regular[]={"31","28","31","30","31","30","31","31","30","31","30","31"};
  char *leap[]=   {"31","29","31","30","31","30","31","31","30","31","30","31"};
  if (yy%4!=0) val = atoi(regular[mt%12]);
  if (yy%4==0) val = atoi(leap[mt%12]);
  return val;
} /*** END of function int timesteps(int yy, int mt) ************************/ 



