/*     File:          append_temp.c
     Created:       9.26.2000 by EdM, revised 10/29/2000.

       This program reads the output from the script preproc_append.scr
       and formats the daily min or max temperature so the regrid program can read them
       Only the output files from the preproc_append.scr script (daily data)
       and station info files are needed.
      
     *****************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define MAXLINE 4000
#ifndef FALSE
#define FALSE 0
#define TRUE !FALSE
#endif

int isaleap (int iyr);
int julday( int imonth, int ileap );
float ***get_3d_mem_float(int dim1, int dim2, int dim3);
float **get_2d_mem_float(int dim1, int dim2);
char **get_2d_mem_char(int dim1, int dim2);
int *get_1d_mem_int(int dim1);

int main(int argc, char **argv)
{
  FILE *fpinfo, *fpappf, *fpoutf;

  int i, lines, apprecs,year,LP,nyrs;
  int st_yr, end_yr, st_mo, end_mo, st_day, end_day;
  int idx,rec,k;
  int coop_id, append_id;

  char str1[MAXLINE],infofile[80],cid[6];
  char appfile[80], outfile[80];
  char junkstr[MAXLINE];
  char **apid;
  float ***temp;
  float **junk;
  int *apyr;

  int match;
      
  /*     leap year always included in input data */
      
  printf("\nBe sure that the void number used is -99\n");

  /*     read command line arguments -- filenames
         1) station info file from CD data
         2) daily data (preprocessed) to be appended to record
         3) start year for data to be appended
         4) start month for data to be appended
         5) end year for data to be appended
         6) end month for data to be appended
         7) output file of data to be appended
         8) number of records (years * stations) in .daily file
  */
  if (argc!=9) {    /* Must be exactly 8 arguments behind the program name */
    printf("Incorrect number of commandline arguments \n");
    exit(EXIT_FAILURE);
  }

  strcpy(infofile,argv[1]);      printf("infofile %s \n",infofile);
  strcpy(appfile,argv[2]);       printf("appfile %s \n",appfile);
  st_yr=atoi(argv[3]);         printf("st_yr %d \n",st_yr);
  st_mo=atoi(argv[4]);         printf("st_mo %d \n",st_mo);
  end_yr=atoi(argv[5]);         printf("end_yr %d \n",end_yr);
  end_mo=atoi(argv[6]);         printf("end_mo %d \n",end_mo);
  strcpy(outfile,argv[7]);        printf("outfile %s \n",outfile);
  apprecs=atoi(argv[8]);         printf("apprecs %d \n",apprecs);

  if((fpinfo = fopen(infofile,"r"))==NULL){
    printf("Cannot open file %s \n",infofile);exit(0);}
  if((fpappf = fopen(appfile,"r"))==NULL){
    printf("Cannot open file %s \n",appfile);exit(0);}
  if((fpoutf = fopen(outfile,"w"))==NULL){
    printf("Cannot open file %s \n",outfile);exit(0);}
  
  fgets(str1,MAXLINE,fpinfo);
  sscanf(str1,"%d",&lines);
  printf("No. of stations in info file: %d\n",lines);
  nyrs=end_yr-st_yr+1;

  /* allocate memory */
  temp=get_3d_mem_float(nyrs+1,lines+1,366+1);
  junk=get_2d_mem_float(apprecs+1,366+1);
  apid=get_2d_mem_char(apprecs+1,6);
  apyr=get_1d_mem_int(apprecs+1);

  /***read in ".daily" preprocessed append data file into memory***/

  for(rec=1;rec<=apprecs;rec++){
    fscanf(fpappf,"%s %d", apid[rec], &apyr[rec]);
    for(k=1;k<=366;k++) fscanf(fpappf,"%f", &junk[rec][k]);
  }

  /************* begin loop for each station *******************/

  for(i=1;i<=lines;i++) {
    fgets(str1,MAXLINE,fpinfo);
    strncpy(cid,&str1[55],6);
    printf("coopid= %6s ... ",cid);

    /* begin loop for current year */

    for(year=st_yr;year<=end_yr;year++){
      idx=year-st_yr+1;
      /* search append data for current station and year */
      match=FALSE;
      for(rec=1;rec<=apprecs;rec++){
	if(! match) {
	  if ((strcmp(apid[rec],cid)==0) && apyr[rec]==year ){
	    for(k=1;k<=366;k++){
	      temp[idx][i][k]=junk[rec][k];
	    }
	    match=TRUE;
	  }
	}
      }
      if (! match) {
	for(k=1;k<=366;k++){
	  temp[idx][i][k]=-99;
	}
      }
    }
    printf("completed station %d of %d\n",i,lines);
  }

  /************* end loop for each station *******************/
  /* now to write it out in the same format as the existing .fmt file */
  /* since array always is 1-366, pass 1 instead of LP to always
     adjust julian days */

  st_day = julday(st_mo-1,1);
  end_day = julday(end_mo,1)-1;

  for(year=st_yr;year<=end_yr;year++){
    idx=year-st_yr+1;
    LP = isaleap(year);

    printf("checking data for possible bad values and writing data\n");
    for(k=1;k<=366;k++){
      for(i=1;i<=lines;i++){
 	if(temp[idx][i][k] != -99) {
 	  temp[idx][i][k] = 5.0*(temp[idx][i][k]-32.0)/9.0;
 	if(temp[idx][i][k] > 60.0 || temp[idx][i][k] < -60.0)  {
 	  printf("extreme temp at sta %d year %d line %d\n",i,year,k);
 	  printf("temp (C): %.2f\n", temp[idx][i][k]);
	}
	}
      }
      /* write out the daily values for all stations one day at a time
	 except for feb 29th - unless leap year */
      if(LP != 0 || k != 60) {
	if((year > st_yr || k >= st_day) &&
	   (year < end_yr || k <= end_day)) {
	  for(i=1;i<=lines;i++){
 	    if(temp[idx][i][k] <= -99.00) temp[idx][i][k] = -99.00;
 	    fprintf(fpoutf,"%7.2f ",temp[idx][i][k]);
	  }
	  fprintf(fpoutf,"\n");
	}
      }
    }

  }
  printf("closing files\n");
  fclose(fpinfo);
  fclose(fpappf);
  fclose(fpoutf);
  return(0);
}

/****************************************************************
*                      FUNCTIONS                               *
****************************************************************/
/*     returns 1 for leap year, 0 otherwise */
int isaleap (int iyr)
{
  int leap;
  /* return 1 if a leap yr else 0 */

  if( ((!(iyr%4)) && (iyr%100) ) || (!(iyr%400)) )
    leap = 1;
  else
    leap = 0;
  return leap;
}
/****************************************************************/
/*     returns julian day of first day of month passed to function */
/*     for month 13 it returns 366 */
int julday( int imonth, int ileap )
{
  int juldy;
  int mnth[]= {0,31,59,90,120,151,181,212,243,273,304,334,365};

  juldy=mnth[imonth]+1;
  if(imonth >= 2) juldy=juldy+ileap;

  return juldy;
}
/***********************************************************/
float ***get_3d_mem_float(int dim1, int dim2, int dim3)
{
  float ***array_3d;
  int k,j;

  if(!(array_3d = (float***) calloc(dim1,sizeof(float**)))) {
    printf("Cannot allocate memory in get_3d_mem_float\n");
    exit(1); }
  for(k=0; k<dim1;k++) {
    if(!(array_3d[k] = (float**) calloc(dim2,sizeof(float*)))) {
      printf("Cannot allocate memory in get_3d_mem_float\n");
      exit(1); }
    for(j=0; j<dim2;j++) {
      if(!(array_3d[k][j] = (float*) calloc(dim3,sizeof(float)))) {
        printf("Cannot allocate memory in get_3d_mem_float\n");
        exit(1); }
    }
  }
  return array_3d;
}
/***********************************************************/
float **get_2d_mem_float(int dim1, int dim2)
{
  float **array_2d;
  int k;

  if(!(array_2d = (float**) calloc(dim1,sizeof(float*)))) {
    printf("Cannot allocate memory in get_2d_mem_float\n");
    exit(1); }
  for(k=0; k<dim1;k++) {
    if(!(array_2d[k] = (float*) calloc(dim2,sizeof(float)))) {
      printf("Cannot allocate memory in get_2d_mem_float\n");
      exit(1); }
  }
  return array_2d;
}
/***********************************************************/
/***********************************************************/
char **get_2d_mem_char(int dim1, int dim2)
{
  char **array_2d;
  int k;

  if(!(array_2d = (char**) calloc(dim1,sizeof(char*)))) {
    printf("Cannot allocate memory in get_2d_mem_char\n");
    exit(1); }
  for(k=0; k<dim1;k++) {
    if(!(array_2d[k] = (char*) calloc(dim2,sizeof(char)))) {
      printf("Cannot allocate memory in get_2d_mem_char\n");
      exit(1); }
  }
  return array_2d;
}
/***********************************************************/
int *get_1d_mem_int(int dim1)
{
  int *array_1d;

  if(!(array_1d = (int*) calloc(dim1,sizeof(int)))) {
    printf("Cannot allocate memory in get_1d_mem_int\n");
    exit(1); }
  return array_1d;
}
