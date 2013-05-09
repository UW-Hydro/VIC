/* File:             vicinput.c                                              */
/* Programmer :      Bernt Viggo Matheussen                                  */
/* Date:             02.10.98                                                */
/* Version:          1.0                                                     */

/* This program reads the output from A.H.s regridding program and           */
/* generates the vicinput files. The input needed is prcp, tmin,             */
/* tmax, maskfile, outp_directory, binary_flag                               */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#define stp50 (50)   /* Number of timesteps to run before write to files */
#define VALID_CELLS (10000) /* If number of valid cells in maskfile is bigger
			       than this this number must be changed */
void make_flist(char msk[BUFSIZ+1]);  /* Function that makes a filelist */
void make_outfiles(char out_dir[BUFSIZ+1],int binflag); /* function that make outputfiles */
int valid_cells(char maskf[BUFSIZ+1]);  /* returns number of valid cells*/
int nr_timesteps(char pp[BUFSIZ+1],int cc, int binflag);  

/* This function reads one timestep into three arrays --CHANGED TO SHORT INT FOR BINARY */
void one_tstp(FILE *fppp, FILE *fpta, FILE *fpti, 
              short int precip[VALID_CELLS], short int tempmax[VALID_CELLS],
              short int tempmin[VALID_CELLS], int cl);
void one_tstp_ascii(FILE *fppp, FILE *fpta, FILE *fpti, 
              float precip[VALID_CELLS], float tempmax[VALID_CELLS],
              float tempmin[VALID_CELLS], int cl);
const char flist[BUFSIZ+1]="latlong99.txt";
/* Temporary file that will be deleted. If the file exist it is overwritten */

void main(int argc, char **argv)     /* commandline arguments */
{
  FILE *fpprcp, *fptmax, *fptmin, *fpflist, *fpout;   /* filepointers used */
  char prcp[BUFSIZ+1];  char tmaxc[BUFSIZ+1];  char tminc[BUFSIZ+1]; /* strings holding filenames */
  char elev_mask[BUFSIZ+1];   /* string holding filename */
  char out_dir[BUFSIZ+1];     /* string holding out directory */
  char old_out_dir[BUFSIZ+1]; /* string holding the original output directory */
  char s1[BUFSIZ+1];              /* dummy string used to read strings */
  int cells, steps, tempsteps;
  int i,j,k,l, binflag;                    /* counters used in loops */
  float *prec,*tmax,*tmin,***BIGA;
  short int *precip,*tempmax,*tempmin,***BIG;
  float BIG_temp;

  if (argc<7) {           /* Must be exactly 6 arguments behind the program name */
  printf("Not correct number of commandline arguments \n");
  printf("vicinput \"prcp\" \"tmax\" \"tmin\" \"elev_mask.txt\" \"out_dir\" \n"); 
  exit(EXIT_FAILURE); }
  strcpy(prcp,argv[1]);        printf("%s \n",prcp);
  strcpy(tmaxc,argv[2]);        printf("%s \n",tmaxc);
  strcpy(tminc,argv[3]);        printf("%s \n",tminc);
  strcpy(elev_mask,argv[4]);   printf("%s \n",elev_mask);
  strcpy(out_dir,argv[5]); strcpy(old_out_dir,out_dir);printf("%s \n",old_out_dir);
  binflag=atoi(argv[6]);        printf("%d \n",binflag);

    /*--------------------------------------------------*/
    /* ALLOCATE MEMORY TO ARRAYS BASED ON BINARY STATUS */
    /*--------------------------------------------------*/
  if(binflag!=1) {
    if(!(prec = (float*) calloc(VALID_CELLS,sizeof(float)))) {
      printf("Memory Allocation error: PREC\n"); exit(8); }
    if(!(tmax = (float*) calloc(VALID_CELLS,sizeof(float)))) {
      printf("Memory Allocation error: TMAX\n"); exit(8); }
    if(!(tmin = (float*) calloc(VALID_CELLS,sizeof(float)))) {
      printf("Memory Allocation error: TMIN\n"); exit(8); }
    if(!(BIGA = (float***) calloc(3,sizeof(float**)))) {
      printf("Memory Allocation Error: BIGA\n"); exit(8); }
    for(j=0; j<3;j++) {
      if(!(BIGA[j] = (float**) calloc(stp50,sizeof(float*)))) {
	printf("Memory Allocation Error: BIGA\n");exit(8); }
      for(i=0; i<stp50;i++) {
	if(!(BIGA[j][i] = (float*) calloc(VALID_CELLS,sizeof(float)))) {
	  printf("Memory Allocation Error: BIGA\n"); exit(8); } } }
  }
  else {
    if(!(precip = (short int*) calloc(VALID_CELLS,sizeof(short int)))) {
      printf("Memory Allocation error: PREC\n"); exit(8); }
    if(!(tempmax = (short int*) calloc(VALID_CELLS,sizeof(short int)))) {
      printf("Memory Allocation error: TMAX\n"); exit(8); }
    if(!(tempmin = (short int*) calloc(VALID_CELLS,sizeof(short int)))) {
      printf("Memory Allocation error: TMIN\n"); exit(8); }
    if(!(BIG = (short int***) calloc(3,sizeof(short int**)))) {
      printf("Memory Allocation Error: BIGA\n"); exit(8); }
    for(j=0; j<3;j++) {
      if(!(BIG[j] = (short int**) calloc(stp50,sizeof(short int*)))) {
	printf("Memory Allocation Error: BIGA\n");exit(8); }
      for(i=0; i<stp50;i++) {
	if(!(BIG[j][i] = (short int*) calloc(VALID_CELLS,sizeof(short int)))) {
	  printf("Memory Allocation Error: BIGA\n"); exit(8); } } }
  }
  /*------------------------------------------------------------*/

  if (valid_cells(elev_mask)>VALID_CELLS) {   /* Check if array is big enough */
    printf("Number of valid cells = %d \n",valid_cells(elev_mask));
    printf("Change VALID_CELLS to %d in source\n",(valid_cells(elev_mask)+1));
    printf("Then recompile\n"); }
  else printf("Number of valid cells = %d \n",valid_cells(elev_mask));

  make_flist(elev_mask);    /* Generates a filelist that will be used later  */
  make_outfiles(old_out_dir,binflag); /* Open files in write mode, closes them again */
  cells = valid_cells(elev_mask); /* Calculates the number of valid cells in basin */
  printf("cells  = %d \n",cells);
  steps = nr_timesteps(prcp,cells,binflag);  
  printf("timesteps = %d \n",steps);

  if((fpflist = fopen(flist,"r"))==NULL) {             /* latlong99.txt  */
    printf("Cannot open file %s \n",flist);exit(0);} 
  if (binflag!=1) {
    if((fpprcp = fopen(prcp,"r"))==NULL){       /* prcp.grd after rescale */
      printf("Cannot open file %s \n",prcp);exit(0);}
    if((fptmax = fopen(tmaxc,"r"))==NULL){             /* tmax.grd      */
      printf("Cannot open file %s \n",tmaxc);exit(0);}
    if((fptmin = fopen(tminc,"r"))==NULL){              /* tmin.grd    */
      printf("Cannot open file %s \n",tminc);exit(0);} 
  }
  else { /* OPEN FILES FOR BINARY READ */
  if((fpprcp = fopen(prcp,"rb"))==NULL) {       /* prcp.grd after rescale */
    printf("Cannot open file %s \n",prcp);exit(0);}
  if((fptmax = fopen(tmaxc,"rb"))==NULL) {           /* tmax.grd      */
    printf("Cannot open file %s \n",tmaxc);exit(0);} 
  if((fptmin = fopen(tminc,"rb"))==NULL) {           /* tmin.grd    */
    printf("Cannot open file %s \n",tminc);exit(0);} 
  }

/***** RUNNING the remainder of steps/50     ****************************/  
  tempsteps = steps%50; /* running timesteps then write to files */
  printf("timestep %d \n",tempsteps);

  for (j=0;j<tempsteps;j++)
    {
      if(binflag!=1) {
	one_tstp_ascii(fpprcp,fptmax,fptmin,prec,tmax,tmin,cells);/* 1 timestep */ 
	for (i=0;i<cells;i++) {
	  BIGA[0][j][i]=prec[i];
	  BIGA[1][j][i]=tmax[i];
	  BIGA[2][j][i]=tmin[i]; }
      }
      else {
	one_tstp(fpprcp,fptmax,fptmin,precip, tempmax, tempmin,cells);/* 1 timestep */ 
      for (i=0;i<cells;i++) {           /* putting the numbers into BIG array */
	BIG[0][j][i]=precip[i];
	BIG[1][j][i]=tempmax[i];
	BIG[2][j][i]=tempmin[i]; }
      }
    }

    for(k=0;k<cells;k++)    /* How many files to run */
      {
      fscanf(fpflist,"%s",s1);    
      strcpy(out_dir, old_out_dir);
      strcat(out_dir, s1);

      if(binflag!=1) {if((fpout = fopen(out_dir,"a"))==NULL){
	printf("Cannot open file %s \n",out_dir);exit(0);}}
      else  {if((fpout = fopen(out_dir,"ab"))==NULL){
	printf("Cannot open file %s \n",out_dir);exit(0);} }

      for (i=0;i<tempsteps;i++) {  /* BIG [unit][timestep][gridcell] */
	if(binflag!=1) BIG_temp = BIGA[0][i][k];
	else BIG_temp = (float) BIG[0][i][k];
	if(BIG_temp<0.00) {
	  fprintf(stderr,"neg. prcp reset to 0 \n");
	  if(binflag!=1) BIGA[0][i][k] = 0;
	  else BIG[0][i][k] = 0; }
	if(binflag!=1)
	  fprintf(fpout,"%4.2f  %4.2f  %4.2f\n",
		  BIGA[0][i][k],BIGA[1][i][k],BIGA[2][i][k]);
	else { /* WRITE OUT IN BINARY FORMAT */
	  fwrite(&BIG[0][i][k],sizeof(short int),1,fpout);
	  fwrite(&BIG[1][i][k],sizeof(short int),1,fpout);
	  fwrite(&BIG[2][i][k],sizeof(short int),1,fpout);
	}
      }
      fclose(fpout);
      }
  rewind(fpflist);
/****** END running the remainder of steps/50  ********************/

/***** RUNNING steps/50     ****************************/  
  if (steps>50) 
    {
    tempsteps = steps/stp50;
    for (l=0;l<tempsteps;l++) {
      printf(" timestep %d \n",(steps%50 + (1+l)*50));
      for (j=0;j<stp50;j++)
        {
      if(binflag!=1) {
	one_tstp_ascii(fpprcp,fptmax,fptmin,prec, tmax, tmin,cells);
	for (i=0;i<cells;i++) {    /* putting the numbers into BIG array */
	  BIGA[0][j][i]=prec[i];
	  BIGA[1][j][i]=tmax[i];
	  BIGA[2][j][i]=tmin[i]; }
      } 
      else {
	one_tstp(fpprcp,fptmax,fptmin,precip, tempmax, tempmin,cells);
      for (i=0;i<cells;i++) {    /* putting the numbers into BIG array */
	BIG[0][j][i]=precip[i];
	BIG[1][j][i]=tempmax[i];
	BIG[2][j][i]=tempmin[i]; }
      }
    }

      for(k=0;k<cells;k++)    /* How many files to run */
        {
      fscanf(fpflist,"%s",s1);    
      strcpy(out_dir, old_out_dir);
      strcat(out_dir, s1);

      if(binflag!=1) fpout = fopen(out_dir,"a");
      else fpout = fopen(out_dir,"ab");

      for (i=0;i<stp50;i++) {  /* BIG [unit][timestep][gridcell] */
	/*reconvert precip to actual values for check */
	if(binflag!=1) BIG_temp = BIGA[0][i][k];
	else BIG_temp = (float) (BIG[0][i][k]);
	if(BIG_temp<0.00) {
	  fprintf(stderr,"neg. prcp reset to 0 \n");
	  if(binflag!=1) BIGA[0][i][k] = 0;
	  else BIG[0][i][k] = 0; }
	if(binflag!=1)
	  fprintf(fpout,"%4.2f  %4.2f  %4.2f\n",
		  BIGA[0][i][k],BIGA[1][i][k],BIGA[2][i][k]);
	else { /* WRITE OUT IN BINARY FORMAT */
	  fwrite(&BIG[0][i][k],sizeof(short int),1,fpout);
	  fwrite(&BIG[1][i][k],sizeof(short int),1,fpout);
	  fwrite(&BIG[2][i][k],sizeof(short int),1,fpout);
	}
      }

      fclose(fpout);
      }
     rewind(fpflist);
     }
    }
/****** END running steps/50 **********************************************/

  fclose(fpprcp);  fclose(fptmax);  fclose(fptmin);  fclose(fpflist);
}
/*** END main    *************/ 

/**************************************************************************/
void make_flist(char msk[BUFSIZ+1])
{
/* This function reads a mask file and generates a text file         */
/* that contains the latitude and longitude to each gridcell         */
/* This will be the names of all the vicinput files                  */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

  FILE *fptrmask, *fptrlatlong;                /* Filepointers for to files*/
  extern const char flist[BUFSIZ+1];
  char maskfile[BUFSIZ+1];
  char latlong[BUFSIZ+1];
  char str1[BUFSIZ+1];
  char void_nr[BUFSIZ+1];
  int cols, rows,i,j; 
  float uplftlat, uplftlong, llftlat,llftlong,resolution,longitude,value;
  float high, low;

  strcpy(latlong,flist);
  printf("latlong %s \n",latlong);
  strcpy(maskfile,msk);
  printf("mask %s \n",maskfile);

  if((fptrmask = fopen(maskfile,"r"))==NULL){        /* Opens maskfile */
    printf("Cannot open file maskfile %s \n",maskfile);  exit(0);}
  if((fptrlatlong = fopen(latlong,"w"))==NULL){      /* Opens latlong file */
    printf("Cannot open file %s \n",latlong);  exit(0);}

  fscanf(fptrmask,"%*s %s",str1);        /* Reads number of cols */
  printf("%s\n",str1);  cols = atoi(str1);  /* Assigns number of cols to cols*/
  fscanf(fptrmask,"%*s %s",str1);        /* Reads number of rows */
  printf("%s\n",str1); rows = atoi(str1);  /* Assigns number of rows to rows*/
  fscanf(fptrmask,"%*s %s", str1);        /* Reads llftlong gridposition */
  printf("%s\n",str1);  llftlong = atof(str1);  /* Assigns gridposition to llftlong */
  fscanf(fptrmask,"%*s %s", str1);        /* Reads llftlat gridposition */
  printf("%s\n",str1); llftlat = atof(str1);  /* Assigns gridposition to llftlat*/  
  fscanf(fptrmask,"%*s %s",str1);        /* Reads resolution in degrees */
  printf("%s\n",str1);resolution= atof(str1);     /* Assigns the resolution */
  fscanf(fptrmask,"%*s %s",void_nr);     /* Reads the void number */
  printf("%s \n",void_nr);

  /* place nodes at center of grid cells */
  uplftlong = llftlong + (resolution/2);
  uplftlat  = llftlat + ((float)(rows-1)*resolution) + (resolution/2);
  printf("uplftlong %.4f \n",uplftlong);
  printf("uplftlat %.4f \n",uplftlat);

  for(i=0;i<rows;i++)               /* Write lat long to latlong.txt */
    {
      longitude = uplftlong;  
      for(j=0;j<cols;j++)
        {
         fscanf(fptrmask,"%s", str1);
         value = atof(str1);
	 high=atof(void_nr)+0.001;
	 low=atof(void_nr)-0.001;
	 if (value>high || value<low ){
          fprintf(fptrlatlong,"data_%.4f_",uplftlat);
          fprintf(fptrlatlong,"%.4f\n",longitude);}
        longitude = longitude + resolution;     /* keep track of position */
	}
    uplftlat = uplftlat - resolution;           /* keep track of position */
    }
  fclose(fptrmask);
  fclose(fptrlatlong);
} /* END of function void make_files(FILE *fp, char dir[400]) *********/
/***********************************************************************/
void make_outfiles(char out_dir[BUFSIZ+1], int binflag)
{
  FILE *fpflist, *fpout;           /* filepointers to the files */
  char out[BUFSIZ+1];                   /* name of outputfile */
  extern const char flist[BUFSIZ+1];    /* name of filelist   "latlong99.txt"   */   
  char str1[BUFSIZ+1];                   /* dummy string  */

  if((fpflist = fopen(flist,"r"))==NULL){       /* opens filelist   */
    printf("Cannot open file %s \n",flist);  exit(0);}

  while(fscanf(fpflist,"%s", str1)!=EOF)    /* read the filelist to EOF */
    {
      /*      fscanf(fpflist,"%s", str1); */
      strcpy(out,out_dir);        /* copy directory name */
      strcat(out,str1);           /* add filename   */
      if(binflag!=1){
	fpout = fopen(out,"w");
      }
      else {   /* open file in write mode -- BINARY*/
	fpout = fopen(out,"wb");
      }
      fclose(fpout);  /* close the file. It will be reopend in append mode */
    }
  fclose(fpflist);
}/** END function void make_outfiles(char out_dir[400] *********/
/************************************************************************/
int valid_cells(char maskf[BUFSIZ+1])  /* function returns number of valid cells in mask */
{
  FILE *fpmaskf;
  char str1[BUFSIZ+1],voidstr[BUFSIZ+1];
  int i;
  float high, low, valmask;
  if((fpmaskf = fopen(maskf,"r"))==NULL){
     printf("Cannot open file %s \n",maskf);  exit(0);}
  for (i=0;i<12;i++)
    fscanf(fpmaskf,"%s",voidstr); /* Skipping the first lines in mask*/
  i = 0;
  while(fscanf(fpmaskf,"%s", str1)!=EOF) {
    valmask=atof(voidstr);
    high=valmask+0.001;
    low=valmask-0.001;
    if (atof(str1)>high || atof(str1)<low) i++; /* if not no_data */
  }
  return i;
  fclose(fpmaskf);
}/*** END of function void valid_cells(char maskf[400]) ********************/
/***************************************************************************/
/* This function reads one timestep into three arrays - BINARY */
void one_tstp(FILE *fppp, FILE *fpta, FILE *fpti, 
              short int *precip, short int *tempmax,
	      short int *tempmin, int cl)
{
 int i;
 short int ttemp;

 for (i=0;i<cl;i++)
   {
   if (fread(&precip[i],2,1,fppp)  !=1)  printf("Error prcp\n");
   if (fread(&tempmax[i],2,1,fpta) !=1)  printf("Error tmax\n");
   if (fread(&tempmin[i],2,1,fpti) !=1)  printf("Error tmin\n");
   if (tempmax[i]<tempmin[i]) {
     ttemp = tempmax[i];
     tempmax[i] = tempmin[i];
     tempmin[i] = ttemp;}
   }
} /*** END void one_tstp ***/
/**************************************************************************/
void one_tstp_ascii(FILE *fppp, FILE *fpta, FILE *fpti, 
		    float *precip, float *tempmax,
		    float *tempmin, int cl)
{
 int i;
 float ttemp;

 for (i=0;i<cl;i++)
   {
   if (fscanf(fppp,"%f",&precip[i])!=1)  printf("Error prcp\n"); 
   if (fscanf(fpta,"%f",&tempmax[i])!=1) printf("Error tmax\n");
   if (fscanf(fpti,"%f",&tempmin[i])!=1) printf("Error tmin\n");
   if (tempmax[i]<tempmin[i]) {
     ttemp = tempmax[i];
     tempmax[i] = tempmin[i];
     tempmin[i] = ttemp;}
   }
}/*** END void one_tstp_ascii ***/
/**************************************************************************/
int nr_timesteps(char pp[BUFSIZ+1], int cc, int binflag)
{
  /* This function returns the number of timesteps in the three inputfiles */
  FILE *fp1;
  char str1[BUFSIZ+1];
  short int dummy;
  long i;
  int r;
  if(binflag!=1) {
    if((fp1 = fopen(pp,"r"))==NULL){      /* opens file  */
    printf("Cannot open file %s \n",pp);  exit(0);} 
  }
  else {
    if((fp1 = fopen(pp,"rb"))==NULL){   /* NEW FOR BINARY READ */
    printf("Cannot open file %s \n",pp);  exit(0);} 
  }
  i =0;
  if(binflag!=1) {
  while (fscanf(fp1,"%s",str1)!=EOF)i++;
  }
  else{
    while(!feof(fp1)) {
      fread(&dummy,2,1,fp1);
      i++;
    }
  }
  r =(int)(i/cc);
  fclose(fp1);
  return r;
}/*** END function int nr_timesteps(char pp[400]) ********/
/*********************************************************/
