/* File:             mk_monthly.c                                             */
/* Programmer :      Bernt Viggo Matheussen                                   */
/* Modified:         4/2000 by EdM to use binary flag and multiplier          */
/* Date:             02.06.98                                                 */
/* Version:          1.0                                                      */

/* This program reads the output file from the regridding program (test.grd). */
/* It then generates monthly values for jan,feb,mar,etc that will be used     */
/* together with PRISM data to scale the timeseries in each gridcell          */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#define MAX_NR_OF_CELLS (10000)  /*Columbia River  0.25dg 64*56 = 3584 cells */ 
#define MAXLINE 80
typedef struct str_big {float jan[MAX_NR_OF_CELLS];float feb[MAX_NR_OF_CELLS];
                        float mar[MAX_NR_OF_CELLS];float apr[MAX_NR_OF_CELLS];
                        float may[MAX_NR_OF_CELLS];float jun[MAX_NR_OF_CELLS];
                        float jul[MAX_NR_OF_CELLS];float aug[MAX_NR_OF_CELLS];
                        float sep[MAX_NR_OF_CELLS];float okt[MAX_NR_OF_CELLS];
                        float nov[MAX_NR_OF_CELLS];float dec[MAX_NR_OF_CELLS];
}STR_BIG;

int valid_cells(FILE *pp);  /* function returns the number of cells which are not void in maskfile */ 
int number_of_timesteps(FILE *fp1, int valcell, int binflag); /* function returns number of timesteps */
void set_0(STR_BIG *pk);              /* sets all values to 0 in STR_BIG *pk     */
void read_one_timestep(FILE *fptgrd,float timestep[], int nr_steps, float mult,
		       int binflag);   /* Reads one timestep in data */
void write_data(STR_BIG *write, FILE *mask,int ys,char out_dr[]);
int tstp_to_start(int start_yy_data, int start_yy_prd);  

void main(int argc, char **argv)
{ 
  FILE *fpmask, *fpgrd;                /* Filepointers for to files*/
  char maskfile[150];
  char grdfile[150];  /*=".grd";*/
  char outp_dir[250]; /* ="out/";*/
  char str1[12];
  char str2[12];
  int nr_of_valid_cells,nr_tstep, i,start_year,start_month;
  int present_tstep, x;
  int run_years, yy, binflag;
  int to_start, start_period, end_period;
  float prec, mult;
  float tstp[MAX_NR_OF_CELLS];   /* Not using pointer could be dangerous */

  STR_BIG *pt;
  pt = malloc(sizeof(STR_BIG));    /* allocates memory to the structure */

  if (argc<9) {        /* Must be exactly 6 arguments behind the program name */
  printf("Not correct number of commandline arguments \n");
  exit(EXIT_FAILURE); }
  strcpy(maskfile,argv[1]);
  strcpy(grdfile,argv[2]);
  strcpy(outp_dir,argv[3]);
  strcpy(str1,argv[4]); start_year=atoi(str1);
  strcpy(str1,argv[5]); start_period=atoi(str1);
  strcpy(str1, argv[6]); end_period=atoi(str1);
  strcpy(str1, argv[7]); binflag=atoi(str1); /*binary flag */
  strcpy(str1, argv[8]); mult=atof(str1); /*data multiplier*/
  printf("Commandline arguments\n");
  printf("mk_monthly ");
  printf("%s %s %s %d %d %d\n",maskfile,grdfile,outp_dir,start_year, start_period,end_period);
          
  if((fpmask = fopen(maskfile,"r"))==NULL){
    printf("Cannot open file %s \n",maskfile);exit(0);}

  if(binflag!=1) {
    if((fpgrd = fopen(grdfile,"r"))==NULL){
      printf("Cannot open file %s \n",grdfile); exit(0);}
  }
  else {
    if((fpgrd = fopen(grdfile,"rb"))==NULL){   
      printf("Cannot open file %s \n",grdfile); exit(0);}
  }

  printf("Maximum number of cells in this program is %d\n",MAX_NR_OF_CELLS);

  nr_of_valid_cells = valid_cells(fpmask);
  printf("Number of valid cells are %d\n",nr_of_valid_cells);  
  nr_tstep = number_of_timesteps(fpgrd, nr_of_valid_cells, binflag);
  printf("Number of timesteps are %d\n",nr_tstep);

  set_0(pt);
  
  printf("start_year = %d\n",start_year);
  printf("start_period = %d\n",start_period);
  printf("end_period = %d\n",end_period);
  to_start = tstp_to_start(start_year,start_period);
  printf("Timesteps to start = %d\n",to_start);

  for (i=0;i<to_start;i++)read_one_timestep(fpgrd,tstp,nr_of_valid_cells,mult,binflag); /* Moving down */

  run_years = tstp_to_start(start_period,(1+end_period));
  printf("Timesteps to run = %d\n",run_years);

  run_years = end_period - start_period +1;

  printf("Number of years to run = %d \n",run_years);
  prec =0.0;
  for(yy=0;yy<run_years;yy++)
    {
    for (i = 0;i<31;i++){  /* January */
      read_one_timestep(fpgrd,tstp,nr_of_valid_cells,mult,binflag);
      for (x=0;x<nr_of_valid_cells;x++) pt->jan[x] = pt->jan[x]+tstp[x]; }

    if (start_period%4!=0) {   /*February */
      for (i =0;i<28;i++) {  
        read_one_timestep(fpgrd,tstp,nr_of_valid_cells,mult,binflag);
        for (x=0;x<nr_of_valid_cells;x++) pt->feb[x] = pt->feb[x]+tstp[x]; }
        printf("Regular year  %d \n",start_period); }

    if (start_period%4==0) {  /* February*/
      for (i =0;i<29;i++) {  
        read_one_timestep(fpgrd,tstp,nr_of_valid_cells,mult,binflag);
        for (x=0;x<nr_of_valid_cells;x++) pt->feb[x] = pt->feb[x]+tstp[x]; }
        printf("Leap year  %d \n",start_period); }

    for (i = 0;i<31;i++){  /* March */
      read_one_timestep(fpgrd,tstp,nr_of_valid_cells,mult,binflag);
      for (x=0;x<nr_of_valid_cells;x++) pt->mar[x] = pt->mar[x]+tstp[x]; }

    for (i = 0;i<30;i++){  /* April */
      read_one_timestep(fpgrd,tstp,nr_of_valid_cells,mult,binflag);
      for (x=0;x<nr_of_valid_cells;x++) pt->apr[x] = pt->apr[x]+tstp[x]; }

   for (i = 0;i<31;i++){  /* May */
      read_one_timestep(fpgrd,tstp,nr_of_valid_cells,mult,binflag);
      for (x=0;x<nr_of_valid_cells;x++) pt->may[x] = pt->may[x]+tstp[x]; }

   for (i = 0;i<30;i++){  /* June */
      read_one_timestep(fpgrd,tstp,nr_of_valid_cells,mult,binflag);
      for (x=0;x<nr_of_valid_cells;x++) pt->jun[x] = pt->jun[x]+tstp[x]; }

   for (i = 0;i<31;i++){  /* July */
      read_one_timestep(fpgrd,tstp,nr_of_valid_cells,mult,binflag);
      for (x=0;x<nr_of_valid_cells;x++) pt->jul[x] = pt->jul[x]+tstp[x]; }

   for (i = 0;i<31;i++){  /* August */
      read_one_timestep(fpgrd,tstp,nr_of_valid_cells,mult,binflag);
      for (x=0;x<nr_of_valid_cells;x++) pt->aug[x] = pt->aug[x]+tstp[x]; }

   for (i = 0;i<30;i++){  /* September */
      read_one_timestep(fpgrd,tstp,nr_of_valid_cells,mult,binflag);
      for (x=0;x<nr_of_valid_cells;x++) pt->sep[x] = pt->sep[x]+tstp[x]; }

   for (i = 0;i<31;i++){  /* Oktober */
      read_one_timestep(fpgrd,tstp,nr_of_valid_cells,mult,binflag);
      for (x=0;x<nr_of_valid_cells;x++) pt->okt[x] = pt->okt[x]+tstp[x]; }

   for (i = 0;i<30;i++){  /* November */
      read_one_timestep(fpgrd,tstp,nr_of_valid_cells,mult,binflag);
      for (x=0;x<nr_of_valid_cells;x++) pt->nov[x] = pt->nov[x]+tstp[x]; }

   for (i = 0;i<31;i++){  /* December */
      read_one_timestep(fpgrd,tstp,nr_of_valid_cells,mult,binflag);
      for (x=0;x<nr_of_valid_cells;x++) pt->dec[x] = pt->dec[x]+tstp[x]; }

   start_period = start_period +1;
   } /* END of yearly loop  */

  write_data(pt,fpmask,run_years,outp_dir);  /* Function writes the output file */

  free(pt);   /***** Removes memory ****/
  fclose(fpgrd);
  fclose(fpmask);
/*** END main *********************/
}

/*******************************************************************************/
int valid_cells(FILE *pp)  /* function returns the number of cells which are not void in maskfile */
{
  int cnt; 
  char s1[MAXLINE]; 
  char cdum[MAXLINE];
  float void_nr,up,down; 
  for (cnt =0;cnt<6;cnt++) fgets(s1,MAXLINE,pp); /* Moving down to data*/
  sscanf(s1,"%s %f", &cdum, &void_nr);
  up = void_nr + 0.0001;
  down = void_nr - 0.0001;
  printf("Void is now = %.0f\n",void_nr);
  cnt = 0;                         /*  Avoid comparing float using '==' */
  while (fscanf(pp,"%s",&s1)!=EOF){if(atof(s1)>up||atof(s1)<down)cnt++;}
  rewind(pp);
  return cnt;
/*** END function valid_cells ***/
}

/*******************************************************************************/
int number_of_timesteps(FILE *fp1, int valcell, int binflag){
  char s2[10];
  short int dummy_int;
  int j;
  j=0;
  if(binflag!=1){
    while (fscanf(fp1,"%s",&s2)!=EOF)j++; /*used for ascii data */
  }
  else {
  while((fread(&dummy_int,sizeof(short int),1,fp1) != 0)) j++; /*Binary data*/
  }
  rewind(fp1);
  return (j/valcell);}

/*******************************************************************************************/
void set_0(STR_BIG *pk){
 int k;
 for (k=0;k<MAX_NR_OF_CELLS;k++) pk->jan[MAX_NR_OF_CELLS] = 0.0;}

/**********************************************************************************************/
void read_one_timestep(FILE *fptgrd,float timestep[],int nr_steps,float mult,int 
		       binflag)
{
  int i;
  short int stmp;
  char ss[10];
  for (i = 0;i<nr_steps;i++) {
    if(binflag!=1) {
      fscanf(fptgrd,"%s",&ss); /* ascii data */
      timestep[i] = atof(ss);
    }
    else {
      fread(&stmp,sizeof(short),1,fptgrd); /*Binary data */
      timestep[i]=(float)stmp/mult; 
    }
  }
}


/*******************************************************************************************/

void write_data(STR_BIG *write, FILE *mask,int ys,char out_dr[])
{
  FILE *fpjan, *fpfeb, *fpmar, *fpapr, *fpmay, *fpjun;
  FILE *fpjul, *fpaug, *fpsep, *fpokt, *fpnov, *fpdec;

  char s1[MAXLINE];  
  char cdum[MAXLINE];
  float void_nr,up,down;  /* Directory name and filename must not be more than 200 char */
  int i;

  char jan[200]="monthly.jan";  char feb[200]="monthly.feb";  char mar[200]="monthly.mar";
  char apr[200]="monthly.apr";  char may[200]="monthly.may";  char jun[200]="monthly.jun";
  char jul[200]="monthly.jul";  char aug[200]="monthly.aug";  char sep[200]="monthly.sep";
  char okt[200]="monthly.oct";  char nov[200]="monthly.nov";  char dec[200]="monthly.dec";

  char void_str[10];
  char old_dir[200];
  int k, cols, rows, cells;

  printf("Writing monthly values \n");

  strcpy(old_dir,out_dr);
  
  strcat(out_dr,jan);  strcpy(jan,out_dr); strcpy(out_dr,old_dir);
  strcat(out_dr,feb);  strcpy(feb,out_dr); strcpy(out_dr,old_dir);
  strcat(out_dr,mar);  strcpy(mar,out_dr); strcpy(out_dr,old_dir);
  strcat(out_dr,apr);  strcpy(apr,out_dr); strcpy(out_dr,old_dir);
  strcat(out_dr,may);  strcpy(may,out_dr); strcpy(out_dr,old_dir);
  strcat(out_dr,jun);  strcpy(jun,out_dr); strcpy(out_dr,old_dir);
  strcat(out_dr,jul);  strcpy(jul,out_dr); strcpy(out_dr,old_dir);
  strcat(out_dr,aug);  strcpy(aug,out_dr); strcpy(out_dr,old_dir);
  strcat(out_dr,sep);  strcpy(sep,out_dr); strcpy(out_dr,old_dir);
  strcat(out_dr,okt);  strcpy(okt,out_dr); strcpy(out_dr,old_dir);
  strcat(out_dr,nov);  strcpy(nov,out_dr); strcpy(out_dr,old_dir);
  strcat(out_dr,dec);  strcpy(dec,out_dr); strcpy(out_dr,old_dir);

  fpjan = fopen(jan,"w");   fpfeb = fopen(feb,"w");   fpmar = fopen(mar,"w"); 
  fpapr = fopen(apr,"w");   fpmay = fopen(may,"w");   fpjun = fopen(jun,"w"); 
  fpjul = fopen(jul,"w");   fpaug = fopen(aug,"w");   fpsep = fopen(sep,"w"); 
  fpokt = fopen(okt,"w");   fpnov = fopen(nov,"w");   fpdec = fopen(dec,"w"); 
  printf("output files opened successfully\n");

  rewind(mask);

  for(i=0;i<6;i++){
  fgets(s1,MAXLINE,mask);
  fprintf(fpjan,"%s",s1);  fprintf(fpfeb,"%s",s1);
  fprintf(fpmar,"%s",s1);  fprintf(fpapr,"%s",s1);
  fprintf(fpmay,"%s",s1);  fprintf(fpjun,"%s",s1);
  fprintf(fpjul,"%s",s1);  fprintf(fpaug,"%s",s1);
  fprintf(fpsep,"%s",s1);  fprintf(fpokt,"%s",s1); 
  fprintf(fpnov,"%s",s1);  fprintf(fpdec,"%s",s1);
  }
  rewind(mask);

  fgets(s1,MAXLINE,mask);
  sscanf(s1,"%s %d", &cdum, &cols);
  fgets(s1,MAXLINE,mask);
  sscanf(s1,"%s %d", &cdum, &rows);
  fgets(s1,MAXLINE,mask);
  fgets(s1,MAXLINE,mask);
  fgets(s1,MAXLINE,mask);
  fgets(s1,MAXLINE,mask);
  sscanf(s1,"%s %f", &cdum, &void_nr);
  sprintf(void_str, "%.2f", void_nr);
  printf("cols, %d\trows, %d\tvoid, %.0f\n", cols, rows, void_nr);

  up = void_nr + 0.0001;   
  down = void_nr - 0.0001;   
  k = 0;
  cells = 0;
  while (fscanf(mask,"%s",&s1)!=EOF)  
    {
    cells++;
    if(atof(s1)>up||atof(s1)<down){  /*  Avoid comparing float using '==' */
      fprintf(fpjan,"%4.2f ",(write->jan[k]/ys)); fprintf(fpfeb,"%4.2f ",(write->feb[k]/ys));
      fprintf(fpmar,"%4.2f ",(write->mar[k]/ys)); fprintf(fpapr,"%4.2f ",(write->apr[k]/ys));
      fprintf(fpmay,"%4.2f ",(write->may[k]/ys)); fprintf(fpjun,"%4.2f ",(write->jun[k]/ys));
      fprintf(fpjul,"%4.2f ",(write->jul[k]/ys)); fprintf(fpaug,"%4.2f ",(write->aug[k]/ys));
      fprintf(fpsep,"%4.2f ",(write->sep[k]/ys)); fprintf(fpokt,"%4.2f ",(write->okt[k]/ys));
      fprintf(fpnov,"%4.2f ",(write->nov[k]/ys)); fprintf(fpdec,"%4.2f ",(write->dec[k]/ys));
      k++; }
    else  {
      fprintf(fpjan,"%s ",void_str);  fprintf(fpfeb,"%s ",void_str);
      fprintf(fpmar,"%s ",void_str);  fprintf(fpapr,"%s ",void_str);
      fprintf(fpmay,"%s ",void_str);  fprintf(fpjun,"%s ",void_str);
      fprintf(fpjul,"%s ",void_str);  fprintf(fpaug,"%s ",void_str);
      fprintf(fpsep,"%s ",void_str);  fprintf(fpokt,"%s ",void_str);
      fprintf(fpnov,"%s ",void_str);  fprintf(fpdec,"%s ",void_str); }
    if (cells==cols) {
      fprintf(fpjan,"\n");cells=0;  fprintf(fpfeb,"\n");cells=0;     
      fprintf(fpmar,"\n");cells=0;  fprintf(fpapr,"\n");cells=0;     
      fprintf(fpmay,"\n");cells=0;  fprintf(fpjun,"\n");cells=0;     
      fprintf(fpjul,"\n");cells=0;  fprintf(fpaug,"\n");cells=0;     
      fprintf(fpsep,"\n");cells=0;  fprintf(fpokt,"\n");cells=0;     
      fprintf(fpnov,"\n");cells=0;  fprintf(fpdec,"\n");cells=0; } 
    }

  fclose(fpjan);  fclose(fpfeb); fclose(fpmar); 
  fclose(fpapr);  fclose(fpmay); fclose(fpjun); 
  fclose(fpjul);  fclose(fpaug); fclose(fpsep); 
  fclose(fpokt);  fclose(fpnov); fclose(fpdec); 
}/*END function write_data*/
/************************************************************************************************/

int tstp_to_start(int start_yy_data, int start_yy_prd)
{
  int value,i,diff;
  value =0;
  if (start_yy_data > start_yy_prd) printf("ERROR not correct commandline arguments \n");
  diff = start_yy_prd - start_yy_data;
  if (start_yy_data < start_yy_prd) {
    for (i=0;i<diff;i++) {
      if (start_yy_data%4==0) value = value +366;
      if (start_yy_data%4!=0) value = value +365;
      start_yy_data++; }     }
  return value;
}/** End of function int tstp_to_start(int start_yy_data, int start_yy_prd)*****/
/*******************************************************************************/

