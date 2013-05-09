/*
 * SUMMARY:      This program reads monthly winddata from Surface Airways
 *               CDs. The rawdata have to be formatted with format_wind.scr 
 *               before this program can read the windfields.
 *               The VIC inputdata (Prcp, Tmax,Tmin) is read and 
 *               a daily timeseries of wind is added as a fourth 
 *               coloumn to the metdata.
 * USAGE:        Used together with the VIC model.
 *
 * AUTHOR:       Bernt Viggo Matheussen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       tempbvm@hydro.washington.edu
 * ORIG-DATE:    24-Aug-98 at 15:34:08
 * LAST-MOD:     Tue Feb 23 17:29:02 1999 by Bernt Viggo Matheussen
 * DESCRIPTION:  
 * DESCRIP-END.
 * FUNCTIONS:    read_wind(), read_metdata(), 
 * COMMENTS:     
 */
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <float.h>
#include <math.h>
#include <string.h>

#define CELLS 4858
#define STEPS 4749
#define WINDSTATIONS 65
#define STARTYEAR 1980     /* Start year of input data. Must start 1 january */ 
#define STOPYEAR 1992      /* Stop year must end 31 december */



/* Reads in the wind reanalasys data from file */
void read_wind(char windfile[200],float WIND[WINDSTATIONS][15]);

/* Reads met data Prcp, Tmax, Tmin */
void read_metdata(char file[200], float METDATA[STEPS][4]);

/* Returns days in month */
int days_in_month(int YY, int MM);

/* Calculates the distance between to pairs of latitude and longitude */
/* This function is the same as latlong.f function in the regrid programs */
float calc_distance(float latcell,float loncell,float latstat,float lonstat);

/* Extracts latitude from filename */
float find_latitude(char filename[200]);

/* Extracts longitude from filename */
float find_longitude(char filename[200]);

/* Finds the two nearest neighbours and puts data into array */
void sort_stations(float latcell, float loncell, 
                   float WIND[WINDSTATIONS][15],   float CELLWIND[12]);

/* Writes daily windseries into array */
void make_windseries( float METDATA[STEPS][4], float CELLWIND[12]);

/* Write metdata to new location */
void write_data(char file[200], float METDATA[STEPS][4]);


void main(void)
{
  FILE *fp;
  int i;
  char windfile[200]="WIND_CBR8.TXT";
  char filelist[200]="filelist_data80-92.txt";
  char metfile[200];
  char indir[200]="/home/tempbvm/cbr8/metdata/data80-92/";
  char outdir[200]="/home/tempbvm/cbr8/metdata/wind80-92/";
  float WIND[WINDSTATIONS][15]; /* Array contains gridcells and wind data + distance */
  float METDATA[STEPS][4]; /* Array stores prcp, tmax, tmin, wind  */
  char readfile[200];
  char writefile[200];
  float latcell, loncell;

  float CELLWIND[12]; /* Monthly wind data for the computational gridcell */

  read_wind(windfile,WIND); /* Reading wind into array */


   if((fp = fopen(filelist,"r"))==NULL)  
      {printf("Cannot open file %s \n",filelist);exit(0);}  

    for (i=0;i<CELLS;i++)   
      {  
      fscanf(fp,"%s",metfile);   /* Reads filename */ 
      strcpy(readfile,indir);    /* Make output filename */ 
      strcat(readfile, metfile);    
      printf("%s\n",readfile); 
      strcpy(writefile,outdir);    
      strcat(writefile, metfile);    
      printf("%s\n",writefile); 
      read_metdata(readfile,METDATA); /* Reads metdata */
      latcell = find_latitude(metfile);  printf("latcell %6.4f\n",latcell); 
      loncell = find_longitude(metfile); printf("loncell %6.4f\n",loncell); 
      /* Sorting stations and calculating the weighted mean */
      sort_stations(latcell,loncell,WIND, CELLWIND); 

      /* Adding the fourth coloumn, wind to timeseries */
      make_windseries(METDATA,CELLWIND);      
      write_data(writefile, METDATA); 
      }  



    fclose(fp); 
}/* END MAIN ******************************************************************/

void read_wind(char windfile[200],float WIND[WINDSTATIONS][15])
{
  FILE *fp;
  int i,j;
  if((fp = fopen(windfile,"r"))==NULL) 
     {printf("Cannot open file %s \n",windfile);exit(0);} 
  for (i=0;i<WINDSTATIONS;i++) for (j=0;j<14;j++) fscanf(fp,"%f",&WIND[i][j]);
  fclose(fp);
}/* END function void read_wind() */

void read_metdata(char file[200], float METDATA[STEPS][4])
{
  FILE *fp;
  int i;
  if((fp = fopen(file,"r"))==NULL) 
     {printf("Cannot open file %s \n",file);exit(0);} 
  for (i=0;i<STEPS;i++) 
    {
    fscanf(fp,"%f %f %f",&METDATA[i][0], &METDATA[i][1], &METDATA[i][2]);
    }
  fclose(fp);
}/* END function  void read_metdata()  */		  

int days_in_month(int YY, int MM)
{
  int days;
  if (MM ==  1) days = 31;
  if (MM ==  2 && YY%4 == 0) days = 29; 
  if (MM ==  2 && YY%4 != 0) days = 28; 
  if (MM ==  3) days = 31;
  if (MM ==  4) days = 30;
  if (MM ==  5) days = 31;
  if (MM ==  6) days = 30;
  if (MM ==  7) days = 31;
  if (MM ==  8) days = 31;
  if (MM ==  9) days = 30;
  if (MM == 10) days = 31;
  if (MM == 11) days = 30;
  if (MM == 12) days = 31;
  return days;
}

float calc_distance(float latcell,float loncell,float latstat,float lonstat)
{
  /*  function latlong(lat1,long1,lat2,long2)            */
  /*  find distance between pair of lat long coordinates */
  /*  south and west are negative                        */
  float lat1, long1, lat2, long2;
  float pi=3.141593;
  float radius=6378.0;
  float dtor=2.0*pi/360.0; /* Degrees to radians */
  float theta1, theta2, phi1, phi2, temp;
  float term1, term2, term3;

  lat1=latcell;
  long1=loncell;
  lat2=latstat;
  long2=lonstat;

  theta1 = dtor*long1;
  phi1 = dtor*lat1;
  theta2 = dtor*long2;
  phi2 = dtor*lat2;
  term1 = cos(phi1)*cos(theta1)*cos(phi2)*cos(theta2);
  term2 = cos(phi1)*sin(theta1)*cos(phi2)*sin(theta2);
  term3 = sin(phi1)*sin(phi2);
  temp = term1+term2+term3;
  if (temp<1.0) temp=temp;
  if(temp>1.0)  temp=1.0;
  return (radius*acos(temp));
}/* END function float calc_distance() */


float find_latitude(char filename[200])
{
  int i;
  char name[200];
  char junk[25];
  float return_lat;
  strcpy(name,filename);    /* Replacing all _ with space */
  for(i=0;i<((strlen(filename))-1);i++) {if(name[i] == '_') name[i] = ' '; } 
  sscanf(name,"%s %f %s",junk, &return_lat,junk); /* Scanning the string */
  return return_lat;
}/* float find_latitude() */

float find_longitude(char filename[200])
{
  int i;
  char name[200];
  char junk[25];
  float return_lon;
  strcpy(name,filename);    /* Replacing all _ with space */
  for(i=0;i<((strlen(filename))-1);i++) {if(name[i] == '_') name[i] = ' '; } 
  sscanf(name,"%s %s %f",junk,junk,&return_lon); /* Scanning the string */
  return return_lon;
}/* END function float find_longitude() */

/* Finds the two nearest neighbours and puts data into array */
void sort_stations(float latcell, float loncell, 
                   float WIND[WINDSTATIONS][15],   float CELLWIND[12])
{
  int i;
  float distance;
  float mindist1,mindist2; /* two smallest distances */
  int nr1,nr2; /* station number */

 /* Calculate distance to all stations */
  for(i=0;i<WINDSTATIONS;i++) 
    {
    distance = calc_distance(latcell,loncell,WIND[i][0], WIND[i][1]);
    WIND[i][14] = distance;
    }

  mindist1=mindist2 = 10000000;
  /* Find first the smallest distance. */
  /* Eliminate this station and then find the next smallest distance */
  for(i=0;i<WINDSTATIONS;i++) 
    {
    if (WIND[i][14] < mindist1) 
      {
      mindist1 = WIND[i][14]; 
      nr1 = i;
      }
    }

  /* eliminate the smallest one and find the next smallest */
  for(i=0;i<WINDSTATIONS;i++) 
    {
    if (i != nr1) 
      {
      if (WIND[i][14] < mindist2) 
        {
        mindist2 = WIND[i][14]; 
        nr2 = i;
        }
      }
    }
  /*  printf("mindist1 = %6.3f mindist2 = %6.3f\n", mindist1, mindist2);*/
  /*  printf("nr1 = %d  nr2 %d\n",nr1,nr2);*/

  for (i=0;i<12;i++) CELLWIND[i] = 0.0; /* Setting initial values */
  /* Calculate final wind value with simpel distance weighing */
  for (i=0;i<12;i++)
    {
    CELLWIND[i] = WIND[nr1][i+2] * (mindist1/(mindist1+mindist2))
                + WIND[nr2][i+2] * (mindist2/(mindist1+mindist2));
    }
}/* END function void sort_stations() */


void make_windseries( float METDATA[STEPS][4], float CELLWIND[12])
{
  int i,j,k;
  int total_steps=0;
  printf("make_windseries\n");
  for(i=STARTYEAR; i < (STOPYEAR+1); i++)
    {
    for(j=0;j<12;j++)
      {
      for (k=0;k < days_in_month(i,j+1); k++) 
         { 
         METDATA[total_steps][3] = CELLWIND[j]; 
         total_steps++; 
         } 
      } /*   for(j=0;j<12;j++)  */
    }/* END for(i=STARTYEAR;i<STOPYEAR+1; i++) */
  printf("total_steps %d\n",total_steps);
}/* END function void make_windseries() */

void write_data(char file[200], float METDATA[STEPS][4])
{
  FILE *fp;
  int i;
  if((fp = fopen(file,"w"))==NULL) 
     {printf("Cannot open file %s \n",file);exit(0);} 
  for (i=0;i<STEPS;i++) 
    {
    fprintf(fp,"%6.2f %6.2f %6.2f %6.2f\n",
            METDATA[i][0],METDATA[i][1],METDATA[i][2],METDATA[i][3]);
    }
  fclose(fp);
}/* END function void write_data() */


