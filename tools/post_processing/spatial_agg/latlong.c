/* File:        latlong.c                                            */
/* Programmer:  Bernt Viggo Matheussen                               */
/* Date:        01.17.98                                             */
/* Version:     1.0                                                  */

/* This program reads an ArcInfo asciigrid file and generates        */
/* and prints out latitude and longitude and the cell value          */
/* for each valid gridcell.                                          */
/* The longitude and latitude will be used to generate input files   */
/* to the VIC model.                                                 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/* Welcome and introduction */
void welcome(void);


void main(int argc, char **argv)
{ 

  FILE *fp;                                     /* Filepointers for to files*/
  char maskfile[400];                /* Declaration of variables */
  char junk[40];
  int i,j;
  int rows,cols;
  int void_nr;
  float xllcorner,yllcorner,resolution;
  float uplftlat, uplftlon;
  float lat,lon;
  float value;

  if (argc<2) {  /* Must be exactly 1 argument behind the program name */
  welcome();
  printf("\nNot correct number of commandline arguments \n");
  printf("Usage:   latlong inputfile > newfile \n");
  exit(EXIT_FAILURE); }

  welcome();

  strcpy(maskfile,argv[1]);  

  if((fp = fopen(maskfile,"r"))==NULL){        /* Opens maskfile */
    printf("Cannot open file maskfile %s \n",maskfile);
    exit(0);}

  /* Read header and assign values to rows,cols,resolution, etc  */
  fscanf(fp,"%s",junk);   fscanf(fp,"%d",&cols);
  fprintf(stderr,"%-15s %d\n",junk,cols);
  fscanf(fp,"%s",junk);   fscanf(fp,"%d",&rows);
  fprintf(stderr,"%-15s %d\n",junk,rows);
  fscanf(fp,"%s",junk);        fprintf(stderr,"%-15s",junk);
  fscanf(fp,"%f",&xllcorner);  fprintf(stderr,"%.4f\n",xllcorner);
  fscanf(fp,"%s",junk);        fprintf(stderr,"%-15s",junk);
  fscanf(fp,"%f",&yllcorner);  fprintf(stderr,"%.4f\n",yllcorner);
  fscanf(fp,"%s",junk);        fprintf(stderr,"%-15s",junk);
  fscanf(fp,"%f",&resolution); fprintf(stderr,"%.4f\n",resolution);
  fscanf(fp,"%s",junk);   fscanf(fp,"%d",&void_nr);
  fprintf(stderr,"%-15s %d\n",junk,void_nr);

  /* Calculate upper left corner (cell centered) */
  lon = uplftlon = xllcorner + (resolution/2);
  lat = uplftlat = yllcorner + (rows*resolution) - (resolution/2);
  
  /* Read file and find valid cells. Print lat lon if cell is valid */
  for(i=0;i<rows;i++)                 /* Write lat long to latlong.txt */ 
    { 
    lon = uplftlon;   
    for(j=0;j<cols;j++) 
      { 
      fscanf(fp,"%f",&value); 
      if (value > 0.0) printf("%.4f %.4f %.4f\n",lat,lon,value); 
      lon = lon + resolution; 
      } 
    lat = lat - resolution; 
    } 

  fclose(fp);
}/* END main */


void welcome(void)
{
  fprintf(stderr,"\nWelcome to latlong.c                           \n");
  fprintf(stderr,"This program reads a mask file and generates     \n");
  fprintf(stderr,"latitude and longitude for valid cells.          \n");
  fprintf(stderr,"The inputfile have to be an ArcInfo asciigridfile \n");
  fprintf(stderr,"NODATA_value have to be an integer,               \n");
  fprintf(stderr,"-9999, -99 or 0, etc..                            \n");
}/* END function void welcome() */

