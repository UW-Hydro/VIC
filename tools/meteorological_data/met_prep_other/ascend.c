/* File:        ascend.c                                             */
/* Programmer:  Bernt Viggo Matheussen                               */
/* Date:        02.27.99                                             */
/* Version:     1.0                                                  */

/* This program reads a mask file and generates a new maskfile       */
/* with unique integers in each valid cell.                          */
/* The input files have to be ArcInfo asciigrids.                    */
/* NODATA_value have to be an integer, -9999, -99 or 0, etc..        */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


void welcome(void);


void main(int argc, char **argv)
{ 

  FILE *fp;         
  char maskfile[400];
  int i,j;
  int rows,cols;
  int void_nr;   /* NODATA_value */
  char junk[20];
  float value;
  int cells =1;

  if (argc<2) {  /* Must be exactly 1 argument behind the program name */
  welcome();
  printf("\nNot correct number of commandline arguments \n");
  printf("Usage:   ascend maskfile > newmaskfile\n");
  exit(EXIT_FAILURE); }


  welcome();


  strcpy(maskfile,argv[1]);  
  
  fprintf(stderr,"Inputfile:  %s\n",maskfile);

  if((fp = fopen(maskfile,"r"))==NULL){        /* Opens maskfile */
    printf("Cannot open file maskfile %s \n",maskfile);
    exit(0);}


  /* Read header of file */
  fscanf(fp,"%s %d",junk,&cols);
  printf("%-15s %d\n",junk,cols);
  fprintf(stderr,"%-15s %d\n",junk,cols);

  fscanf(fp,"%s %d",junk,&rows);  
  printf("%-15s %d\n",junk,rows);
  fprintf(stderr,"%-15s %d\n",junk,rows);


  fscanf(fp,"%s %f",junk,&value); 
  printf("%-15s %.4f\n",junk,value);
  fprintf(stderr,"%-15s %.4f\n",junk,value);

  fscanf(fp,"%s %f",junk,&value); 
  printf("%-15s %.4f\n",junk,value);
  fprintf(stderr,"%-15s %.4f\n",junk,value);


  fscanf(fp,"%s %f",junk,&value); 
  printf("%-15s %.4f\n",junk,value);
  fprintf(stderr,"%-15s %.4f\n",junk,value);

  fscanf(fp,"%s %d",junk,&void_nr); 
  printf("%-15s %d\n",junk,void_nr);
  fprintf(stderr,"%-15s %d\n",junk,void_nr);

  for(i=0;i<rows;i++)                 /* Write lat long to latlong.txt */
    {
      for(j=0;j<cols;j++)
        {
        fscanf(fp,"%f",&value);
        if((int)(value) != void_nr) {printf("%d ",cells); cells++;}
        if((int)(value) == void_nr) printf("%d ",void_nr);

	}
    printf("\n");
    }
  fprintf(stderr,"Number of valid cells in mask is %d\n",(cells-1));
  fclose(fp);
} /* END main */


void welcome(void)
{
  fprintf(stderr,"\nWelcome to ascend.c                             \n");
  fprintf(stderr,"This program reads a mask file and generates a new\n"); 
  fprintf(stderr,"maskfile with unique integers in each valid cell. \n"); 
  fprintf(stderr,"The inputfile have to be an ArcInfo asciigridfile \n");
  fprintf(stderr,"NODATA_value have to be an integer,               \n");
  fprintf(stderr,"-9999, -99 or 0, etc..                            \n");
  fprintf(stderr,"Usage:    ascend maskfile > newmaskfile           \n");
}/* END function void welcome() */











