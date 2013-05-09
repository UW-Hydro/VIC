/*
 * SUMMARY:      Calculates the accumulated flow from a direction file 
 *
 * AUTHOR:       Bernt Viggo Matheussen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       tempbvm@hydro.washington.edu
 * ORIG-DATE:    22-Apr-99 at 10:22:02
 * LAST-MOD: Fri Jul  2 17:53:19 1999 by Keith Aric Cherkauer <cherkaue@u.washington.edu>
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


#define MAX       500  /* maximum number of rows/columns in direction file */
#define MIN_ACCUM 0    /* minimum accumulation value to print */


/* Reads in the outputfile from the flowgen program */
int read_file(char filename[400], int DATA[MAX][MAX], int rows, int cols);

/* Makes an accumulation file */
void make_accumulation(int DATA[MAX][MAX], int ACC[MAX][MAX], int rows,int cols);

/* Writes the accumulation data to screen. */
void write_accumulation(int ACC[MAX][MAX],int rows,int cols);

void main(int argc, char **argv) { 
  FILE *fp;
  char direc[400];                /* Declaration of variables */
  int DATA[MAX][MAX]; /* Contains the flownetwork */
  int ACC[MAX][MAX];  /* Contains the accumulated values from the flownetwork */
  int cols, rows,i,j;
  float uplftlat, uplftlon;
  float lat,lon,xllcorner,yllcorner;
  float resolution;
  int valid_cells;
  char junk[40];
  int void_nr;

  if (argc!=2) {           /* Must be exactly 1 arguments behind the program name */
    fprintf(stderr,"Usage: %s <direction file>\n",argv[0]);
    fprintf(stderr,"\tThis program computes the flow accumulation for each grid cell based on the given routing direction file.\n");
    fprintf(stderr,"\tThe accumulation file is output to stdout.\n");
    exit(0); 
  }

  strcpy(direc,argv[1]);

  fprintf(stderr,"Welcome to make_accumulation.c \n");

  if((fp = fopen(direc,"r"))==NULL) {        /* Opens maskfile */
    printf("Cannot open file maskfile %s \n",direc);
    exit(0);
  }

  /* Read header and assign values to rows,cols,resolution, etc  */
  fscanf(fp,"%s",junk);   fscanf(fp,"%d",&cols);
  printf("%-15s %d\n",junk,cols);

  fscanf(fp,"%s",junk);   fscanf(fp,"%d",&rows);
  printf("%-15s %d\n",junk,rows);

  fscanf(fp,"%s",junk);        printf("%-15s",junk);
  fscanf(fp,"%f",&xllcorner);  printf("%.4f\n",xllcorner);
  fscanf(fp,"%s",junk);        printf("%-15s",junk);
  fscanf(fp,"%f",&yllcorner);  printf("%.4f\n",yllcorner);
  fscanf(fp,"%s",junk);        printf("%-15s",junk);
  fscanf(fp,"%f",&resolution); printf("%.4f\n",resolution);
  fscanf(fp,"%s",junk);   fscanf(fp,"%d",&void_nr);
  printf("%-15s %d\n",junk,void_nr);


  /* Calculate upper left corner (cell centered) */
  lon = uplftlon = xllcorner + (resolution/2);
  lat = uplftlat = yllcorner + (rows*resolution) - (resolution/2);
  fprintf(stderr,"uplftlon %.4f uplftlat %.4f\n",uplftlon,uplftlat);

  fclose(fp);


  /* Set initial data */
  for(i=0;i<rows+2;i++)   for(j=0;j<cols+2;j++) DATA[i][j] = 9;
  for(i=0;i<rows+2;i++)   for(j=0;j<cols+2;j++) ACC[i][j] = 0;


  fprintf(stderr,"READING FILE\n");
  valid_cells = read_file(direc,DATA,rows,cols);

  fprintf(stderr,"valid_cells = %d\n",valid_cells);

  make_accumulation(DATA,ACC,rows,cols); 

  write_accumulation(ACC,rows,cols);  


}/* END main  ********************************************************************/

int read_file(char filename[400], int DATA[MAX][MAX], int rows, int cols)
{
  FILE *fp;                /* Filepointers for to files*/
  int i,j;
  int cells=0;
  int value;
  int nines=0;
  char junk[50];
  if((fp = fopen(filename,"r"))==NULL)
    {printf("Cannot open file  %s \n",filename); exit(0);}

  /* Putting in one row of 9 around the whole grid */
  /* Example                                       */
  /* 335                                           */
  /* 135                                           */
  /* 335                                           */

  /* 99999                                         */
  /* 93359                                         */
  /* 91359                                         */
  /* 93359                                         */
  /* 99999                                         */ 
     
  /* This is to make the program not go into       */
  /* unexpected loops. It will always end          */

  /* Scipping header */
  for (j=0;j<12;j++) fscanf(fp,"%s",junk);

  for (j=0;j<cols+2;j++) DATA[0][j] = 9;      /* First row */
  for (j=0;j<cols+2;j++) DATA[rows+1][j] = 9; /* Last row */
  for (j=0;j<rows+2;j++) DATA[j][0] = 9; /* First coloumn */
  for (j=0;j<rows+2;j++) DATA[j][cols+1] = 9; /* Last coloumn */
  for(i=1;i<rows+1;i++) for(j=1;j<cols+1;j++)
    {
    fscanf(fp,"%d",&value);
    if(value > 0) 
      {
      DATA[i][j] = value;
      cells++;
      }
    if(value == 9) nines++;
    }

  fprintf(stderr,"nines %d\n",nines);
  fclose(fp);  
  return cells;
}/* END function *********************************************/

void make_accumulation(int DATA[MAX][MAX], int ACC[MAX][MAX], int rows,int cols)
{
  int i,j;
  int no_more_data=0;
  int valid_cell=0;
  int local_row, local_col;
  int new_row, new_col;
  int cells_in_flowpath=0;

  int LOOP[MAX][MAX];
  int k,l;

  for (i=1;i<rows+1;i++)
    {
    for(j=1;j<cols+1;j++)
      { /****** START MAIN LOOP *************/
      valid_cell = 0;
      if(DATA[i][j] == 9) valid_cell = 1; /* Checks to see if cell have value 9 */
      if(valid_cell == 0)
        {
        local_row = new_row = i;
        local_col = new_col = j;
        no_more_data=0; 
        cells_in_flowpath =0;

        for(k=0;k<rows+2;k++)   for(l=0;l<cols+2;l++) LOOP[k][l] = 0;

        while (no_more_data == 0) {
          LOOP[local_row][local_col] = LOOP[local_row][local_col] + 1; 
          ACC[local_row][local_col] = ACC[local_row][local_col] + 1; 
	  /* Adds 1 because the cell is not a 9 */
          cells_in_flowpath++;
          if (DATA[local_row][local_col] == 1) {new_row--;}
          if (DATA[local_row][local_col] == 2) {new_row--;new_col++;}
          if (DATA[local_row][local_col] == 3) {new_col++;}
          if (DATA[local_row][local_col] == 4) {new_row++;new_col++;}
          if (DATA[local_row][local_col] == 5) {new_row++;}
          if (DATA[local_row][local_col] == 6) {new_row++;new_col--;}
          if (DATA[local_row][local_col] == 7) {new_col--;}
          if (DATA[local_row][local_col] == 8) {new_row--;new_col--;}
        
          if(DATA[new_row][new_col] == 9) no_more_data = 1;

          if(LOOP[local_row][local_col] > 1) { 
            fprintf(stderr,"LOOP at local_row %d local_col %d\n",
		    local_row,local_col); 
            fprintf(stderr,"in gridcell i= %d j= %d\n",i,j); 
            fprintf(stderr,"DATA[i][j] = %d\n",DATA[i][j]); 
	    
            no_more_data = 1; 
            exit(0); /* Aborts the program */
	  } 
          local_row = new_row;
          local_col = new_col;
          }
        }
        /****** END MAIN LOOP ************/
      }/* for(j=0;j<cols+2;j++) */
    } /*  for (i=0;i<rows+2;i++) */
}/* END function make_accumulation() */

void write_accumulation(int ACC[MAX][MAX],int rows,int cols)
{
  int i,j;
  for(i=1;i<rows+1;i++) 
    {
    for(j=1;j<cols+1;j++) 
      {
      if(ACC[i][j] > MIN_ACCUM) printf("%d ",ACC[i][j]);
      else printf("0 ");
      }
    printf("\n"); 
    }
}/* END function void write_accumulation() */

