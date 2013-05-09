/*
 * SUMMARY:      This program reads a direction file and produces a new direction
                 file with the double resolution i.e. 1/8 -> 1/4,  1/4 -> 1/2, 
 * USAGE:        Used to produce a new VIC model at all resolutions from the 
 *               1/8th degree model in a consistant way. The accumulation file is also needed.
 *
 * AUTHOR:       Bernt Viggo Matheussen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       tempbvm@hydro.washington.edu
 * ORIG-DATE:    23-Apr-99 at 09:57:10
 * LAST-MOD: Fri Jul  2 16:15:25 1999 by Keith Aric Cherkauer <cherkaue@u.washington.edu>
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


#define MAX 500

/* Reads data into array */
void read_dir(char file[400], int DIR1[MAX][MAX], int rows, int cols);

/* Makes new maskfile with twice the resolution */
void make_mask(int DIR1[MAX][MAX], int NEW[MAX][MAX], int rows, int cols);

/* Counts lines in file */
int line_count(char file[400]);

/* Counts number of coloumns in file */
int cols_count(char file[400], int rows);

void main(int argc, char **argv)
{
  FILE *fp;
  char dir1[400];
  char acc1[400]; 
  char junk[50];

  int DIR1[MAX][MAX];
  int ACC1[MAX][MAX];
  int NEW[MAX][MAX]; 
  

  int rows; 
  int cols; 

  int i,j; 
  int AI,BI,CI,DI; 
  int AJ,BJ,CJ,DJ; 
  int max_acc;
  int new_dir;

  float xllcorner,yllcorner,resolution;
  int void_nr;

  if (argc!=3) { /* Must be exactly 2 arguments behind the program name */
    fprintf(stderr,"Usage: %s <direction file> <accumulation file>\n",argv[0]);
    fprintf(stderr,"\tThis program builds a new routing direction file at twice the resolution (e.g. 1/8 to 1/4) of the provided direction file.  The accumulation file for the current direction file must also be provided.\n");
    fprintf(stderr,"\tThe new direction file is output to stdout.\n");
    exit(0); 
  }

  strcpy(dir1,argv[1]);
  strcpy(acc1,argv[2]);

   /*Read the header of the direction file */
  if((fp = fopen(dir1,"r"))==NULL){        /* Opens maskfile */
    printf("Cannot open file %s \n",dir1);
    exit(0);
  }
  
  fscanf(fp,"%s %d",junk,&cols);       
  printf("%-15s %d\n",junk,cols/2);
  fscanf(fp,"%s %d",junk,&rows);       
  printf("%-15s %d\n",junk,rows/2);
  fscanf(fp,"%s %f",junk,&xllcorner);  
  printf("%-15s %.4f\n",junk,xllcorner);
  fscanf(fp,"%s %f",junk,&yllcorner);  
  printf("%-15s %.4f\n",junk,yllcorner);
  fscanf(fp,"%s %f",junk,&resolution); 
  printf("%-15s %.4f\n",junk,resolution*2.0);
  fscanf(fp,"%s %d",junk,&void_nr);    
  printf("%-15s %d\n",junk,void_nr);
  fclose(fp);

  read_dir(dir1,DIR1,rows,cols); /* Reads direction file */
  read_dir(acc1,ACC1,rows,cols);  /* reads accumulation file */
  
  make_mask(DIR1,NEW,rows,cols); /* Makes new mask with new resolution */

  for(i=0;i<(rows/2);i++) {
    for(j=0;j<(cols/2);j++) {
      AJ = j*2;
      AI = i*2;
      BJ = j*2 + 1;
      BI = i*2;
      CJ = j*2;
      CI = i*2 + 1;
      DJ = j*2 + 1;
      DI = i*2 + 1;
      max_acc = new_dir = 0;

      if(NEW[i][j] > 0) {
        if(DIR1[AI][AJ] == 2) DIR1[AI][AJ] = 1; 
        if(DIR1[AI][AJ] == 6) DIR1[AI][AJ] = 7; 
        if(DIR1[BI][BJ] == 4) DIR1[BI][BJ] = 3; 
        if(DIR1[BI][BJ] == 8) DIR1[BI][BJ] = 1; 
        if(DIR1[CI][CJ] == 4) DIR1[CI][CJ] = 5; 
        if(DIR1[CI][CJ] == 8) DIR1[CI][CJ] = 7; 
        if(DIR1[DI][DJ] == 2) DIR1[DI][DJ] = 3; 
        if(DIR1[DI][DJ] == 6) DIR1[DI][DJ] = 5; 
	
        if(ACC1[AI][AJ] > max_acc) {
	  max_acc = ACC1[AI][AJ]; 
	  new_dir = DIR1[AI][AJ];
	}
        if(ACC1[BI][BJ] > max_acc) {
	  max_acc = ACC1[BI][BJ]; 
	  new_dir = DIR1[BI][BJ];
	}
        if(ACC1[CI][CJ] > max_acc) {
	  max_acc = ACC1[CI][CJ]; 
	  new_dir = DIR1[CI][CJ];
	}
        if(ACC1[DI][DJ] > max_acc) {
	  max_acc = ACC1[DI][DJ]; 
	  new_dir = DIR1[DI][DJ];
	}

/*         printf("DIR1[%d][%d] %d  ACC1[%d][%d] %d\n",AI,AJ,DIR1[AI][AJ],AI,AJ,ACC1[AI][AJ]); */
/*         printf("DIR1[%d][%d] %d  ACC1[%d][%d] %d\n",BI,BJ,DIR1[BI][BJ],BI,BJ,ACC1[BI][BJ]); */
/*         printf("DIR1[%d][%d] %d  ACC1[%d][%d] %d\n",CI,CJ,DIR1[CI][CJ],CI,CJ,ACC1[CI][CJ]); */
/*         printf("DIR1[%d][%d] %d  ACC1[%d][%d] %d\n",DI,DJ,DIR1[DI][DJ],DI,DJ,ACC1[DI][DJ]); */
/*         printf("\n"); */

        /* Sets the new direction */
        NEW[i][j] = new_dir;
      }
    }
  }


/*    printf("ncols           %d\n",cols/2);  */
/*    printf("nrows           %d\n",rows/2);  */
/*    printf("xllcorner      \n");  */
/*    printf("yllcorner      \n");  */
/*    printf("cellsize       \n");  */
/*    printf("NODATA_value    0\n");  */

  for(i=0;i<(rows/2);i++) {  
    for(j=0;j<(cols/2); j++) printf("%d ",NEW[i][j]);  
    printf("\n");  
  }  

}/* END MAIN ******************************************************************/


void read_dir(char file[400], int DIR1[MAX][MAX], int rows, int cols)
{
  FILE *fp;
  int i,j;
  char junk[50];
  if((fp = fopen(file,"r"))==NULL){        /* Opens maskfile */
    printf("Cannot open file maskfile %s \n",file);
    exit(0);}
  for(i=0;i<12;i++) fscanf(fp,"%s",junk);
  for(i=0;i<rows;i++) for(j=0;j<cols;j++)  fscanf(fp,"%d",&DIR1[i][j]);
  fclose(fp);
}

void make_mask(int DIR1[MAX][MAX], int NEW[MAX][MAX], int rows, int cols)
{

  int i,j;
  int AI,BI,CI,DI;
  int AJ,BJ,CJ,DJ;
  int OK;

  for(i=0;i<(rows/2);i++) for(j=0;j<(cols/2); j++) NEW[i][j] = 0;

  for(i=0;i<(rows/2);i++)
    {
    for(j=0;j<(cols/2); j++)
      {
      OK = 0;

      AJ = j*2;
      AI = i*2;

      BJ = j*2 + 1;
      BI = i*2;

      CJ = j*2;
      CI = i*2 + 1;

      DJ = j*2 + 1;
      DI = i*2 + 1;

/*       printf("%d %d, %d %d, %d %d, %d %d\n",AI,AJ, BI,BJ, CI,CJ, DI,DJ); */
      if(DIR1[AI][AJ] > 0) OK = 1;
      if(DIR1[BI][BJ] > 0) OK = 1;
      if(DIR1[CI][CJ] > 0) OK = 1;
      if(DIR1[DI][DJ] > 0) OK = 1;

      if(OK == 1) NEW[i][j] = 1;
      }
    }
}/* END function */


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

int cols_count(char file[400], int rows)
{
  FILE *fp;
  int c=0;
  char junk[40];
  if((fp = fopen(file,"r"))==NULL){
    printf("Cannot open file %s \n",file);exit(0);}
  while(fscanf(fp,"%s",junk) !=EOF) c++;
  fclose(fp);
  return (c/rows);
}


