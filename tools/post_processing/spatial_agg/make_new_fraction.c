/*
 * SUMMARY:      This program reads a fraction file and produces a new fraction
                 file with the double resolution i.e. 1/8 -> 1/4,  1/4 -> 1/2, 
 * USAGE:        make_new_fraction <fraction file>
 *
 * AUTHOR:       Bernt Viggo Matheussen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       tempbvm@hydro.washington.edu
 * ORIG-DATE:    23-Apr-99 at 09:57:10
 * LAST-MOD: Fri May 28 15:13:39 1999 by Keith Aric Cherkauer <cherkaue@u.washington.edu>
 * DESCRIPTION:   Used to produce a new VIC model at all resolutions from the 
 *               1/8th degree model in a consistant way.
 * DESCRIP-END.
 * FUNCTIONS:    
 * COMMENTS:     5-28-99 Added fraction file name to command line KAC
 */

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <float.h>
#include <math.h>
#include <string.h>


#define MAX 500

/* Reads fraction data into array */
void read_fraction(char file[400], float FR1[MAX][MAX], int rows, int cols);

/* Makes new maskfile with twice the resolution */
void make_mask(float FR1[MAX][MAX], float NEW[MAX][MAX], int rows, int cols);


void main(int argc, char *argv[])
{
  FILE *fp;
  char frac1[400];

  float FR1[MAX][MAX];
  float NEW[MAX][MAX];


  int rows,cols,void_nr;
  int i,j; 

  char junk[50];

  float xllcorner,yllcorner,resolution;

  if(argc!=2) {
    fprintf(stderr,"Usage: %s <fraction file>\n",argv[0]);
    fprintf(stderr,"\t<fraction file> is the routing fraction file for the current model resolution WITH AN ARC/INFO HEADER.\n");
    fprintf(stderr,"\tOutput is to stdout.\n");
    exit(0);
  }

  strcpy(frac1,argv[1]);

  if((fp = fopen(frac1,"r"))==NULL){        /* Opens maskfile */
    printf("Cannot open file maskfile %s \n",frac1);
    exit(0);}
  fscanf(fp,"%s %d",junk,&cols);       printf("%s      %d\n",junk,cols/2);     
  fscanf(fp,"%s %d",junk,&rows);       printf("%s      %d\n",junk,rows/2);
  fscanf(fp,"%s %f",junk,&xllcorner);  printf("%s      %.4f\n",junk,xllcorner);
  fscanf(fp,"%s %f",junk,&yllcorner);  printf("%s      %.4f\n",junk,yllcorner);
  fscanf(fp,"%s %f",junk,&resolution); printf("%s      %.4f\n",junk,resolution*2.0); 
  fscanf(fp,"%s %d",junk,&void_nr);    printf("%s      %d\n",junk,void_nr);
  fclose(fp);

  read_fraction(frac1,FR1,rows,cols); /* Reads fractionfile */

  make_mask(FR1,NEW,rows,cols);  /* Makes new mask with new resolution */ 


  for(i=0;i<(rows/2);i++)
    {
    for(j=0;j<(cols/2); j++) printf("%.3f ",NEW[i][j]);
    printf("\n");
    }



}/* END MAIN ******************************************************************/


void read_fraction(char file[400], float FR1[MAX][MAX], int rows, int cols)
{
  FILE *fp;
  int i,j;
  char junk[50];
  if((fp = fopen(file,"r"))==NULL){        /* Opens maskfile */
    printf("Cannot open file maskfile %s \n",file);
    exit(0);}
  for(j=0;j<12;j++)  fscanf(fp,"%s",junk); /* Scipping header */
  for(i=0;i<rows;i++) for(j=0;j<cols;j++)  fscanf(fp,"%f",&FR1[i][j]);
  fclose(fp);
}




void make_mask(float FR1[MAX][MAX], float NEW[MAX][MAX], int rows, int cols)
{
  int i,j;
  int AI,BI,CI,DI;
  int AJ,BJ,CJ,DJ;
  int okA,okB,okC,okD;
  float frA,frB,frC,frD;
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
      frA = frB = frC = frD = 0.0;
      okA = okB = okC = okD = 0;

/*       printf("%d %d, %d %d, %d %d, %d %d\n",AI,AJ, BI,BJ, CI,CJ, DI,DJ); */
      if(FR1[AI][AJ] > 0)  {okA = 1; frA = FR1[AI][AJ];}
      if(FR1[BI][BJ] > 0)  {okB = 1; frB = FR1[BI][BJ];}
      if(FR1[CI][CJ] > 0)  {okC = 1; frC = FR1[CI][CJ];}
      if(FR1[DI][DJ] > 0)  {okD = 1; frD = FR1[DI][DJ];}
      OK = okA + okB + okC + okD;
      if(OK > 0) 
        {
        NEW[i][j] = FR1[AI][AJ]*0.25 + FR1[BI][BJ]*0.25 + FR1[CI][CJ]*0.25 + FR1[DI][DJ]*0.25;
        }

      }
    }
}/* END function */

