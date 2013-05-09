/* File:             get_prism.c                                                    */
/* Programmer :      Bernt Viggo Matheussen                                         */
/* Date:             02.10.98                                                       */
/* Version:          1.0                                                            */

/* This program reads the prism data (C.Daly et al 1993) and converts them into     */
/* monthly gridded data at a chosen resolution. The program needs a maskfile */
/* as input together with the prism rawdata.                                        */

/* The rawdata used in this program is located at /can-usa-prism/can_us_prism.jan   */
/* These files contains prism data for Canada and USA. If other prism rawdatafiles  */
/* are used, then lines 28 to 38 should be changed in this program                  */

/* If there are gridcells forexample in Canada that the prism data is not covering  */
/* The no_data value will be written in the outputfiles corresponding gridcell      */
/* The rescaling program will check for this void value at a later stage and the    */
/* precip value will not be scaled for this gridcell                                */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


/* IF OTHER PRISM RAWDATA ARE USED DOUBLECHECK THE NEXT LINES OF CODE  ****/
/**************************************************************************/
/*** Header of prism_us.mnth files, must be changed if header changes   ***/
#define NORTH_P (53.979166)          /*   north: 53:58:45N   */
                                     /*   south: 24:03:45N   */
                                     /*   east: 64:58:45W    */
#define WEST_P (126.020833)          /*   west: 126:01:15W   */
#define ROWS_P (718)                 /*   rows: 718          */
#define COLS_P (1465)                /*   cols: 1465         */
#define RES_P  (0.041666)            /*   1/24               */
/***************************************************************************/

int one_cell(int big[(int)ROWS_P][(int)COLS_P], float res_m, float res_p, float
	     lat, float lon, float N, float W, float void_nr);

void main(int argc, char **argv) 
 { 
  FILE *fpmask, *fpprism,*fpout;        
  char maskfile[400];     /* Total length of directory and filename < 400 */
  char prism_in[400]; 
  char out[400];
  char outfile[400];
  char str1[40];
  char void_nr_mask[10];
  char cdum[400];
  /*         rows   cols  */
  int all_data[(int)ROWS_P][(int)COLS_P];

  float lat_mask, long_mask, res_mask,north_prsm,west_prsm, res_prsm, north_mask;
  int  rows_mask, cols_mask, rows_prsm, cols_prsm;
  int i,k,j;
  int rows, cols;
  float elevation, latitude,longitude,void_nr, high, low;

  if (argc<4) {           /* Must be exactly 3 arguments behind the program name */
  printf("Not correct number of commandline arguments \n");
  exit(EXIT_FAILURE); }
  strcpy(maskfile,argv[1]);  printf("%s\n",maskfile);
  strcpy(prism_in,argv[2]);  printf("%s\n",prism_in);
  strcpy(outfile,argv[3]);   printf("%s\n",outfile);

  if((fpmask = fopen(maskfile,"r"))==NULL){
    printf("Cannot open file %s \n",maskfile);exit(0);}

  if((fpprism = fopen(prism_in,"r"))==NULL){
    printf("Cannot open file %s \n",prism_in);exit(0);}

  fpout = fopen(outfile,"w");

  /* Read header of maskfile lower left corner, not center */    
   fgets(str1, 40, fpmask);
   sscanf(str1,"%s %d", &cdum, &cols_mask);
   fprintf(fpout,"%s",str1);

   fgets(str1, 40 , fpmask);
   sscanf(str1,"%s %d", &cdum, &rows_mask);
   fprintf(fpout,"%s",str1);

   fgets(str1, 40 , fpmask);
   sscanf(str1,"%s %f", &cdum, &long_mask);
   fprintf(fpout,"%s",str1);
   long_mask *= -1.0;

   fgets(str1, 40 , fpmask);
   sscanf(str1,"%s %f", &cdum, &lat_mask);
   fprintf(fpout,"%s",str1);

   fgets(str1, 40 , fpmask);
   sscanf(str1,"%s %f", &cdum, &res_mask);
   fprintf(fpout,"%s",str1);

   fgets(str1, 40 , fpmask);
   sscanf(str1,"%s %s", &cdum, &void_nr_mask);
   fprintf(fpout,"%s",str1);
   void_nr = atof(void_nr_mask);
   low = void_nr - 0.001;
   high = void_nr + 0.001;

   /* Header of prism_us.mnth files, must be changed if header changes  */
   north_prsm = (float)NORTH_P;                  /*   north: 53:58:45N   */
                                                 /*   south: 24:03:45N   */
                                                 /*   east: 64:58:45W    */
   west_prsm = (float)WEST_P;                    /*   west: 126:01:15W   */
   rows_prsm = (int)ROWS_P;                      /*   rows: 718          */
   cols_prsm = (int)COLS_P;                      /*   cols: 1465         */
   res_prsm  = (float)RES_P;                     /*   1/24               */

   printf("north_prsm = %4.5f \n",north_prsm);
   printf("west_prsm = %4.5f \n",west_prsm);
   printf("rows_prsm = %d \n",rows_prsm);
   printf("cols_prsm = %d \n",cols_prsm);
   printf("res_prsm  = %4.5f \n",res_prsm);
   printf("res_mask %f \n",res_mask);
   north_mask = lat_mask + (res_mask*rows_mask) - (res_mask/2);
   printf("north_mask = %4.4f \n",north_mask);
   long_mask = long_mask- (res_mask/2);
   printf("long_mask = %4.4f\n",long_mask);


  for (i=0;i<12;i++) fscanf(fpprism,"%s",&str1);/*   moving down to data in prism.mnth */

  for (k=0;k<(int)ROWS_P;k++){ /** Sets all the data in big array */
    for (i=0;i<(int)COLS_P;i++)fscanf(fpprism,"%d",&all_data[k][i]); }
  

  latitude = north_mask;
  longitude = long_mask;

  for (j=0;j<rows_mask;j++)
    {
    for (i=0;i<cols_mask;i++)
      {
      fscanf(fpmask,"%f",&elevation);
      /*if (elevation>0.01)*/
      if(elevation>high || elevation<low)
        {
        printf("%d ",one_cell(all_data,res_mask,res_prsm,latitude,longitude,north_prsm,west_prsm, void_nr)); 
        fprintf(fpout,"%d ",one_cell(all_data,res_mask,res_prsm,latitude,longitude,north_prsm,west_prsm,void_nr));
        } 
      /*if (elevation <0.01)*/
      /*if(elevation<high && elevation>low)*/
      else
        {
        fprintf(fpout,"%.0f ",void_nr);
        printf("%.0f ",void_nr);
        }
      longitude = longitude - (res_mask);
      }
    longitude = long_mask;
    latitude = latitude - (res_mask);
    printf("\n");
    fprintf(fpout,"\n");
    }

   fclose(fpout);
   fclose(fpprism);
   fclose(fpmask);

   /* Header of prism_us.mnth files, must be changed if header changes  */
   printf("NORTH_P  %4.6f\n",NORTH_P);
   printf("WEST_P   %4.6f\n",WEST_P);
   printf("ROWS_P   %d\n",ROWS_P);
   printf("COLS_P   %d\n",COLS_P);
   printf("RES_P    %4.6f\n",RES_P);

    
}/*** END of main *************************************************************************/


int one_cell(int big[(int)ROWS_P][(int)COLS_P], float res_m, float res_p, float
	     lat, float lon, float N, float W, float void_nr)
{
  /* This function returns the average value of the prism values in one maskfile gridcell */
  int ll, cc;
  int val, i,j, nr;
  int ret;
  int valid;

  lat = lat + ((res_m/2.0)-(res_p/2.0));
  lon = lon + ((res_m/2.0)-(res_p/2.0));

  nr = (int)(res_m/res_p);
  ll = (int)((N-lat)/res_p);
  cc = (int)((W-lon)/res_p);
  val = 0;
  valid=0;
  for (i=ll;i<(ll+nr);i++) 
    {
    for(j=cc;j<(cc+nr);j++) 
     {
      if (big[i][j]>0){
	val += big[i][j]; 
	valid++;
      }
     }
    }
  if (val<1) ret= (int) void_nr;
  if (val>=1) ret = (int)(val/(valid));
  return ret;
}/*** END function */
