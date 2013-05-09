/* Programmer: Laura Bowling                                   
   Usage: elevband <soil file> <resolution> <elevation file> <outfile>   
   Compile with: gcc elevband.c -lm -o elevband
   Purpose: This program reads in a VIC-NL soil file, the resolution of 
   the VIC run,and a finer resolution elevation file (with arc/info header) 
   in geographic projection which encompasses the entire area.  It will
   output the VIC-NL snow band elevation file.                 

   Elevation band widths are calculated as 2*max(xmax-mean,
   mean-xmin)/NUMBANDS.  Elevations are then binned:
   bin1: mean - max(xmax-mean,mean-xmin) through
   mean - max(xmax-mean,mean-xmin)*((2-NUMBANDS)/NUMBANDS);
   bin2: mean - max(xmax-mean,mean-xmin)*((2-NUMBANDS)/NUMBANDS)
   through mean - max(xmax-mean,mean-xmin)*((4-NUMBANDS)/NUMBANDS);
   etc.
   
   The mean elevation and area fraction are then calculated for 
   each bin.

   Edited 8/10/98 to merge bands less than MINDELTA meters apart in elevation
   and eliminate bands of less than MINFRACTION fractional area. 

   Edited 8/16/98 to calculate precipitation in each band based on the 
   elevation of the band relative to the pixel mean.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAXSTRING 200
#define NUMBANDS 5.
#define MINFRACTION .01
/* #define MINDELTA 153.85 */
#define MINDELTA 38.46

int main(int argc, char *argv[])
{
  FILE *fo, *fdem, *fm;
  char demfile[80], soilfile[80], outfile[80];
  char tempstr[MAXSTRING];
  int i,j,m, i_corner, j_corner;
  int **dem;
  int columns, rows;
  float xmin, ymin, size, nodata;
  float sum, min, max, count;
  float mean, difference, *increment;
  float resolution;
  int *incsum, *inccount, flag;
  float lat, lon, avg, fract;
  int testcount, numcells, cellnum;
  float corner_lat, corner_lon;
  float pixelmean;
  int numremove =0;
  int numzero=0;

  if(argc != 5) {
    printf("Usage: elevband <soil file> <resolution> <elevation file> <outfile>\n\n");
    printf("\t\tsoil file: VIC-NL soil parameter file;\n");
    printf("\t\tresolution: (float) resolution of the VIC model run;\n");
    printf("\t\televation file: DEM with arcinfo header;\n");
    printf("\t\toutfile: the elevation band input file for VUC.\n\n");
    exit(0);
  }

  strcpy(soilfile,argv[1]);
  strcpy(demfile,argv[3]);
  strcpy(outfile,argv[4]);
  resolution = atof(argv[2]);
  
  /*--------------------------------------------*/
  /* OPEN FILES                                 */
  /*--------------------------------------------*/

  if((fdem = fopen(demfile,"r")) == NULL) {
    printf("Cannot open/read dem file, %s\n",demfile);
    exit(8); }

  if((fm = fopen(soilfile,"r")) == NULL) {
    printf("Cannot open/read soil file, %s\n",soilfile);
    exit(8); }

  if((fo = fopen(outfile, "w")) == NULL) {
    printf("Cannot open/read outfile, %s\n",outfile);
    exit(8); }
 
  printf("Files successfully opened...\n\n");

  /*--------------------------------------------*/
  /* Scan DEM header.                           */
  /*--------------------------------------------*/

  fscanf(fdem,"%s %d",tempstr,&columns);
  fscanf(fdem,"%s %d",tempstr,&rows);
  fscanf(fdem,"%s %f",tempstr,&xmin);
  fscanf(fdem,"%s %f",tempstr,&ymin);
  fscanf(fdem,"%s %f",tempstr,&size);
  fscanf(fdem,"%s %f",tempstr,&nodata);
  printf("%d %d %f %f\n",columns,rows,xmin,ymin);

  /* ALLOCATE MEMORY TO ARRAYS */

  printf("Allocating memory to arrays...\n");


  if(!(dem = (int**) calloc(rows,sizeof(int*)))) {
    printf("Cannot allocate memory to first record: BASIN\n");
    exit(8); }
  for(i=0; i<rows;i++) {
    if(!(dem[i] = (int*) calloc(columns,sizeof(int)))) {
      printf("Cannot allocate memory to first record: BASIN\n");
      exit(8); }
  }

  if(!(inccount = (int*) calloc(NUMBANDS,sizeof(int)))) {
    printf("Cannot allocate memory to first record: BASIN\n");
    exit(8); }

  if(!(incsum = (int*) calloc(NUMBANDS,sizeof(int)))) {
    printf("Cannot allocate memory to first record: BASIN\n");
    exit(8); }

  if(!(increment = (float*) calloc(NUMBANDS + 1,sizeof(float)))) {
    printf("Cannot allocate memory to first record: BASIN\n");
    exit(8); }


  printf("Done...\n\n");


  
  /*--------------------------------------------*/
  /* READ DEM FILES                         */
  /*--------------------------------------------*/

  printf("Reading in files\n");
  for(i=0; i<rows;i++) {
    for(j=0; j<columns; j++) {
      fscanf(fdem,"%d",&dem[i][j]);
      if(dem[i][j] < 0)
	dem[i][j] = nodata;
    }
  }

  fclose(fdem);

	
  /*--------------------------------------------*/
  /* FIND INITIAL DISTRIBUTION                  */
  /*--------------------------------------------*/

  fgets(tempstr,MAXSTRING,fm);
  
  numcells = (int)(resolution/size);
  while((fgets(tempstr,MAXSTRING,fm)) != NULL) {
    sscanf(tempstr,"%d %d %f %f", &flag, &cellnum, &lat, &lon);
    
    /* Calculate position in dem file. */
    corner_lat = lat + resolution/2.0;
    i_corner = (int)rows - ((corner_lat - ymin)/size) - 1;
    corner_lon = lon - resolution/2.0;
    j_corner = (int)(corner_lon - xmin)/size;
    
    min = 9999.;
    max = 0.;
    sum=0.;
    count =0.;
    for(i=i_corner; i<(i_corner+numcells);i++) {
      for(j=j_corner; j<(j_corner+numcells);j++) {

	if(dem[i][j] != nodata ) {
	  sum += dem[i][j];
	  count ++;
	  if(dem[i][j] < min)
	    min = dem[i][j];
	  if(dem[i][j] > max)
	    max = dem[i][j];
	}
      }
    }
    if(count==0)
      {
	printf("ERROR: No data in current cell. \n");
	printf("lat = %.2f, lon = %.2f\n",lat,lon);
      }
    else
      {
	mean = (sum/count);
	if((max-mean) > (mean-min))
	  difference = max-mean;
	else
	  difference = (mean-min);
	
	/*	printf("mean=%f max =%f min=%f \n",mean,max,min);*/
	/* Find boundaries of elevations bins. */
	increment[0] = mean - difference;
	if(min < increment[0]) {
	  printf("range higher than min\n");
	  exit(0);
	}
	
	for(m=1;m<NUMBANDS+1;m++) {
	  increment[m]=increment[m-1] + (2.*difference)/NUMBANDS;
	  inccount[m-1]=0;                
	  incsum[m-1]=0;
	}
	if(max>(increment[(int)NUMBANDS]+1)) {
	  printf("range less than max\n");
	  exit(0);
	}
	increment[(int)NUMBANDS]+=1;   /* To account for round off errors. */
	
	/*	printf("Begin allocating elevations to bins.\n");*/
	for(i=i_corner; i<(i_corner+numcells);i++) {
	  for(j=j_corner; j<(j_corner+numcells);j++) {
	    if(dem[i][j] != nodata) {
	      for(m=0;m<NUMBANDS;m++) {
		if(dem[i][j] >= increment[m] && dem[i][j] < increment[m+1]) {
		  inccount[m]++;
		  incsum[m] += dem[i][j];
		}
	      }
	    }
	  }
	}

	/*********************************************/
	/* Check bins. */
	/*********************************************/
	testcount = 0;
	for(m=0; m<NUMBANDS;m++)
	  testcount += inccount[m];
    
	if(count!=testcount) {
	  printf("Not all elevations allocated to bins.\n");
	  printf("Total values: %d\t Allocated values: %d\n",count,testcount);
	}
	
	/* Check for bins within MINDELTA. */

	for(m=1;m<NUMBANDS;m++) {
	  if(inccount[m]>0 && inccount[m-1] >0&&((incsum[m]/inccount[m]) - (incsum[m-1]/inccount[m-1])) < MINDELTA)
	    { /*Combine bins. m and m-1*/
	      inccount[m]+=inccount[m-1];
	      incsum[m]+=incsum[m-1];
	      inccount[m-1]=0;
	      incsum[m-1]=0.0;
	      numremove += 1;
	    }
	}

	/* Check for bins less than MINFRACTION. */
	for(m=0;m<NUMBANDS-1;m++) {
	  fract=(float)inccount[m]/(float)count;
	  if(fract <= MINFRACTION && inccount[m] > 0.0 && inccount[m+1] >0.0) {
	    inccount[m+1] += inccount[m];
	    incsum[m+1] += incsum[m];
	    inccount[m]=0;
	    incsum[m]=0.0;
	    numremove +=1;
	  }
	}
	if(((float)inccount[(int)NUMBANDS-1]/(float)count) < MINFRACTION &&
	   inccount[(int)NUMBANDS-1]>0 && inccount[(int)NUMBANDS-2]>0) {
	  inccount[(int)NUMBANDS-2] += inccount[(int)NUMBANDS-1];
	  incsum[(int)NUMBANDS-2] += incsum[(int)NUMBANDS-2];
	  inccount[(int)NUMBANDS-1]=0;
	  incsum[(int)NUMBANDS-1]=0.0;
	  numremove +=1;
	}
	  
	
	/*--------------------------------------------*/
	/* OUTPUT RESULTS                      */
	/*--------------------------------------------*/
    
	/*	printf("Output...\n"); */

    
	fprintf(fo,"%d\t",cellnum);
	
	/* Output fractional area of each band. */
	for(m=0;m<(int)NUMBANDS;m++) {
	  fract = (float)inccount[m]/(float)count;
	  fprintf(fo,"%.3f\t",fract);
	}
    
	/* Output mean elevation of each band. */
	pixelmean=0.0;
	for(m=0;m<(int)NUMBANDS;m++) {
	  fract = (float)inccount[m]/(float)count;
	  if(inccount[m] > 0)
	    avg = (float)incsum[m]/(float)inccount[m];
	  else {
	    avg = 0.0;
	    numzero+=1;
	  }
	  pixelmean+=avg*fract;
	  fprintf(fo,"%.2f\t",avg);
	}

	/* Output fractional precip for each band. */
	for(m=0;m<(int)NUMBANDS;m++) {
	  avg = (float)incsum[m]/(float)inccount[m];
	  fract = (float)inccount[m]/(float)count;
	  if(fract >0.0)
	    fprintf(fo,"%.3f ",fract*avg/pixelmean);      
	  else
	    fprintf(fo,"0.000 ");
	}
	fprintf(fo,"\n");
      }
  }
  printf("Eliminated %d snow bands.\n",numremove);
  printf("Total of %d zero bands.\n",numzero);
  fclose(fo);
  fclose(fm);
}

