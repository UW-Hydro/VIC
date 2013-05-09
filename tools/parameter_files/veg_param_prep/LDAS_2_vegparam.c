#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <math.h>

#define MAXSTRING 1024

typedef struct {
  int   cellnum;
  int   Nveg;
  float fract[14];
} VEG_STRUCT;

float read_arcinfo_value(char *, float, float);
int read_arcinfo_info(char    *, float **, float **, int **, int *, int *,
		      float *, float *, float *, int *);
void     quick(VEG_STRUCT *item, int count);
void     qs(VEG_STRUCT *item,int left,int right);

main(int argc, char *argv[]) {
/**********************************************************************
  LDAS_2_vegparam.c          Keith Cherkauer        March 28, 1999

  This program was written to create a VIC vegetation parameter file
  (without root zones) from the ASCII file of UMD LDAS vegetation
  cover categories (14 in all).

  Category   Vegetation
  0          Water
  1          Evergreen needleleaf
  2          Evergreen broadleaf
  3          Deciduous needleleaf
  4          Deciduous broadleaf
  5          Mixed cover
  6          Woodland
  7          Wooded grasslands
  8          Closed shrublands
  9          Open shrublands
  10         Grassland
  11         Cropland
  12         Bare ground
  13         Urban

**********************************************************************/

  FILE  *fin, *fout;
  char   cellname[256];
  char   tmpstr[1024];
  float  lat, lng;
  int    Nfreq, freq[14];
  float *cell_lat;
  float *cell_lng;
  int   *cell_num;
  int    Ncell;
  int    cell;
  int    Nzones;
  int    zone;
  int    nonzero;
  int    i, j;
  int    ncols;
  int    nrows;
  float  ll_lng;
  float  ll_lat;
  float  cellsize;
  int    NODATA;
  VEG_STRUCT *veg_param;

  if(argc!=5) {
    fprintf(stderr,"Usage: %s <LDAS file> <ARCINFO cell number file> <vegparam file> <# root zones>\n",argv[0]);
    fprintf(stderr,"\tNOTE: This program adds generic roots zones.  The user is expected to adjust them before running VIC.\n");
    exit(0);
  }

  if((fin=fopen(argv[1],"r"))==NULL) {
    fprintf(stderr,"ERROR: Unable to open LDAS file %s.\n",argv[1]);
    exit(0);
  }

  strcpy(cellname,argv[2]);

  if((fout=fopen(argv[3],"w"))==NULL) {
    fprintf(stderr,"ERROR: Unable to open vegetation parameter file %s.\n",argv[3]);
    exit(0);
  }

  Nzones = atoi(argv[4]);

  Ncell = read_arcinfo_info(cellname,&cell_lat,&cell_lng,&cell_num,
			    &ncols,&nrows,&ll_lng,&ll_lat,&cellsize,
			    &NODATA);

  printf("Ncel = %i\n",Ncell);

  veg_param = (VEG_STRUCT *)calloc(Ncell,sizeof(VEG_STRUCT));

  fscanf(fin,"%*s %*s %f %f %i", &lat, &lng, &Nfreq);

  while(!feof(fin)) {

    for(i=0;i<14;i++) fscanf(fin,"%i", &freq[i]);

    for(cell=0;cell<Ncell;cell++) {
      if(lat == cell_lat[cell] && lng == cell_lng[cell]) {

	nonzero = 0;
	for(i=0;i<14;i++) 
	  if(freq[i] > 0) nonzero++;
	veg_param[cell].cellnum = cell_num[cell];
	veg_param[cell].Nveg = nonzero;
	for(i=0;i<14;i++) {
	  veg_param[cell].fract[i] = (float)freq[i]/(float)Nfreq;
	}
      }
    }
    fscanf(fin,"%*s %*s %f %f %i", &lat, &lng, &Nfreq);
  }

  quick(veg_param,Ncell);

  for(cell=0;cell<Ncell;cell++) {

    fprintf(fout,"%i\t%i\n",veg_param[cell].cellnum,veg_param[cell].Nveg);
    for(i=0;i<14;i++) {
      if(veg_param[cell].fract[i] > 0) {
	fprintf(fout,"\t%i\t%f",i,veg_param[cell].fract[i]);
	for(j=0;j<Nzones;j++) fprintf(fout,"\t0.50\t%f",1.0/(float)Nzones);
	fprintf(fout,"\n");
      }
    }
  }

}

float read_arcinfo_value(char *filename,
			  float lat,
			  float lng) {
/**********************************************************************
  read_arcinfo_value           Keith Cherkauer           May 5, 1998

  This subroutine reads a single data value from an ARC/INFO ASCII 
  output grid file.  The latitude and longitude of the center of the
  grid cell of interest is provided to this routine.

**********************************************************************/

  FILE   *farc;
  int    i, j;
  int    ncols;
  int    nrows;
  float ll_lat;
  float ll_lng;
  float cellsize;
  float NODATA;
  float value;
  float tmp_lat;
  float tmp_lng;
  float tmpvalue;
  char   errstr[MAXSTRING];

  if((farc=fopen(filename,"r"))==NULL) {
    fprintf(stderr,"Unable to open ARC/INFO soil file %s",filename);
    exit(0);
  }

  /***** Read ARC/INFO Header *****/
  fscanf(farc,"%*s %i",&ncols);
  fscanf(farc,"%*s %i",&nrows);
  fscanf(farc,"%*s %f",&ll_lng);
  fscanf(farc,"%*s %f",&ll_lat);
  fscanf(farc,"%*s %f",&cellsize);
  fscanf(farc,"%*s %f",&NODATA);

  /***** Check for Valid Location *****/
  if(lat<ll_lat || lat>ll_lat+cellsize*(float)nrows) {
    fprintf(stderr,"Given latitude %f does not fall within ARC file %s",
	    lat,filename);
    exit(0);
  }

  if(lng<ll_lng || lng>ll_lng+cellsize*(float)ncols) {
    fprintf(stderr,"Given longitude %f does not fall within ARC file %s",
	    lng,filename);
    exit(0);
  }

  for(j=0;j<nrows;j++) {
    tmp_lat = ll_lat+(float)(nrows-j-0.5)*cellsize;
    for(i=0;i<ncols;i++) {
      tmp_lng = ll_lng + (float)(i+0.5)*cellsize;
      fscanf(farc,"%f",&tmpvalue);
      if(((int)(tmp_lat*1.e4+0.5) == (int)(lat*1.e4+0.5)) 
	 && ((int)(tmp_lng*1.e4+0.5) == (int)(lng*1.e4+0.5)))
	value = tmpvalue;
    }
  }
  fclose(farc);

  if(value==NODATA) value = -999.;
  return value;

}

int read_arcinfo_info(char    *filename,
		      float **lat,
		      float **lng,
		      int   **cellnum,
		      int    *ncols,
		      int    *nrows,
		      float  *ll_lng,
		      float  *ll_lat,
		      float  *cellsize,
		      int    *NODATA) {
/**********************************************************************
  read_arcinfo_info           Keith Cherkauer           May 5, 1998

  This subroutine reads data from an ARC/INFO ascii grid and returns 
  the central latitude and longitude of all grid cells that are not 
  assigned the value of NODATA.  Routine returns the number of defined 
  grid cell.

**********************************************************************/

  FILE   *farc;
  int    i, j;
  float value;
  float tmp_lat;
  float tmp_lng;
  int    cell;
  int    Ncells;
  int    tmpvalue;

  if((farc=fopen(filename,"r"))==NULL) {
    fprintf(stderr,"ERROR: Unable to open LDAS file %s.\n",filename);
    exit(0);
  }

  /***** Read ARC/INFO Header *****/
  fscanf(farc,"%*s %i",&ncols[0]);
  fscanf(farc,"%*s %i",&nrows[0]);
  fscanf(farc,"%*s %f",&ll_lng[0]);
  fscanf(farc,"%*s %f",&ll_lat[0]);
  fscanf(farc,"%*s %f",&cellsize[0]);
  fscanf(farc,"%*s %i",&NODATA[0]);

  /***** Allocate Latitude and Longitude Arrays for maximum size *****/
  Ncells  = ncols[0]*nrows[0];
  lat[0]  = (float *)calloc(Ncells,sizeof(float));
  lng[0]  = (float *)calloc(Ncells,sizeof(float));
  cellnum[0] = (int *)   calloc(Ncells,sizeof(int));

  /***** Check for Valid Location *****/
  cell = 0;
  for(j=0;j<nrows[0];j++) {
    tmp_lat = ll_lat[0]+(float)(nrows[0]-j-0.5)*cellsize[0];
    for(i=0;i<ncols[0];i++) {
      tmp_lng = ll_lng[0] + (float)(i+0.5)*cellsize[0];
      fscanf(farc,"%i",&tmpvalue);
      if(tmpvalue!=NODATA[0]) {
	lat[0][cell]  = tmp_lat;
	lng[0][cell]  = tmp_lng;
	cellnum[0][cell] = tmpvalue;
	cell++;
      }
    }
  }
  fclose(farc);
  Ncells = cell;

  if(value==NODATA[0]) value = -999.;
  return Ncells;

}

void quick(VEG_STRUCT *item, int count)
/**********************************************************************
        this subroutine starts the quick sort
**********************************************************************/
{
  qs(item,0,count-1);
}
 
void qs(VEG_STRUCT *item,int left,int right)
/**********************************************************************
        this is the quick sort subroutine - it returns the values in
	an array from high to low.
**********************************************************************/
{
  register int i,j;
  VEG_STRUCT x,y;
 
  i=left;
  j=right;
  x=item[(left+right)/2];
 
  do {
    while(item[i].cellnum<x.cellnum && i<right) i++;
    while(x.cellnum<item[j].cellnum && j>left) j--;
 
    if (i<=j) {
      y=item[i];
      item[i]=item[j];
      item[j]=y;
      i++;
      j--;
    }
  } while (i<=j);
 
  if(left<j) qs(item,left,j);
  if(i<right) qs(item,i,right);
}

