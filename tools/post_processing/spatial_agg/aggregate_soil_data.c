#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define MAXSTR 512
#define D_PLACES 4

int read_arcinfo_grid(char    *, double **, double **, int     *,
		      int     *, double  *, double  *, double  *,
		      int     *, double **);
int write_arcinfo_grid(char   *, double *, double *, double *,
		       int     , int     , double  , double  ,
		       double  , int     , int     );

main(int argc, char *argv[]) {

  char lowresmask[MAXSTR];
  char infile[MAXSTR];
  char outfile[MAXSTR];

  double  store_values[4];
  double *hr_values, *hr_lat, *hr_lng;
  double *lr_values, *lr_lat, *lr_lng;
  double  hr_ll_lat, hr_ll_lng, hr_cellsize;
  double  lr_ll_lat, lr_ll_lng, lr_cellsize;
  int     hr_ncols, hr_nrows;
  int     lr_ncols, lr_nrows;
  int     NODATA, storecnt;
  int     Ncells_lr, Ncells_hr;
  int     i, j;

  if(argc!=4) {
    fprintf(stderr,"Usage: %s <high res grid> <low res mask> <low res output gird>\n",argv[0]);
    fprintf(stderr,"\t<high res grid> is the high resolution ARC/INFO grid.\n");
    fprintf(stderr,"\t<low res mask> is an ARC/INFO grid mask file at the lower resolution.\n");
    fprintf(stderr,"\tThis program aggregates the high resolution ARCINFO file to the resolution of the low resolution mask file (only aggregates by a factor of 2, e.g. 1/8 to 1/4, 1/4 to 1/2) and outputs a new ARCINFO grid file to <low res output grid>.\n");
    exit(0);
  }

  Ncells_lr = read_arcinfo_grid(argv[2],&lr_lat,&lr_lng,&lr_ncols,
				&lr_nrows,&lr_ll_lat,&lr_ll_lng,
				&lr_cellsize,&NODATA,&lr_values);

  Ncells_hr = read_arcinfo_grid(argv[1],&hr_lat,&hr_lng,&hr_ncols,
				&hr_nrows,&hr_ll_lat,&hr_ll_lng,
				&hr_cellsize,&NODATA,&hr_values);

  for(i=0;i<Ncells_lr;i++) {
    j = 0;
    storecnt = 0;
    while(j<Ncells_hr && storecnt<4) {
      if(rint(hr_lat[j]*pow(10,(double)D_PLACES)) 
	 == rint((lr_lat[i]+hr_cellsize/2)*pow(10,D_PLACES)) 
	 && rint(hr_lng[j]*pow(10,(double)D_PLACES)) 
	 == rint((lr_lng[i]+hr_cellsize/2)*pow(10,D_PLACES))) {
	if(hr_values[j] != NODATA) {
	  store_values[storecnt] = hr_values[j];
	  storecnt++;
	}
      }
      if(rint(hr_lat[j]*pow(10,(double)D_PLACES)) 
	 == rint((lr_lat[i]-hr_cellsize/2)*pow(10,D_PLACES)) 
	 && rint(hr_lng[j]*pow(10,(double)D_PLACES)) 
	 == rint((lr_lng[i]+hr_cellsize/2)*pow(10,D_PLACES))) {
	if(hr_values[j] != NODATA) {
	  store_values[storecnt] = hr_values[j];
	  storecnt++;
	}
      }
      if(rint(hr_lat[j]*pow(10,(double)D_PLACES)) 
	 == rint((lr_lat[i]+hr_cellsize/2)*pow(10,D_PLACES)) 
	 && rint(hr_lng[j]*pow(10,(double)D_PLACES)) 
	 == rint((lr_lng[i]-hr_cellsize/2)*pow(10,D_PLACES))) {
	if(hr_values[j] != NODATA) {
	  store_values[storecnt] = hr_values[j];
	  storecnt++;
	}
      }
      if(rint(hr_lat[j]*pow(10,(double)D_PLACES)) 
	 == rint((lr_lat[i]-hr_cellsize/2)*pow(10,D_PLACES)) 
	 && rint(hr_lng[j]*pow(10,(double)D_PLACES)) 
	 == rint((lr_lng[i]-hr_cellsize/2)*pow(10,D_PLACES))) {
	if(hr_values[j] != NODATA) {
	  store_values[storecnt] = hr_values[j];
	  storecnt++;
	}
      }
      j++;
    }
    if(storecnt>0) {
      lr_values[i] = 0;
      for(j=0;j<storecnt;j++) lr_values[i] += store_values[j];
      lr_values[i] /= (double)storecnt;
    }
    else lr_values[i] = NODATA;
  }

  Ncells_hr = write_arcinfo_grid(argv[3],lr_values,lr_lat,lr_lng,
				 lr_ncols,lr_nrows,lr_ll_lat,lr_ll_lng,
				 lr_cellsize,NODATA,Ncells_lr);

}


int read_arcinfo_grid(char    *filename,
		      double **lat,
		      double **lng,
		      int     *ncols,
		      int     *nrows,
		      double  *ll_lat,
		      double  *ll_lng,
		      double  *cellsize,
		      int     *NODATA,
		      double **values) {
/**********************************************************************
  read_arcinfo_info           Keith Cherkauer           May 5, 1998

  This subroutine reads data from an ARC/INFO ascii grid and returns 
  the central latitude and longitude of all grid cells that are not 
  assigned the value of NODATA.  Routine returns the number of defined 
  grid cell.

**********************************************************************/

  FILE   *farc;
  int    i, j;
  double tmp_lat;
  double tmp_lng;
  int    cell;
  int    Ncells;
  double tmpvalue;

  if((farc=fopen(filename,"r"))==NULL) {
    fprintf(stderr,"ERROR: Unable of open ARCINFO grid %s for reading.\n",filename);
    exit(0);
  }

  /***** Read ARC/INFO Header *****/
  fscanf(farc,"%*s %i",&ncols[0]);
  fscanf(farc,"%*s %i",&nrows[0]);
  fscanf(farc,"%*s %lf",&ll_lng[0]);
  fscanf(farc,"%*s %lf",&ll_lat[0]);
  fscanf(farc,"%*s %lf",&cellsize[0]);
  fscanf(farc,"%*s %i",&NODATA[0]);

  /***** Allocate Latitude and Longitude Arrays for maximum size *****/
  Ncells    = ncols[0]*nrows[0];
  lat[0]    = (double *)calloc(Ncells,sizeof(double));
  lng[0]    = (double *)calloc(Ncells,sizeof(double));
  values[0] = (double *)   calloc(Ncells,sizeof(double));

  /***** Check for Valid Location *****/
  cell = 0;
  for(j=0;j<nrows[0];j++) {
    tmp_lat = ll_lat[0]+(double)(nrows[0]-j-0.5)*cellsize[0];
    for(i=0;i<ncols[0];i++) {
      tmp_lng = ll_lng[0] + (double)(i+0.5)*cellsize[0];
      fscanf(farc,"%lf",&tmpvalue);
      lat[0][cell]  = tmp_lat;
      lng[0][cell]  = tmp_lng;
      values[0][cell] = tmpvalue;
      cell++;
    }
  }
  fclose(farc);
  Ncells = cell;

  return Ncells;

}

int write_arcinfo_grid(char   *filename,
		       double *values,
		       double *lat,
		       double *lng,
		       int     ncols,
		       int     nrows,
		       double  ll_lat,
		       double  ll_lng,
		       double  cellsize,
		       int     NODATA,
		       int     Ncells) {
/**********************************************************************
  write_arcinfo_info       Keith Cherkauer           April 14, 1999

  This subroutine reads data from an ARC/INFO ascii grid and returns 
  the central latitude and longitude of all grid cells that are not 
  assigned the value of NODATA.  Routine returns the number of defined 
  grid cell.

**********************************************************************/

  FILE   *farc;
  int    i, j;
  int    cell;
  int    tmpvalue;

  if((farc=fopen(filename,"w"))==NULL) {
    fprintf(stderr,"ERROR: Unable of open ARCINFO grid %s for writing.\n",filename);
    exit(0);
  }

  /***** Write ARC/INFO Header *****/
  fprintf(farc,"ncols\t%i\n",ncols);
  fprintf(farc,"nrows\t%i\n",nrows);
  fprintf(farc,"xllcorner\t%lf\n",ll_lng);
  fprintf(farc,"yllcorner\t%lf\n",ll_lat);
  fprintf(farc,"cellsize\t%lf\n",cellsize);
  fprintf(farc,"NODATA_value\t%i\n",NODATA);

  /***** Check for Valid Location *****/
  cell = 0;
  for(j=0;j<nrows;j++) {
    for(i=0;i<ncols;i++) {
      if(values[cell] != NODATA) 
	fprintf(farc,"%lf",values[cell]);
      else
	fprintf(farc,"%i",NODATA);
      cell++;
      if(i<ncols-1) fprintf(farc,"\t");
      else fprintf(farc,"\n");
    }
  }
  fclose(farc);
  if(Ncells != cell) {
    fprintf(stderr,"ERROR: number of cells printed does not equal number of cells defined.\n");
    exit(0);
  }

  return cell;

}
