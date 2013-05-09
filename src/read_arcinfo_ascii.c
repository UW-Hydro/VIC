#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
#include <math.h>

static char vcid[] = "$Id$";

double read_arcinfo_value(char *filename,
			  double lat,
			  double lng) {
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
  double ll_lat;
  double ll_lng;
  double cellsize;
  double NODATA;
  double value;
  double tmp_lat;
  double tmp_lng;
  double tmpvalue;
  char   errstr[MAXSTRING];

  if((farc=fopen(filename,"r"))==NULL) {
    sprintf(errstr,"Unable to open ARC/INFO soil file %s",filename);
    nrerror(errstr);
  }

  /***** Read ARC/INFO Header *****/
  fscanf(farc,"%*s %d",&ncols);
  fscanf(farc,"%*s %d",&nrows);
  fscanf(farc,"%*s %lf",&ll_lng);
  fscanf(farc,"%*s %lf",&ll_lat);
  fscanf(farc,"%*s %lf",&cellsize);
  fscanf(farc,"%*s %lf",&NODATA);

  /***** Check for Valid Location *****/
  if(lat<ll_lat || lat>ll_lat+cellsize*(double)nrows) {
    sprintf(errstr,"Given latitude %f does not fall within ARC file %s",
	    lat,filename);
    nrerror(errstr);
  }

  if(lng<ll_lng || lng>ll_lng+cellsize*(double)ncols) {
    sprintf(errstr,"Given longitude %f does not fall within ARC file %s",
	    lng,filename);
    nrerror(errstr);
  }

  value = NODATA;
  for(j=0;j<nrows;j++) {
    tmp_lat = ll_lat+(double)(nrows-j-0.5)*cellsize;
    for(i=0;i<ncols;i++) {
      tmp_lng = ll_lng + (double)(i+0.5)*cellsize;
      fscanf(farc,"%lf",&tmpvalue);
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
		      double **lat,
		      double **lng,
		      int    **cellnum) {
/**********************************************************************
  read_arcinfo_info           Keith Cherkauer           May 5, 1998

  This subroutine reads data from an ARC/INFO ascii grid and returns 
  the central latitude and longitude of all grid cells that are not 
  assigned the value of NODATA.  Routine returns the number of defined 
  grid cell.

**********************************************************************/

  FILE   *farc;
  int    i, j;
  int    ncols;
  int    nrows;
  double ll_lat;
  double ll_lng;
  double cellsize;
  double tmp_lat;
  double tmp_lng;
  int    cell;
  int    Ncells;
  int    NODATA;
  int    tmpvalue;

  farc=open_file(filename,"r");

  /***** Read ARC/INFO Header *****/
  fscanf(farc,"%*s %d",&ncols);
  fscanf(farc,"%*s %d",&nrows);
  fscanf(farc,"%*s %lf",&ll_lng);
  fscanf(farc,"%*s %lf",&ll_lat);
  fscanf(farc,"%*s %lf",&cellsize);
  fscanf(farc,"%*s %d",&NODATA);

  /***** Allocate Latitude and Longitude Arrays for maximum size *****/
  Ncells  = ncols*nrows;
  lat[0]  = (double *)calloc(Ncells,sizeof(double));
  lng[0]  = (double *)calloc(Ncells,sizeof(double));
  cellnum[0] = (int *)   calloc(Ncells,sizeof(int));

  /***** Check for Valid Location *****/
  cell = 0;
  for(j=0;j<nrows;j++) {
    tmp_lat = ll_lat+(double)(nrows-j-0.5)*cellsize;
    for(i=0;i<ncols;i++) {
      tmp_lng = ll_lng + (double)(i+0.5)*cellsize;
      fscanf(farc, "%d", &tmpvalue);
      if(tmpvalue != NODATA) {
	lat[0][cell]     = tmp_lat;
	lng[0][cell]     = tmp_lng;
	cellnum[0][cell] = tmpvalue;
	cell++;
      }
    }
  }
  fclose(farc);
  Ncells = cell;

  return Ncells;

}
