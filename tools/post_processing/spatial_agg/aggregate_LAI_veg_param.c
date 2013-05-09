#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <math.h>

#define MAXSTR 1024
#define MAXVEGTYPES 14
#define MAXROOTZONES 3
#define DPLACES 4

typedef struct {
  int   cellnum;
  int   Nveg;
  int   cnt[MAXVEGTYPES];
  float fract[MAXVEGTYPES];
  float rootzone[MAXVEGTYPES][2][MAXROOTZONES];
  float LAI[MAXVEGTYPES][12];
} VEG_STRUCT;

int read_arcinfo_grid(char    *, float **, float **, int **, int *, int *,
		      float *, float *, float *, int *);

main(int argc, char *argv[]) {
/**********************************************************************
  aggragate_LAI_veg_param.c         Keith Cherkauer        July 6, 1999

  This program aggregates the resolution of a vegetation parameter file 
  by a factor of 2 (e.g. 1/8 to 1/4, 1/4 to 1/2).  It requires as input
  the high resolution vegetation parameter file, the high resolution
  basin ARCINFO grid cell number file, the name of the new low 
  resolution vegetation parameter file, the low resolution ARCINFO ASCII
  grid cell number file, and the number of root zones defined in the high
  resolution parameter file.  All parameters are aggregated to the 
  new resolution by finding their means.

  Modified from aggregate_veg_param.c on July 16, 1999 so that it can
  handle the new vegetation parameter files which contain monthly LAI
  values.

**********************************************************************/

  FILE  *fin, *fout;
  char   hr_maskname[MAXSTR];
  char   lr_maskname[MAXSTR];
  char   tmpstr[MAXSTR];
  float  lat, lng;
  float  tmpzone[2], tmpfract;
  float  tmpLAI;
  float *hr_lat, *hr_lng;
  float *lr_lat, *lr_lng;
  float  hr_ll_lat, hr_ll_lng, hr_cellsize;
  float  lr_ll_lat, lr_ll_lng, lr_cellsize;
  float  tmpsum;
  int    *hr_values, *lr_values;
  int    hr_ncols, hr_nrows;
  int    lr_ncols, lr_nrows;
  int    NODATA, storecnt;
  int    Ncells_lr, Ncells_hr;
  int    Nzones, tmpveg;
  int    zone, cell;
  int    Nfreq, Nfract;
  int    hr_cellnum;
  int    i, j, month;
  int    donecell;
  VEG_STRUCT *veg_param;

  if(argc!=6) {
    fprintf(stderr,"Usage: %s <high res vegparam file> <high res cellnum file> <low res vegparam file> <low res cellnum file> <# root zones>\n",argv[0]);
    fprintf(stderr,"\tThis program aggregates the given vegetation parameter file by a factor of 2 (e.g. 1/8 to 1/4, 1/4 to 1/2, etc.).\n");
    fprintf(stderr,"\t<high res vegparam file> is the vegetation parameter file for the high resolution model simualtion.\n");
    fprintf(stderr,"\t<high res cellnum file> is an ARCINFO ASCII grid cell number file for the high resolution model basin.\n");
    fprintf(stderr,"\t<low res vegparam file> is the vegetation parameter file for the low resolution (output) model simualtion.\n");
    fprintf(stderr,"\t<low res cellnum file> is an ARCINFO ASCII grid cell number file for the low resolution model basin.\n");
    exit(0);
  }

  /** Process command line arguments **/

  if((fin=fopen(argv[1],"r"))==NULL) {
    fprintf(stderr,"ERROR: Unable to open high resolution vegetation parameter file %s.\n",
	    argv[1]);
    exit(0);
  }

  strcpy(hr_maskname,argv[2]);

  if((fout=fopen(argv[3],"w"))==NULL) {
    fprintf(stderr,"ERROR: Unable to open new low resolution vegetation parameter file %s.\n",argv[3]);
    exit(0);
  }

  strcpy(lr_maskname,argv[4]);

  Nzones = atoi(argv[5]);

  /** Read mask files **/

  Ncells_lr = read_arcinfo_grid(lr_maskname,&lr_lat,&lr_lng,&lr_values,
				&lr_ncols,&lr_nrows,&lr_ll_lat,
				&lr_ll_lng,&lr_cellsize,&NODATA);

  Ncells_hr = read_arcinfo_grid(hr_maskname,&hr_lat,&hr_lng,&hr_values,
				&hr_ncols,&hr_nrows,&hr_ll_lat,
				&hr_ll_lng,&hr_cellsize,&NODATA);

  /** Allocate data arrays **/

  veg_param = (VEG_STRUCT *)calloc(Ncells_lr,sizeof(VEG_STRUCT));

  /** Process vegetation parameter file **/

  fscanf(fin,"%i %i", &hr_cellnum, &Nfract);
  lat = hr_lat[hr_cellnum-1];
  lng = hr_lng[hr_cellnum-1];

  fprintf(stdout,"Proccessed 000%%");
  donecell = 0;

  while(!feof(fin) && donecell<Ncells_hr) {

    for(cell=0;cell<Ncells_lr;cell++) {
      if((rint(lat * pow(10,(float)DPLACES)) 
	  == rint((lr_lat[cell]-hr_cellsize/2.) * pow(10,(float)DPLACES))
	  && rint(lng * pow(10,(float)DPLACES)) 
	  == rint((lr_lng[cell]-hr_cellsize/2.) * pow(10,(float)DPLACES)))
	 || (rint(lat * pow(10,(float)DPLACES)) 
	     == rint((lr_lat[cell]+hr_cellsize/2.) * pow(10,(float)DPLACES))
	     && rint(lng * pow(10,(float)DPLACES)) 
	     == rint((lr_lng[cell]-hr_cellsize/2.) * pow(10,(float)DPLACES)))
	 || (rint(lat * pow(10,(float)DPLACES)) 
	     == rint((lr_lat[cell]-hr_cellsize/2.) * pow(10,(float)DPLACES))
	     && rint(lng * pow(10,(float)DPLACES)) 
	     == rint((lr_lng[cell]+hr_cellsize/2.) * pow(10,(float)DPLACES)))
	 || (rint(lat * pow(10,(float)DPLACES)) 
	     == rint((lr_lat[cell]+hr_cellsize/2.) * pow(10,(float)DPLACES))
	     && rint(lng * pow(10,(float)DPLACES)) 
	     == rint((lr_lng[cell]+hr_cellsize/2.) * pow(10,(float)DPLACES)))
	 ) {

	for(i=0;i<Nfract;i++) {
	  fscanf(fin,"%i", &tmpveg);
	  veg_param[cell].cnt[tmpveg]++;
	  fscanf(fin,"%f", &tmpfract);
	  veg_param[cell].fract[tmpveg] += tmpfract;
      
	  for(zone=0;zone<Nzones;zone++) {
	    fscanf(fin,"%f %f",&tmpzone[0],&tmpzone[1]);
	    veg_param[cell].rootzone[tmpveg][0][zone] += tmpzone[0];
	    veg_param[cell].rootzone[tmpveg][1][zone] += tmpzone[1];
	  }

	  for(month=0;month<12;month++) {
	    fscanf(fin,"%f",&tmpLAI);
	    veg_param[cell].LAI[tmpveg][month] += tmpLAI;
	  }
	}
      }
    }
    fscanf(fin,"%i %i", &hr_cellnum, &Nfract);
    lat = hr_lat[hr_cellnum-1];
    lng = hr_lng[hr_cellnum-1];
    
    donecell++;
    fprintf(stdout,"\b\b\b\b%03i%%",donecell/Ncells_hr);
    fflush(stdout);
    
  }

  fprintf(stdout,"\b\b\b\b100%%\n");

  /** Verify that vegetation coverage for new cell sums to 
      that of old cell **/
  for(cell=0;cell<Ncells_lr;cell++) {
    tmpsum = 0;
    for(i=0;i<MAXVEGTYPES;i++) {
      if(veg_param[cell].cnt[i] > 0) {
	veg_param[cell].fract[i] /= (float)veg_param[cell].cnt[i];
	tmpsum += veg_param[cell].fract[i];
      }
    }
    if(tmpsum > 1) {
      for(i=0;i<MAXVEGTYPES;i++) {
	if(veg_param[cell].cnt[i] > 0) {
	  veg_param[cell].fract[i] /= tmpsum;
	}
      }
    }
  }

  /** Output new vegetation parameter files **/
  
  for(cell=0;cell<Ncells_lr;cell++) {

    Nfract = 0;
    for(i=0;i<MAXVEGTYPES;i++) {
      if(veg_param[cell].cnt[i] > 0) Nfract++;
    }
    fprintf(fout,"%i %i\n",cell+1,Nfract);
    for(i=0;i<MAXVEGTYPES;i++) {
      if(veg_param[cell].cnt[i] > 0) {
	fprintf(fout,"\t%i %f",i,
		veg_param[cell].fract[i]);
	for(zone=0;zone<Nzones;zone++) {
	  fprintf(fout," %f %f", veg_param[cell].rootzone[i][0][zone] 
		  / (float)veg_param[cell].cnt[i], 
		  veg_param[cell].rootzone[i][1][zone] 
		  / (float)veg_param[cell].cnt[i]);
	}
	fprintf(fout,"\n");
	for(month=0;month<12;month++) {
	  fprintf(fout," %.3f", veg_param[cell].LAI[i][month] 
		  / (float)veg_param[cell].cnt[i]);
	}
	fprintf(fout,"\n");
      }
    }
  }
}

int read_arcinfo_grid(char    *filename,
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


