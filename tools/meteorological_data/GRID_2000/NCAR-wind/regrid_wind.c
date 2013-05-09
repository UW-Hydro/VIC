/* regrid_wind program by EDM, May 1999, with the regridding routine */
/* taken from Laura's original reanalysis program */
/* Compile with the -lm flag to include the math library */
/* August 1999 revised interpolation to weight data propoerly */

/**************************************************************************/
/* This is the second step in processing the reanalysis wind data         */
/* for VIC input. The data is first retrieved with "getwind.scr".         */
/*                                                                        */
/* This program takes the yearly u-wind and v-wind NCEP/NCAR reanalysis   */
/* data and regrids it from the t62-gaussian grid to a grid corresponding */
/* to the input mask file. The required inputs are those listed in the    */
/* "run_regrid_wind.scr" file, and the files "gauss_t62_lat.list" and     */
/* "gauss_t62_lon.list" must reside in the same directory as this file.   */
/* These two files have the global lat/long coordinates of the reanalysis */
/* data. The output is one wind timeseries file for each grid cell active */
/* in the mask. This program also dumps two other files as output:        */
/*  1) wind_latlong.txt which is a list of file names for output          */
/*  2) cell_avg.out which is an average for each of 365 days for each cell*/
/* These two files are used in the fill_data.c program, used to fill the  */
/* wind data files with average daily wind speeds for each cell, for the  */
/* years of VIC data preceding the reanalysis data.                       */
/**************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define LEAPYR(y) (!((y)%400) || (!((y)%4) && ((y)%100)))
#define SCALE 0.01              /* scaling factor for reanalysis wind data */
#define OFFSET 225.45           /* offset of reanalysis wind data */
#define PI 3.1415
#define RADIUS 6371228.
double get_dist(double lat1, double long1, double lat2, double long2);
void make_flist(char msk[200]);            /* Function that makes a filelist */
void make_outfiles(char out_dir[200]);     /* function that make outputfiles */
int num_days(int year, int index_year, int index_mo); /* days to month in year */
const char flist[200]="wind_latlong.txt";    /* Created file with lat_long filenames */
const char ncar_lat[200]="./gauss_t62_lat.list"; /* latitudes of ncar data  */
const char ncar_lon[200]="./gauss_t62_lon.list"; /* longitudes of ncar data */

void main(int argc, char *argv[])
{

FILE *fu, *fv, *fposlat, *fposlon, *fm, *fpout,*fpflist,*fpmean;
char ufile[100], vfile[100],str1[200],out_dir[200],mask[200];
char old_out_dir[200], s1[200], maskfile[200],in_dir[200],old_in_dir[200];
char uftemp[100], vftemp[100];
int MASKROWS, MASKCOLS;
float MASKRES, void_nr,high,low, xll, yll, maskulong, maskulat;
int i,j,k,l,t;
long int count;
int lon, lat, min_lat, max_lat, min_lon, max_lon, skipdays;
float **uwnd, **vwnd,sum_wind, **maskx;
int ROWS,COLUMNS,SIZEM;
int twind, ii,cell_no;
int year, start_year, start_mo, end_year, end_mo, days_per_year,active_cells;
double ***wspeed,**cell_wind;
float *latitude, *longitude, uinterp, vinterp, dum;
float reflat, reflon,avg_wind;
double vertdist, reldist, hordist;
float newuleft, newvleft, newuright, newvright;
short int vic_wind, cellavg;

  if(argc!=12)
    {
      printf("\nUsage: %s <start_year> <start_mo> <end_year> <end_mo> <min_lat> <max_lat> <min_lon> <max_lon> <mask> <in_dir> <out_dir>\n",argv[0]);
      printf("See code for more details\n");
      exit(8);
    }

  printf("Assigning command line inputs to variables...\n\n");

   start_year = atoi(argv[1]); printf("start_year %d \n",start_year);
   start_mo = atoi(argv[2]);   printf("start_mo %d \n",start_mo);
   end_year = atoi(argv[3]);   printf("end_year %d \n",end_year);
   end_mo = atoi(argv[4]);     printf("end_mo %d \n",end_mo);
   min_lat = atoi(argv[5]);    printf("gauss_t62 coords:\nmin_lat %d \n",min_lat);
   max_lat = atoi(argv[6]);    printf("max_lat %d \n",max_lat);
   min_lon = atoi(argv[7]);    printf("min_lon %d \n",min_lon);
   max_lon = atoi(argv[8]);    printf("max_lon %d \n",max_lon);
   strcpy(mask,argv[9]);       printf("mask %s \n",mask);
   strcpy(in_dir,argv[10]); strcpy(old_in_dir,in_dir);
   strcpy(out_dir,argv[11]); strcpy(old_out_dir,out_dir); 
                printf("out_dir %s \n\n",old_out_dir);

  printf("Creating output file list %s\n",flist);
  strcpy(maskfile,mask);
  make_flist(mask);     /* Generates a filelist that will be used later  */
  printf("Making Output Files\n");
  make_outfiles(old_out_dir);/* Open all files in write mode and closes them again */
  ROWS=max_lat - min_lat +1; /* rows in reanalysis grid */
  COLUMNS=max_lon - min_lon + 1; /* columns in reanalysis grid */
  if((fpflist = fopen(flist,"r"))==NULL) {  
    printf("Cannot open file %s \n",flist);exit(0);} 

  /*------------------------------------------------------*/
  /* FIRST READ IN HEADER FOR MASK FILE (ARC/INFO STYLE)  */
  /*------------------------------------------------------*/
  if((fm = fopen(maskfile,"r")) == NULL) {
    printf("Cannot open/read mask\n"); exit(0); }
    fscanf(fm,"%*s %s",str1); MASKCOLS = atoi(str1);
    fscanf(fm,"%*s %s",str1); MASKROWS = atoi(str1);
    fscanf(fm,"%*s %s",str1); xll = atof(str1);
    fscanf(fm,"%*s %s",str1); yll = atof(str1);
    fscanf(fm,"%*s %s",str1); MASKRES = atof(str1);
    fscanf(fm,"%*s %s",str1); void_nr = atof(str1);
    high=void_nr+0.001;
    low=void_nr-0.001;

    /* convert to deg. east and make cell centered */
    maskulong = (xll+360.0)+MASKRES/2.0;
    /* upper left corner latitude cell center */
    maskulat = yll+(MASKROWS*MASKRES)-(MASKRES/2.0);
    /* initialize days_per_year for memory allocation */
    days_per_year=366;
    active_cells=0;
    SIZEM = MASKCOLS*MASKROWS;

    /*----------------------------*/
    /* ALLOCATE MEMORY TO ARRAYS  */
    /*----------------------------*/

  if(!(longitude = (float*) calloc(COLUMNS+1,sizeof(float)))) {
    printf("Cannot allocate memory to first record: BASIN\n");
    exit(8); }
  
  if(!(latitude = (float*) calloc(ROWS,sizeof(float)))) {
    printf("Cannot allocate memory to first record: BASIN\n");
    exit(8); }
  
  if(!(uwnd = (float**) calloc(ROWS,sizeof(float*)))) {
    printf("Cannot allocate memory to first record: BASIN\n");
    exit(8); }
  for(i=0; i<ROWS;i++) {
  if(!(uwnd[i]= (float*) calloc(COLUMNS,sizeof(float)))) {
    printf("Cannot allocate memory to first record: BASIN\n");
    exit(8); } }

  if(!(maskx = (float**) calloc(MASKROWS,sizeof(float*)))) {
    printf("Cannot allocate memory to first record: BASIN\n");
    exit(8); }
  for(i=0; i<MASKROWS;i++) {
  if(!(maskx[i] = (float*) calloc(MASKCOLS,sizeof(float)))) {
    printf("Cannot allocate memory to first record: BASIN\n");
    exit(8); } }

 if(!(vwnd = (float**) calloc(ROWS,sizeof(float*)))) {
    printf("Cannot allocate memory to first record: BASIN\n");
    exit(8); }
  for(i=0; i<ROWS;i++) {
  if(!(vwnd[i] = (float*) calloc(COLUMNS,sizeof(float)))) {
    printf("Cannot allocate memory to first record: BASIN\n");
    exit(8); } }
  
  if(!(wspeed = (double***) calloc(days_per_year,sizeof(double**)))) {
    printf("Cannot allocate memory to first record: WSPEED\n");
    exit(8); }
  for(j=0; j<days_per_year;j++) {
    if(!(wspeed[j] = (double**) calloc(MASKROWS,sizeof(double*)))) {
      printf("Cannot allocate memory to first record: BASIN\n");
      exit(8); }
    for(i=0; i<MASKROWS;i++) {
      if(!(wspeed[j][i] = (double*) calloc(MASKCOLS,sizeof(double)))) {
	printf("Cannot allocate memory to first record: BASIN\n");
	exit(8); } } }

  if(!(cell_wind = (double**) calloc(SIZEM,sizeof(double*)))) {
    printf("Cannot allocate memory to first record: CELL_WIND\n");
    exit(8); }
  for(j=0; j<SIZEM;j++) {
    if(!(cell_wind[j] = (double*) calloc(days_per_year,sizeof(double)))) {
      printf("Cannot allocate memory to first record: CELL_WIND\n");
      exit(8); } } 


  /*--------------------------------------------*/
  /* READ IN REMAINDER OF MASK FILE             */
  /*--------------------------------------------*/
    printf("reading mask file\n\n");
  for(i=0; i<MASKROWS;i++) {
    for(j=0; j<MASKCOLS; j++) {
      fscanf(fm,"%f",&maskx[i][j]);
      if(maskx[i][j]>high || maskx[i][j]<low ) active_cells++;
    }
  }
  fclose(fm);
  printf("%d active cells\n",active_cells);

  /*------------------------------------------------------*/
  /* BEGIN LOOP FOR EACH YEAR                             */
  /*------------------------------------------------------*/

  for(year=start_year;year<=end_year;year++) {
    printf("Processing year %d\n",year);
    count = 0; sum_wind = 0.0;

  /* NAMES OF REANALYSIS U-WIND AND V-WIND DATA FILES FOR CURRENT YEAR */

  sprintf(uftemp,"uwnd.%d.asc",year);
  sprintf(vftemp,"vwnd.%d.asc",year);
  strcat(in_dir,uftemp);
  strcpy(ufile,in_dir);
  strcpy(in_dir,old_in_dir);
  strcat(in_dir,vftemp);
  strcpy(vfile,in_dir);
  strcpy(in_dir,old_in_dir);
  printf("files %s and %s\n",ufile, vfile);
  days_per_year= num_days(year,end_year,end_mo);
  printf("days_per_year= %d\n",days_per_year);
  /* if (LEAPYR(year)) */
  /*  days_per_year=366; */
  /* else */
  /*  days_per_year=365; */

  /*--------------------------------------------*/
  /* OPEN REANALYSIS WIND FILES                 */
  /*--------------------------------------------*/

  if((fu = fopen(ufile,"r")) == NULL) {
    printf("Cannot open/read %s\n",ufile);
    exit(8); }

  if((fv = fopen(vfile,"r")) == NULL) {
    printf("Cannot open/read %s\n",vfile);
    exit(8); }

  if((fposlat = fopen(ncar_lat, "r")) == NULL) {
    printf("Cannot open/read %s\n",ncar_lat);
    exit(8); }
  if((fposlon = fopen(ncar_lon, "r")) == NULL) {
    printf("Cannot open/read %s\n",ncar_lon);
    exit(8); }
 
/*  printf("Files successfully opened...\n\n"); */

  /*----------------------------------------------------*/
  /* READ LAT/LON FILES -- POSITION IN REANALYSIS GRID  */
  /*----------------------------------------------------*/
  /* move to starting position for min lat and long (using reanalysis numbers 
             for coordinates in each file */
  for(i=1; i<min_lat;i++) fscanf(fposlat,"%f",&dum);
  for(i=1; i<min_lon;i++) fscanf(fposlon,"%f",&dum);

  for(i=0; i<ROWS;i++) 
    fscanf(fposlat,"%f",&latitude[i]);

  for(i=0; i<COLUMNS+1;i++) 
    fscanf(fposlon,"%f",&longitude[i]);
  
  fclose(fposlat);
  fclose(fposlon);


  /*-----------------------------------------------------------------*/
  /* READ IN REANALYSIS WIND GRID -- RUN FOR ALL DAYS IN ONE YEAR    */
  /*-----------------------------------------------------------------*/
  for(t=0; t<days_per_year;t++) {
      for(i=0; i<ROWS;i++) {
	for(j=0; j<COLUMNS;j++) {
	  fscanf(fu, "%d", &twind);
	  uwnd[i][j]= (float)twind*SCALE + OFFSET;
	  fscanf(fv, "%d", &twind);
	  vwnd[i][j]= (float)twind*SCALE + OFFSET;
	}
      }

      for(j=0; j<MASKROWS; j++) {
	reflat = maskulat - (float)j*MASKRES;
	for(i=0; i<MASKCOLS; i++) {
	  reflon = maskulong + (float)i*MASKRES;
	  lat=lon=9999;

      if(maskx[j][i]>high || maskx[j][i]<low) { /* if active cell,interpolate */

	  for(k=0; k<ROWS;k++) {
	    if(reflat<=latitude[k] && reflat>latitude[k+1]){
	      lat = k;
	    }
	  }
	  for(l=0; l<COLUMNS+1; l++) {
	    if(reflon>=longitude[l] && reflon<longitude[l+1]){
	      lon=l;
	    }
	  }
	  if(lat==9999 || lon==9999) {
	    printf("reflat/lon not in range\n");
	    exit(0);
	  }
	  
	  /*-----------------------------------*/
	  /* VERTICAL INTERPOLATION            */
	  /*-----------------------------------*/
	  vertdist=get_dist(latitude[lat], longitude[lon], latitude[lat+1],
			    longitude[lon]);
	  reldist=get_dist(latitude[lat],longitude[lon], reflat, longitude[lon]);

	  newuleft=((uwnd[lat][lon]*(1-(reldist/vertdist)))+(reldist/vertdist)
	    *uwnd[lat+1][lon]);
	  newvleft=((vwnd[lat][lon]*(1-(reldist/vertdist)))+(reldist/vertdist)
	    *vwnd[lat+1][lon]);

	  if(lon==COLUMNS) {
	     newuright=((uwnd[lat][0]*(1-(reldist/vertdist)))+
		     (reldist/vertdist)*uwnd[lat+1][0]);
	     newvright=((vwnd[lat][0]*(1-(reldist/vertdist)))+
		     (reldist/vertdist)*vwnd[lat+1][0]);
	  }
	  else {
	    newuright=((uwnd[lat][lon+1]*(1-(reldist/vertdist)))+
		       (reldist/vertdist)*uwnd[lat+1][lon+1]);
	    newvright=((vwnd[lat][lon+1]*(1-(reldist/vertdist)))+
		       (reldist/vertdist)*vwnd[lat+1][lon+1]);
	  }
	  /*-----------------------------------*/
	  /* HORIZONTAL INTERPOLATION          */
	  /*-----------------------------------*/

	  hordist=get_dist(reflat,longitude[lon],reflat,longitude[lon+1]);
	  reldist=get_dist(reflat,longitude[lon],reflat,reflon);

	  uinterp=newuleft*(1-(reldist/hordist)) + newuright*(reldist/hordist);
	  vinterp=newvleft*(1-(reldist/hordist)) + newvright*(reldist/hordist);
	  wspeed[t][j][i] = sqrt((double) (uinterp*uinterp + vinterp*vinterp));
      } /* end of if statement for interpolation only for active cell */
	} /* end loop for MASKCOLS */ 
      } /* end loop for MASKROWS */
  }    /* end loop for days per year */

  /*--------------------------------------------------*/
  /* OUTPUT TIMESERIES FOR EACH CELL IN MASK          */
  /*--------------------------------------------------*/
  printf("Writing output\n");
  skipdays=num_days( year, start_year, start_mo-1);
  if(year!=start_year) skipdays=0;
  printf("skipdays = %d\n",skipdays);
  cell_no=0;
  for(j=0; j<MASKROWS; j++) {
    for(i=0; i<MASKCOLS; i++) {
      if( maskx[j][i]>high || maskx[j][i]<low ) { /* if active cell, write */
	    fscanf(fpflist,"%s",s1);    
	    strcpy(out_dir,old_out_dir);
	    strcat(out_dir,s1);          /* create file name for grid cell */
	    if((fpout = fopen(out_dir,"ab"))==NULL) {
              printf("error opening output file %s\n",out_dir);exit(0);}
	    /* convert wind to short int values -- multiplying by 100 */
	    for(t=skipdays;t<days_per_year;t++) {
	      sum_wind += wspeed[t][j][i]; count+=1; /* calc sum for gross average */
	      vic_wind= (short int) (wspeed[t][j][i]*100);
	      fwrite(&vic_wind,sizeof(short int),1,fpout);
	    }
	    ii=skipdays;
	    for(t=skipdays;t<days_per_year;t++) {  /*avg daily wind for each cell */
      	      if( !(LEAPYR(year)) || t !=59){ /* skip feb 29 data for daily avg */
		cell_wind[cell_no][ii]+=wspeed[t][j][i]; /* for each active cell, sum for each day */
		ii++;}
	    }
	    cell_no++;
      fclose(fpout);
      }
    }
  }
  avg_wind = sum_wind/(float)count; /* gross average wind speed for year */
  printf("Average Daily Wind Speed Year %d = %.2f m/s\n\n",year,avg_wind);
  rewind(fpflist);

 fclose(fu); /* close reanalysis data files for the current year */
 fclose(fv);
  }  /* end year loop */
 /* write file with 365 daily averages for each cell */
  if((fpmean = fopen("cell_avg.out","wb")) == NULL) {
    printf("Cannot open cell_avg.out\n"); exit(0); }
  for(i=0;i<active_cells;i++){
    for(j=0;j<365;j++){
              cellavg = (short int) (cell_wind[i][j]/(end_year-start_year+1)*100);
	      fwrite(&cellavg,sizeof(short int),1,fpmean);
    }
  }
  fclose(fpmean);

} /* end main program */
/******************************************************************/
double get_dist(double lat1, double long1, double lat2, double long2)
{
  double theta1;
  double phi1;
  double theta2;
  double phi2;
  double dtor;
  double term1;
  double term2;
  double term3;
  double temp;
  double dist;

  dtor = 2.0*PI/360.0;
  theta1 = dtor*long1;
  phi1 = dtor*lat1;
  theta2 = dtor*long2;
  phi2 = dtor*lat2;
  term1 = cos(phi1)*cos(theta1)*cos(phi2)*cos(theta2);
  term2 = cos(phi1)*sin(theta1)*cos(phi2)*sin(theta2);
  term3 = sin(phi1)*sin(phi2);
  temp = term1+term2+term3;
  temp = (double) (1.0 < temp) ? 1.0 : temp;
  dist = RADIUS*acos(temp);

  return dist;
}
/**********************************************************************/
void make_flist(char msk[200])
{
/* This function reads a mask file and generates a text file         */
/* that contains the latitude and longitude to each gridcell         */
/* This will be the names of all the vicinput files                  */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

  FILE *fptrmask, *fptrlatlong;                /* Filepointers for files*/
  extern const char flist[200];
  char maskfile[200];
  char latlong[200];
  char str1[200];
  char void_nr[200];
  int cols, rows,i,j; 
  float uplftlat, uplftlong, llftlat,llftlong,resolution,longitude,value;
  float voidnr,high,low;

  strcpy(latlong,flist);
  strcpy(maskfile,msk);

  if((fptrmask = fopen(maskfile,"r"))==NULL){        /* Opens maskfile */
    printf("Cannot open file maskfile %s \n",maskfile);  exit(0);}
  if((fptrlatlong = fopen(latlong,"w"))==NULL){      /* Opens latlong file */
    printf("Cannot open file %s \n",latlong);  exit(0);}

  fscanf(fptrmask,"%*s %s",str1);        /* Reads number of cols */
  cols = atoi(str1);  /* Assigns number of cols to cols*/
  fscanf(fptrmask,"%*s %s",str1);        /* Reads number of rows */
  rows = atoi(str1);  /* Assigns number of rows to rows*/
  fscanf(fptrmask,"%*s %s", str1);        /* Reads llftlong gridposition */
  llftlong = atof(str1);  /* Assigns gridposition to llftlong */
  fscanf(fptrmask,"%*s %s", str1);        /* Reads llftlat gridposition */
  llftlat = atof(str1);  /* Assigns gridposition to llftlat*/  
  fscanf(fptrmask,"%*s %s",str1);        /* Reads resolution in degrees */
  resolution= atof(str1);     /* Assigns the resolution */
  fscanf(fptrmask,"%*s %s",void_nr);     /* Reads the void number */
  voidnr = atof(void_nr);
  high=voidnr+0.001;
  low=voidnr-0.001;

  /* place nodes at center of grid cells */
  uplftlong = llftlong + (resolution/2);
  uplftlat  = llftlat + ((float)(rows-1)*resolution) + (resolution/2);

  for(i=0;i<rows;i++)                 /* Write lat long to latlong.txt */
    {
      longitude = uplftlong;  
      for(j=0;j<cols;j++)
        {
         fscanf(fptrmask,"%s", str1);
         value = atof(str1);
	 if (value>high || value<low){
          fprintf(fptrlatlong,"wind_%6.4f_",uplftlat);
          fprintf(fptrlatlong,"%6.4f\n",longitude);}
        longitude = longitude + resolution;          /* keep track of position */
	}
    uplftlat = uplftlat - resolution;                /* keep track of position */
    }
  fclose(fptrmask);
  fclose(fptrlatlong);
} /* END of function void make_files(FILE *fp, char dir[400]) *********/
/***********************************************************************************/
void make_outfiles(char out_dir[200])
{
  FILE *fpflistx, *fpoutx;              /* filepointers to the files */
  char out[400];                 /* name of outputfile */
  extern const char flist[200];  /* name of filelist   "wind_latlong.txt"   */   
  char str1[200];                /* dummy string  */
  if((fpflistx = fopen(flist,"r"))==NULL){       /* opens filelist   */
    printf("Cannot open file %s \n",flist); exit(0);}
  while(fscanf(fpflistx,"%s", str1)!=EOF)    /* read the filelist to EOF */
    {
    strcpy(out,out_dir);        /* copy directory name */
    strcat(out,str1);           /* add filename   */
    fpoutx = fopen(out,"wb");    /* BINARY WRITE */
    fclose(fpoutx);       /* close the file. It will be reopened in append mode */
    }
  fclose(fpflistx);
}/** END function void make_outfiles(char out_dir[400] *********/
/***********************************************************************************/
int num_days(int year, int index_year, int index_mo)
  /* num of days to index month */
  /* index_mo is numbered 1-12 */
{
  int ndays,LPP;
  int days_nlp[] = {0,31,59,90,120,151,181,212,243,273,304,334,365};
  int days_lp[] = {0,31,60,91,121,152,182,213,244,274,305,335,366};

  if(LEAPYR(year)) LPP=1;
  else LPP=0;

  if(year==index_year) {
    if(LEAPYR(year)) ndays=days_lp[index_mo];
    else ndays=days_nlp[index_mo];
  }
  else ndays = 365+LPP;

  return ndays;
}
/***********************************************************************************/
