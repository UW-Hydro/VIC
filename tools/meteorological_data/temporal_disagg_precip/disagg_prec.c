/*
 * SUMMARY:      Disaggregates daily VIC precipitation from a 4-column binary 
 *               input data (forcing) file. For each file in the input file
 *               list, it is assumed the first column (2-byte, unsigned short
 *               int) is daily precipitation total, which is disaggregated to 
 *               hourly or other sub-daily increment (2,3,4,6,8, or 12 hours)
 * AUTHOR:       Ed Maurer
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       edm@hydro.washington.edu
 * ORIG-DATE:    8/9/2000
 * DESCRIPTION:
 * COMMENTS:     Uses same basic formulation as GMOD's program for the Ark-Red,
 *               in 1997. But of course is superior.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#define LEAPYR(y) (!((y)%400) || (!((y)%4) && ((y)%100)))
#define ICAT 5
#define MAXLINE 256
#define TRUE 1
#define FALSE 0
#define ASCII_INPUT FALSE  /* FALSE indicates 2-byte binary input  */
#define ASCII_OUTPUT FALSE /* FALSE indicates 2-byte binary output */

float *get_1d_mem_float(int dim1);
unsigned short int *get_1d_mem_usi(int dim1);
float **get_2d_mem_float(int dim1, int dim2);
float ****get_4d_mem_float(int dim1, int dim2, int dim3, int dim4);
float ran1(long *idum);
int find_near(float,float,float *,float *,int);
void get_latlong( char fname[MAXLINE], float *lat, float *lon );
float getdist(float lat1, float long1, float lat2, float long2);

void main(int argc, char **argv)
{
  FILE *fpsta, *fps, *fpd, *fpflist, *fpgrid, *fpout;
 
  float p_tmp[24];
  int days_in_mo[] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };
  int st_num=0,n_sta=0,n_grids=0,n_days=0;
  int d_hr;
  int h,i,j,k,jj;
  int st_yr,st_mo;
  int yr,mo,dy;
  int temp_days, mo_sea;
  int p_mult;
  int seas_flag, seas_fact=1;
  int sta,bin,mon,sta2,bin2,mon2;
  int ibin,idur,istrt,idur_tmp,istrt_tmp,iend_tmp,loop;
  long int idum=123; /* seed for ran1 random number generator */
  unsigned short int usprec;
  short int dummy_int[4];
  float *sta_lat, *sta_lon;
  float p_bin[ICAT] = { 5, 10, 15, 20, 9999 };
  float psum_tmp;
  float *prec;
  float ****dur, ****strt;
  float grid_lat,grid_lon;
  float seed;
  char fname[MAXLINE];
  char sta_file[BUFSIZ+1],flist[BUFSIZ+1];
  char s_file[BUFSIZ+1],d_file[BUFSIZ+1];
  char str1[MAXLINE], temp_fname[MAXLINE];
  char out_dir[BUFSIZ+1], outname[BUFSIZ+1], out[BUFSIZ+1];
  char season[] = { 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 0};
 
  if (argc!=11) {
    fprintf(stderr,"usage: disagg_prec <station info file> <start/hour lookup file> <duration lookup file> <out_dir> <file list> <disaggregation increment> <season flag>\n");
    exit(EXIT_FAILURE); }

  strcpy(sta_file,argv[1]); fprintf(stderr,"station file: %s \n",sta_file);
  strcpy(s_file,argv[2]);   fprintf(stderr,"start/hour lookup file: %s \n",s_file);
  strcpy(d_file,argv[3]);   fprintf(stderr,"duration lookup file: %s \n",d_file);
  strcpy(out_dir,argv[4]);  fprintf(stderr,"out_dir: %s \n",out_dir);
  strcpy(flist,argv[5]);
  d_hr=atoi(argv[6]);         
  seas_flag=atoi(argv[7]);  fprintf(stderr,"season flag: %d \n",seas_flag);
  p_mult=atoi(argv[8]);     fprintf(stderr,"VIC P multiplier: %d \n",p_mult);
  st_yr=atoi(argv[9]);      fprintf(stderr,"Start year of VIC data: %d \n",st_yr);
  st_mo=atoi(argv[10]);     fprintf(stderr,"Start month of VIC data: %d \n",st_mo);

  if((d_hr>=1 && d_hr<=4) || d_hr==6 || d_hr==8 || d_hr==12)
    fprintf(stderr,"daily P to be disaggregated to %d hours\n",d_hr);
  else {printf("Invalid disaggregation hours. Must be 1,2,3,4,6,8 or 12\n");
  exit(1);}

  if(seas_flag==1) {
    fprintf(stderr,"lookup tables contain seasonal values\n");
    seas_fact=3;
  }

  /* open station info file */
  if((fpsta = fopen(sta_file,"r"))==NULL) {
    printf("Cannot open file %s \n",sta_file);exit(1);}
  /* open file list */
  if((fpflist = fopen(flist,"r"))==NULL) {
    printf("Cannot open file %s \n",flist);exit(1);}

  /* count number of hourly data stations, VIC grid cells, days of data */
  while(fgets(str1,MAXLINE,fpsta)!=NULL)n_sta++;
  fprintf(stderr,"total number of hourly met stations: %d\n",n_sta);
  rewind (fpsta);

  while(fgets(str1,MAXLINE,fpflist)!=NULL) n_grids++;
  fprintf(stderr,"total number of VIC grid cells: %d\n",n_grids);
  rewind (fpflist);
  sscanf(str1,"%s",temp_fname);

  if(ASCII_INPUT){
    if((fpgrid = fopen(temp_fname,"r"))==NULL) {
      printf("Cannot open file %s \n",temp_fname);exit(1);}
    while (fgets(str1,MAXLINE,fpgrid) != '\0') n_days++;
  }
  else {
    if((fpgrid = fopen(temp_fname,"rb"))==NULL) {
      printf("Cannot open file %s \n",temp_fname);exit(1);}
    while((fread(&dummy_int,sizeof(short int),4,fpgrid) != 0)) n_days++;
  }
  fprintf(stderr,"total number of days in VIC data files: %d\n",n_days);
  fclose(fpgrid);
  
  /* allocate memory */
  sta_lat=get_1d_mem_float(n_sta);
  sta_lon=get_1d_mem_float(n_sta);
  prec=get_1d_mem_float(n_days);
  dur=get_4d_mem_float(n_sta,12/seas_fact+1,ICAT,24);
  strt=get_4d_mem_float(n_sta,12/seas_fact+1,ICAT,24);

  /* read station location data into array */
  for(i=0;i<n_sta;i++) {
    fgets(str1,MAXLINE,fpsta);
    sscanf(str1,"%f %f",&sta_lat[i], &sta_lon[i]);
  }
  fclose(fpsta);

  /* open S (hourly occurence) file -- lookup table */
  if((fps = fopen(s_file,"r"))==NULL) {
    printf("Cannot open file %s \n",s_file);exit(1);}
  /* open D (duration CDF) file -- lookup table */
  if((fpd = fopen(d_file,"r"))==NULL) {
    printf("Cannot open file %s \n",d_file);exit(1);}
  /* read in all CDF data into arrays */
  fprintf(stderr,"loading data arrays into memory\n");
  for(i=0;i<n_sta;i++){
    for(j=0;j<ICAT;j++){
      for(k=0;k<12/seas_fact;k++){
	fscanf(fps,"%d %d %d",&sta,&bin,&mon);
	fscanf(fpd,"%d %d %d",&sta2,&bin2,&mon2);
	if(sta!=sta2 || bin!=bin2 || mon!=mon2){
	  fprintf(stderr,"lookup tables do not match\n");
	  fprintf(stderr,"sta=%d sta2=%d bin=%d bin2=%d mon=%d mon2=%d\n",
		  sta,sta2,bin,bin2,mon,mon2);
	  exit(1);
	}
	for(h=0;h<24;h++){
	  /* in lookup table, arrays start at 1, adjust here */
	  fscanf(fps,"%f",&strt[sta-1][mon-1][bin-1][h]);
	  fscanf(fpd,"%f",&dur[sta-1][mon-1][bin-1][h]);
	}
      }
    }
  }
  fclose(fps);
  fclose(fpd);

  /*****************************************************/
  /* Loop to disaggregate each grid cell precipitation */
  /*****************************************************/

  for(i=0;i<n_grids;i++){
    fprintf(stderr,"Disaggregating grid cell %d of %d ",i+1,n_grids);

    /* read in name of current grid cell */
    fscanf(fpflist,"%s",fname);
    fprintf(stderr,"%s\n",fname);

    /* read in precipitation for current grid cell */
    if(ASCII_INPUT){
      if((fpgrid = fopen(fname,"r"))==NULL) {
	printf("Cannot open file %s \n",fname);exit(1);}
      for(j=0;j<n_days;j++) {
	fscanf(fpgrid,"%f %*s %*s %*s",&prec[j]);
      }
    }
    else {
      if((fpgrid = fopen(fname,"rb"))==NULL) {
	printf("Cannot open file %s \n",fname);exit(1);}
      for(j=0;j<n_days;j++) {
	if(fread(&usprec,sizeof(unsigned short int),1,fpgrid) != 1)
	  printf("read prec[j] err rec %d %s\n", j,fname);
	if(fread(&dummy_int,sizeof(short int),3,fpgrid) != 3)
	  printf("read dummy_int err rec %d %s\n", j,fname);
	prec[j]=(float) usprec/p_mult;
      }
    }
    fclose(fpgrid);

    /* get lat and long for current grid cell */
    get_latlong(fname, &grid_lat, &grid_lon);
    
    /* find nearest neighbor for current grid cell */
    st_num=find_near(grid_lat, grid_lon, sta_lat, sta_lon, n_sta);

    /* open output file for the current grid cell */
    sprintf(outname,"data_%.4f_%.4f",grid_lat,grid_lon);
    strcpy(out,out_dir);
    strcat(out,outname);
    if(ASCII_OUTPUT){
    if((fpout = fopen(out,"w"))==NULL) {
      printf("Cannot open file %s \n",out);exit(1);}
    }
    else {
    if((fpout = fopen(out,"wb"))==NULL) {
      printf("Cannot open file %s \n",out);exit(1);}
    }

    dy=0;
    mo=st_mo-1; /* adjust so array starts at 0 */
    yr=st_yr;

    while(dy<n_days) {
      temp_days=days_in_mo[mo];
      if(LEAPYR(yr) && mo==1) temp_days++;

      /*** loop for all days in current month ***/
      for(k=0;k<temp_days;k++) {

	/* disaggregate current day for current grid cell */
	/* if appropriate, find season for current day */
	if(seas_flag==1) mo_sea = season[mo];
	else mo_sea=mo;

	/* if no precip, write zero & move to next record */
	if(prec[dy]<0.0001){
	  for(j=0;j<24/d_hr;j++) {
	    if(ASCII_OUTPUT) fprintf(fpout,"%.2f\n",0.0);
	    else {
	      usprec=(unsigned short int) (0*p_mult);
	      fwrite(&usprec,sizeof(unsigned short int),1,fpout);
	    }
	  }
	}
	else {
	  /* find prec bin */
	  for(j=0;j<ICAT;j++){
	    if(prec[dy]<p_bin[j]){
	      ibin=j;
	      break;
	    }
	  }

	  /* find duration */
	  seed=ran1(&idum);
	  for(j=0;j<24;j++){
	    if(seed<dur[st_num][mo_sea][ibin][j]){
	      idur=j+1; /* duration must be between 1 and 24 */
	      break;
	    }
	  }

	  /* find start/hour of center of event */
	  seed=ran1(&idum);
	  for(j=0;j<24;j++){
	    if(seed<strt[st_num][mo_sea][ibin][j]){
	      istrt=j;
	      break;
	    }
	  }

	  /*  fprintf(stderr,"%d %d %d %d %d %f\n",yr, mo+1, dy+1, istrt+1, idur,
	      seed); */

	  /* initialize p_tmp  -- output at desired disaggregation level */
	  for(j=0;j<24;j++) p_tmp[j]=0;

 	  /* center daily precip about the selected hour */
	  seed=ran1(&idum);
	  /* if even number of hours, decide whether extra hour falls
	     at beginning or end of interval */
	  if(!(idur%2)) {
	    idur_tmp=idur-1;
	    istrt_tmp=istrt-idur_tmp/2;
	    iend_tmp=istrt+idur_tmp/2;
	    if(seed>0.5) iend_tmp=iend_tmp+1;
	    else istrt_tmp=istrt_tmp-1;
	  }
	  else {
	    istrt_tmp=istrt-idur/2;
	    iend_tmp=istrt+idur/2;
	  }
	  /* move time interval backward if it extends past 24 hrs */
	  if(iend_tmp>23) {
	    istrt_tmp=istrt_tmp-(iend_tmp-23);
	    iend_tmp=23;
	  }
	  else if(istrt_tmp<0){
	    iend_tmp=iend_tmp-istrt_tmp;
	    istrt_tmp=0;
	  }

	  /* disaggregate to hourly */	
	  loop=0;
	  for(j=istrt_tmp;j<=iend_tmp;j++){
	    loop=loop+1;
	    p_tmp[j]=prec[dy]/idur;
	  }
	  if(loop!=idur){
	    fprintf(stderr,"Error in disaggregation\n");
	    fprintf(stderr,"record %d loop=%d idur=%d\n",dy,loop,idur);
	    fprintf(stderr,"start=%d end=%d\n",istrt_tmp,iend_tmp);
	    exit(1);
	  }

	  /* aggregate back up to desired time resolution and write */
	  for(j=0;j<24/d_hr;j++){
	    psum_tmp=0;
	    for(jj=0;jj<d_hr;jj++) psum_tmp+=p_tmp[j*d_hr+jj];
	    if(ASCII_OUTPUT) fprintf(fpout,"%.2f\n",psum_tmp);
	    else {
	      usprec=(unsigned short int) (psum_tmp*p_mult);
	      fwrite(&usprec,sizeof(unsigned short int),1,fpout);
	    }
	  }

	} /* end of "else" loop before finding prec bin */
	dy++;
      } /* end of loop for k, all days in current month */
      mo++;
      if(mo>11) {mo=0;yr++;}
    } /*end of while loop for all days in VIC grid file */
    fclose(fpout);
  } /* end of loop for i, all grid cells in file list */
  fclose(fpflist);
}
/*************************FUNCTIONS***************/
/* random number generator */
float ran1(long *idum)
{
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

  int j;
  long k;
  float temp;
  static long iy=0;
  static long iv[NTAB];

  if (*idum <= 0 || !iy ) {
    if (-(*idum) < 1) *idum=1;
    else *idum = -(*idum);
    for(j=NTAB+7;j>=0;j--){
      k=(*idum)/IQ;
      *idum=IA*(*idum-k*IQ)-IR*k;
      if (*idum<0) *idum+=IM;
      if (j < NTAB) iv[j]= *idum;
    }
    iy=iv[0];
  }
  k=(*idum)/IQ;
  *idum=IA*(*idum-k*IQ)-IR*k;
  if(*idum<0) *idum+=IM;
  j=iy/NDIV;
  iy=iv[j];
  iv[j]=*idum;
  if ((temp=AM*iy) > RNMX) return RNMX;
  else return temp;
}
/*************************************************************/
/* finds nearest station for current grid cell/VIC file */
/* returns station number associated with nearest hourly station */
int find_near(float grid_lat,float grid_lon,float *sta_lat,float *sta_lon, int n_sta)
{
  float min_dist;
  int i,iflag;
  float *dist;

  dist=get_1d_mem_float(n_sta);
  /* find distances to each station, minimum */
  min_dist=999999;
  for(i=0;i<n_sta;i++){
    dist[i]=getdist(grid_lat,grid_lon,sta_lat[i],sta_lon[i]);
    if(dist[i]<min_dist) {
      min_dist=dist[i];
      iflag=i;
    }
  }
  free(dist);
  return iflag;
}
/************************************************************************/
float getdist(float lat1, float long1, float lat2, float long2)
{
  /* find distance between pair of lat long coordinates 
	 south and west are negative */
#define pi 3.141593
#define radius 6378.0
#define dtor 2.0*pi/360.0
      
  float theta1, theta2, phi1, phi2, temp;
  float term1, term2, term3;
      
  theta1 = dtor*long1;
  phi1 = dtor*lat1;
  theta2 = dtor*long2;
  phi2 = dtor*lat2;
  term1 = cos(phi1)*cos(theta1)*cos(phi2)*cos(theta2);
  term2 = cos(phi1)*sin(theta1)*cos(phi2)*sin(theta2);
  term3 = sin(phi1)*sin(phi2);
  temp = term1+term2+term3;
  temp = temp>1.0 ? 1.0:temp;
  return radius*acos(temp);
}
/***********************************************************/
void get_latlong( char fname[MAXLINE], float *templat, float *templong )
  /* extracts the lat and long from the current cell from its */
  /* filename, e.g. data_35.6875_-110.4375 */
{
  int i,j, delim;
  char latlong[2][BUFSIZ+1];

  for(i=0;i<2;i++){
    delim=strlen(fname)-1;
    while(fname[delim]!='_')delim--;
    fname[delim]='\0';
    delim++;
    j=0;
    while(fname[delim]!='\0'){
      latlong[i][j]=fname[delim];
      fname[delim]='\0';
      delim++;
      j++;
    }
    latlong[i][j]='\0';
  }
  *templong = atof(latlong[0]);
  *templat = atof(latlong[1]);
}
/***********************************************************/
float *get_1d_mem_float(int dim1)
{
  // allocate memory for a 1-d array, pass back pointer
  float *array_1d;

  if (!(array_1d = (float *) calloc(dim1, sizeof(float)))) {
    fprintf(stderr, "unable to allocate memory in get_1d_mem_float\n");
    fprintf(stderr, "dimension %d\n", dim1);
    exit(1);
  }
  return array_1d;
}
/***********************************************************/
unsigned short int *get_1d_mem_usi(int dim1)
{
  // allocate memory for a 1-d array, pass back pointer
  unsigned short int *array_1d;

  if (!(array_1d = (unsigned short int *) calloc(dim1, sizeof(unsigned short int)))) {
    fprintf(stderr, "unable to allocate memory in get_1d_mem_usi\n");
    fprintf(stderr, "dimension %d\n", dim1);
    exit(1);
  }
  return array_1d;
}
/***********************************************************/
float **get_2d_mem_float(int dim1, int dim2)
{
  // allocate memory for a 2-d array, pass back pointer
  float **array_2d;
  int i;

  if (!(array_2d = (float **) calloc(dim1, sizeof(float *)))) {
    fprintf(stderr, "unable to allocate memory in get_2d_mem_float\n");
    fprintf(stderr, "dimension1= %d\n", dim1);
    exit(1);
  }
  for (i = 0; i < dim1; i++) {
    if (!(array_2d[i] = (float *) calloc(dim2, sizeof(float)))) {
      fprintf(stderr, "unable to allocate memory in get_2d_mem_float\n");
      fprintf(stderr, "dimension1= %d, dimension2=%d\n", dim1,dim2);
      exit(1);
    }
  }
  return array_2d;
}
/***********************************************************/
float ****get_4d_mem_float(int dim1, int dim2, int dim3, int dim4)
{
  float ****array_4d;
  int i,j,k;

  if(!(array_4d = (float****) calloc(dim1,sizeof(float***)))) {
    printf("Cannot allocate memory in get_4d_mem_float\n");
    exit(1); }
  for(k=0; k<dim1;k++) {
    if(!(array_4d[k] = (float***) calloc(dim2,sizeof(float**)))) {
      printf("Cannot allocate memory in get_4d_mem_float\n");
      exit(1); }
    for(j=0; j<dim2;j++) {
      if(!(array_4d[k][j] = (float**) calloc(dim3,sizeof(float*)))) {
	printf("Cannot allocate memory in get_4d_mem_float\n");
	exit(1); }
      for(i=0; i<dim3;i++) {
	if(!(array_4d[k][j][i] = (float*) calloc(dim4,sizeof(float)))) {
	  printf("Cannot allocate memory in get_4d_mem_float\n");
	  exit(1); }
      } 
    }
  }
  return array_4d;
}
