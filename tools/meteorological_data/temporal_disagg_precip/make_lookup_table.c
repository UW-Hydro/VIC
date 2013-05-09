/*
 * SUMMARY:      Creates a set of look-up tables (CDFs) for precipitation
 *               duration (D) and start-time/occurence time (S) on an hourly
 *               basis, by month.
 *               Takes as input the output from the preprocessing script 
 *               reformat_prcp.scr, which sets up both input files from
 *               raw Earthinfo CD hourly precipitation data.
 * AUTHOR:       Ed Maurer
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       edm@hydro.washington.edu
 * ORIG-DATE:    8/7/2000
 * DESCRIPTION:  
 * COMMENTS:     Uses same basic formulation at GMOD's program for the Ark-Red,
 *               in 1997. Many stats are calculated for other reasons but are
 *               not used in the program.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#define ICAT 5 /* number of bins into which P is divided */
#define MAXLINE 1024
#define intomm 25.4

void P_stats(float *prcp,int *rain_hrs, int *conseq_hrs,int *st_hr,
	     float *max_p,int *max_hr, int mo, int *pflag);

void main(int argc, char **argv)
{
  FILE *fpsta, *fpdat, *fpouts, *fpoutd;
  char sta_file[BUFSIZ+1], dat_file[BUFSIZ+1];
  char outs_file[BUFSIZ+1],outd_file[BUFSIZ+1];
  char str1[MAXLINE],state[3],stdat[3];
  int i,j,k,ii,st_num;
  int id,iddat;
  int yr, mo, dy;
  int seas_flag, fact;
  int ndays[12];
  int rain_hrs, conseq_hrs, max_hr;
  int pflag[24];
  int strt_hr;
  int sum_dur, sum_hrs;
  int p_dur[12][ICAT][24], p_hrs[12][ICAT][24];
  int temp_dur[4],temp_hrs[4];
  float prcp[24], max_p;
  float lat, lon, elev; 
  float sum_p;
  float acc_dur_old, acc_hrs_old;
  float accum_dur[24], accum_hrs[24];
  float pcats[ICAT+1] = {0, 5, 10, 15, 20, 9999};
  char season[] = { 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 0};
 
  if (argc!=6) {
    printf("usage: make_lookup_table <station info file> <hourly P data file> <output start time file> <output duration file> <season flag>\n");
    exit(EXIT_FAILURE); }

  strcpy(sta_file,argv[1]);        printf("station file: %s \n",sta_file);
  strcpy(dat_file,argv[2]);        printf("data file: %s \n",dat_file);
  strcpy(outs_file,argv[3]);       printf("output S file: %s \n",outs_file);
  strcpy(outd_file,argv[4]);       printf("output D file: %s \n",outd_file);
  seas_flag=atoi(argv[5]);         printf("season flag: %d \n",seas_flag);

  fact=1;
  if(seas_flag==1) {
    printf("\nseasonal values will be produced instead of monthly\n");
    fact=3;
  }

  /* open station info file */
  if((fpsta = fopen(sta_file,"r"))==NULL) {
    printf("Cannot open file %s \n",sta_file);exit(1);}
  /* open hourly data file */
  if((fpdat = fopen(dat_file,"r"))==NULL) {
    printf("Cannot open file %s \n",dat_file);exit(1);}
  /* open output S file -- lookup table */
  if((fpouts = fopen(outs_file,"a"))==NULL) {
    printf("Cannot open file %s \n",outs_file);exit(1);}
  /* open output D file -- lookup table */
  if((fpoutd = fopen(outd_file,"a"))==NULL) {
    printf("Cannot open file %s \n",outd_file);exit(1);}

  /*     initialize arrays */
  for(j=0;j<12;j++){
    ndays[j]=0;
  }

  /* Station data is sequential in data file. Start with first
     station in station information file, and collect all stats, print out
     month, CDF and location data in append mode into lookup table */

  fscanf(fpdat,"%s %d",stdat,&iddat);
  st_num=0;

  /***********get current station data from station info file ************/
  while(fgets(str1,MAXLINE,fpsta)!=NULL) {
    sscanf(str1,"%f %f %f %s %d",&lat, &lon, &elev, state, &id);
    st_num++;    
    fprintf(stderr,"Current Station %d lat/long: %.4f %.4f state/ID: %s %d\n",
	    st_num, lat, lon, state,id);
    /* initialize arrays for current station */
    for(k=0;k<24;k++){
      pflag[k]=0;
      for(i=0;i<12;i++){
	for(j=0;j<ICAT;j++){
	  p_hrs[i][j][k]=0;
	  p_dur[i][j][k]=0;
	}
      }
    }
    /**find all hourly data in data file corresponding to current station **/
    while(iddat==id){
      fscanf(fpdat,"%d %d %d",&yr,&mo,&dy);
      mo--; /* month arrays all start at 0 */
      sum_p=0;
      for(i=0;i<24;i++) {
	fscanf(fpdat,"%f",&prcp[i]);
	if(prcp[i]<0.0)sum_p=-999.0; /* if any hour missing, skip day */
	sum_p+=prcp[i];
      }

      /* if there is rain, calculate statistics for the current day */
      if(sum_p>0.0){
	P_stats(prcp, &rain_hrs, &conseq_hrs, &strt_hr, &max_p, &max_hr,
		mo, pflag);
	
	/* create a lookup table with statistics */

	/* find appropriate bin based on total daily P */
	sum_p*=intomm;
	for(i=0;i<ICAT;i++) {
	  if(sum_p>pcats[i] && sum_p<pcats[i+1]) {
	    p_dur[mo][i][rain_hrs]++;
	    /* rain occurences for each hour */
	    for(j=0;j<24;j++) p_hrs[mo][i][j]+=pflag[j];
	  }
	}
      }
      if(fscanf(fpdat,"%s %d",stdat,&iddat)!=2) {
	fprintf(stderr,"end of data file\n");
	iddat=0;
      }
    } /* end data retrieval for current station */

    /* normalize frequencies to probabilities -- total for current station */

    /* for seasonal output, aggregate monthly */
    if(seas_flag==1){ /* one CDF per season rather than month */
      for(j=0;j<ICAT;j++){
	for(k=0;k<24;k++){
	  for(i=0;i<4;i++){temp_dur[i]=0;temp_hrs[i]=0;}
	  for(i=0;i<12;i++){
	    ii=season[i];
	    temp_dur[ii]+=p_dur[i][j][k];
	    temp_hrs[ii]+=p_hrs[i][j][k];
	  }
	  for(i=0;i<4;i++){
	    p_dur[i][j][k]=temp_dur[i];
	    p_hrs[i][j][k]=temp_hrs[i];
	  }
	}
      }
    }
  
    for(j=0;j<ICAT;j++){
      for(i=0;i<(12/fact);i++){
	sum_dur=0;
	sum_hrs=0;
	for(k=0;k<24;k++){
	  sum_dur += p_dur[i][j][k];
	  sum_hrs += p_hrs[i][j][k];
	  }
	/* convert to cumulative probability */
	acc_dur_old=0;
	acc_hrs_old=0;
	for(k=0;k<24;k++){
	  if(sum_dur==0) accum_dur[k]=acc_dur_old;
	  else accum_dur[k]=acc_dur_old + (float) p_dur[i][j][k]/(float) sum_dur;
	  acc_dur_old=accum_dur[k];
	  if(sum_hrs==0) accum_hrs[k]=acc_hrs_old;
	  else accum_hrs[k]=acc_hrs_old + (float) p_hrs[i][j][k]/(float) sum_hrs;
	  acc_hrs_old=accum_hrs[k];
	}
	fprintf(fpouts,"%d %d %d ", st_num, j+1, i+1);
	fprintf(fpoutd,"%d %d %d ", st_num, j+1, i+1);
	for(k=0;k<24;k++){
	  fprintf(fpouts,"%.4f ",accum_hrs[k]);
	  fprintf(fpoutd,"%.4f ",accum_dur[k]);
	}
	fprintf(fpouts,"\n");
	fprintf(fpoutd,"\n");
      }
    }
  }
  fclose(fpsta);
  fclose(fpdat);
  fclose(fpouts);
  fclose(fpoutd);
}
/***************************************************************/
void P_stats(float *prcp,int *rain_hrs, int *conseq_hrs,int *st_hr,
	       float *max_p,int *max_hr, int mo, int *pflag)
{

  /* prcp          vector of 24 hourly prcp values
     rain_hrs      number of rainy hours for the current day
     conseq_hrs    number of consecutive hours with rain for current day
     st_hr         start hour of rainfall for current day
     max_p         maximum precip rate for current day
     max_hr        hour at which maximum precip rate occurred
     mo            month to which current day belongs
     pflag         vector of 1 and 0 indicating occurence of rain each hour
  */

  int i;
  int tmp_hrs=0;

  *rain_hrs=0;
  *max_hr=0;
  *max_p=0;
  *conseq_hrs=0;
  *st_hr=0;

  /* find max P, max hour, # of rain hours */
  for(i=0;i<24;i++){
    if(prcp[i]>0.0){
      *rain_hrs+=1;
      if(prcp[i]>(*max_p)){
	*max_p=prcp[i];
	*max_hr=i;
      }
    }
  }

  /* find duration and start hour */
  for(i=0;i<24;i++){
    pflag[i]=0;
    if(prcp[i]>0.0){
      tmp_hrs+=1;
      pflag[i]=1;
    }
    else if(tmp_hrs>*conseq_hrs){
      *conseq_hrs=tmp_hrs;
      *st_hr=i-*conseq_hrs;
      tmp_hrs=0;
    }
  }

  /*     catch for t=24 */
  if(tmp_hrs>*conseq_hrs){
    *conseq_hrs=tmp_hrs;
    *st_hr=i-*conseq_hrs;
  }
}

