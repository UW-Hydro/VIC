#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <math.h>

#define NUM_DATE_VALS 4
#define NUM_DATA_VALS 15
#define NUM_LAYERS    3
#define NUM_FRONTS    3  /* maximum number of freezing fronts if FROZEN 
                            SOIL activated */

#define FALSE 0
#define TRUE 1

double juldate(int, int, int, double);
int get_file_time_step(FILE **, int *);
void get_record(FILE **, int **, float **,int,int);

main(int argc, char *argv[])
/**********************************************************************
  extract_daily_means_from_LDAS_binary.c  
                          Keith Cherkauer          September 30, 1999

  This program reads the LDAS binary files and extracts daily average
  or summed values for all model output variables for the specified 
  period of interest.  
  The output file is in ASCII column format: 
  <year> <month> <day> <prec> <evap> .....

  LDAS binary format:
  unsigned short int         year
  char                       month
  char                       day
  (char)                     (hour)
  unsigned short int         prec * 100              sum
  signed short int           evap * 100              sum
  float                      runoff                  sum
  float                      baseflow                sum
  unsigned short int         moist[Nlayers] * 10     mean
  unsigned short int         swq * 100               mean
  signed short int           net_short * 10          mean
  signed short int           in_long * 10            mean
  signed short int           r_net * 10              mean
  signed short int           latent * 10             mean
  signed short int           sensible * 10           mean
  signed short int           grnd_flux * 10          mean
  unsigned short int         albedo * 10000          mean
  signed short int           surf_temp * 100         mean
  unsigned short int         rel_humid * 100         mean
  signed short int           air_temp * 100          mean
  if FROZEN_SOIL activated
    unsigned short int         ice[Nlayers] * 10     mean
    for each NUM_FRONTS:
      unsigned short int         fdepth[] * 100     mean
      unsigned short int         tdepth[] * 100     mean

  Modified:
  11-02-99 modified to determine if file contains frozen soil
           information along with model time step.     KAC 

**********************************************************************/
{
  FILE  *infile, *outfile, *files;
  int    *date;
  int    *tmpdate, *lastdate;
  int    i, l, Nfiles;
  int    dt;
  int    Ncols;
  float  *data, *sum;
  float  lat, lon;
  char   name[600], out_name[600], file_index[35], resultdir[512], outdir[512];
  char   model[20];
  int    FIRSTFILE = TRUE;
  int    DAILY     = FALSE;
  int    FROZEN;
  int    DONE;
  int    datacnt;
  int    *startdate, *enddate;
  double startjday, endjday, tmpjday, lastjday;

  date      = (int *)malloc(NUM_DATE_VALS*sizeof(int));
  tmpdate   = (int *)malloc(NUM_DATE_VALS*sizeof(int));
  lastdate  = (int *)malloc(NUM_DATE_VALS*sizeof(int));
  enddate   = (int *)malloc(NUM_DATE_VALS*sizeof(int));
  startdate = (int *)malloc(NUM_DATE_VALS*sizeof(int));

  if(argc!=8) {
    fprintf(stderr,"Usage: %s <result dir> <result file prefix> ",argv[0]);
    fprintf(stderr,"<output directory> <file extension list> ");
    fprintf(stderr,"<number of files> <start date> <end date>\n");
    fprintf(stderr,"\tThis program reads the LDAS binary files and extracts daily average or summed parameters for the specified period of interest.  The output file is in ASCII column format: \n");
    fprintf(stderr,"\t\t<year> <month> <day> <prec> <evap> ...\n");
    fprintf(stderr,"\n\tAll dates are in the form MMDDYYYY (MM = month, DD = day, YYYY = year) and computations start at midnight the morning of that day.\n");
    fprintf(stderr,"\n\tThis program is currently configured to read files with %i soil layers.\n",NUM_LAYERS);
    exit(0);
  }

  strcpy(resultdir,argv[1]);
  strcpy(model,argv[2]);
  strcpy(outdir,argv[3]);
  files        = fopen(argv[4], "r");
  Nfiles       = atoi(argv[5]);
  startdate[1] = (int)(atof(argv[6])/1000000);
  startdate[2] = (int)(atof(argv[6])/10000) - startdate[1]*100;
  startdate[0] = (int)(atof(argv[6])) - startdate[1]*1000000 
    - startdate[2]*10000;
  startjday    = juldate(startdate[0],startdate[1],startdate[2],0);
  enddate[1]   = (int)(atof(argv[7])/1000000);
  enddate[2]   = (int)(atof(argv[7])/10000) - enddate[1]*100;
  enddate[0]   = (int)(atof(argv[7])) - enddate[1]*1000000 
    - enddate[2]*10000;
  endjday      = juldate(enddate[0],enddate[1],enddate[2],0);

  for (l = 0; l < Nfiles ; l++) {

    /** Read Grid Cell File **/
    fscanf(files,"%s  %f  %f \n", file_index,&lat,&lon);    

    /** Create Grid Cell Flux Filename and Open for Reading **/
    strcpy(name,resultdir);
    strcat(name,model);
    strcat(name,file_index);
    puts(name);
    infile  = fopen(name, "r");
    if (infile == NULL){
      printf("File does not exist\n");
      exit(1);
    }
    printf("Ready for reading\n");

    /** Create Output Filename and Open for Writing **/
    strcpy(out_name,outdir);
    strcat(out_name,model);
    strcat(out_name,file_index);
    puts(out_name);
    outfile = fopen(out_name, "w");
    if (outfile == NULL){
      printf("File does not exist\n");
      exit(2);
    }

    /** Determine if files are Daily or sub-daily **/
    if(FIRSTFILE) {
      FIRSTFILE = FALSE;
      dt = get_file_time_step(&infile,&FROZEN);
      if(dt==24) {
	DAILY = TRUE;
	date[3] = 0;
	fprintf(stdout,"Model time step is daily.\n");
      }
      else {
	DAILY = FALSE;
	fprintf(stdout,"Model time step = %i hours.\n", dt);
      }
      if(FROZEN) {
	fprintf(stdout,"Frozen soil algorithm activated.\n");
	Ncols = NUM_DATA_VALS+2*NUM_LAYERS+2*NUM_FRONTS;
      }
      else Ncols = NUM_DATA_VALS+NUM_LAYERS;
      data = (float *)malloc(Ncols*sizeof(float));
      sum  = (float *)malloc(Ncols*sizeof(float));
      rewind(infile);
    }

    /** Process grid flux file **/
    DONE = FALSE;
    for(i=0;i<4;i++) lastdate[i] = -9999;
    get_record(&infile,&date,&data,DAILY,FROZEN);
    while(!feof(infile) && !DONE) {
      tmpjday = juldate(date[0],date[1],date[2],date[3]);
      if(tmpjday>=startjday && tmpjday<=endjday) {
	if(tmpjday==startjday || lastdate[0] == -9999) {
	  if(tmpjday != startjday) 
	    fprintf(stderr,"WARNING: first time step incomplete\n");
	  lastjday = startjday;
	  for(i=0;i<4;i++) lastdate[i] = startdate[i];
	  for(i=0;i<Ncols;i++) sum[i] = 0;
	  datacnt = 0;
	}
	if(lastjday==floor(tmpjday)) {
	  for(i=0;i<Ncols;i++) sum[i] += data[i];
	  datacnt ++;
	}
	else {
	  fprintf(outfile,"%d %d %d", lastdate[0], lastdate[1],
		  lastdate[2]);
	  for(i=0;i<4;i++) 
	    fprintf(outfile," %f", sum[i]);
	  for(i=4;i<Ncols;i++) fprintf(outfile," %f", sum[i]*((float)dt/24.));
	  fprintf(outfile,"\n");
	  lastjday = tmpjday;
	  for(i=0;i<3;i++) lastdate[i] = date[i];
	  lastdate[3] = 0;
	  for(i=0;i<Ncols;i++) sum[i] = data[i];
	  datacnt = 1;
	}
      }
      if(tmpjday > endjday) DONE = TRUE;
      get_record(&infile,&date,&data,DAILY,FROZEN);
    }
    if(datacnt > 0) {
      if(datacnt != 24 / dt) 
	fprintf(stderr,"WARNING: last time step incomplete.\n");
      fprintf(outfile,"%d %d %d", lastdate[0], lastdate[1],
	      lastdate[2]);
	  for(i=0;i<4;i++) 
	    fprintf(outfile," %f", sum[i]);
	  for(i=4;i<Ncols;i++) 
	    fprintf(outfile," %f", sum[i]*((float)dt/24.));
	  fprintf(outfile,"\n");
    }

    fclose(infile);
    fclose(outfile);
  }

  free((char *)date);
  free((char *)tmpdate);
  free((char *)lastdate);
  free((char *)enddate);
  free((char *)startdate);
  free((char *)data);
  free((char *)sum);

}
  
double juldate(int year,
	       int month,
	       int day,
	       double time)
/***********************************************************************
  juldate                     Keith Cherkauer             June 16, 1999

  This subroutine computes the julian date for the given decimal hour, 
  day, month, and year.  Julian dates from this routine start on 
  January 1, 1900.  Subtract 29220 to start on January 1, 1980, and 
  32873 to start on January 1, 1990.
***********************************************************************/
{
  double date, calc1, calc2, calc3, calc4;
  double hour, min, sec;

  hour = floor(time);
  time = time - (hour);
  min  = floor(time * 60.);
  time = time - (min / 60.);
  sec  = floor(time * 3600);
  time = time - (min / 3600);
  hour = hour + (min/60) + (sec/3600);

  calc1 = (367 * (double)year);
  calc2 = abs(floor(((double)month + 9)/12)); 
  calc2 = abs(floor(7*((double)year + calc2)/4));
  calc3 = abs(floor(275 * (double)month/9));
  calc4 = ((double)day + (hour/24.0) - 694006);
  
  date = calc1 - calc2 + calc3 + calc4;

  return(date);

}

#define FILEFORMATS 4

int get_file_time_step(FILE **fin, int *FROZEN) {

  int dt, idx;
  int syear, smonth, sday, shour;
  int eyear, emonth, eday, ehour;
  int FOUND, DAILY;
  int Nskip[FILEFORMATS] = { 34+2*NUM_LAYERS, 34+2*NUM_LAYERS, 34+4*NUM_LAYERS+4*NUM_FRONTS, 34+4*NUM_LAYERS+4*NUM_FRONTS };
  char tmpdata[34+4*NUM_LAYERS+4*NUM_FRONTS];
  char *cptr;
  unsigned short int *usiptr;
  
  cptr = (char *)malloc(1*sizeof(char));
  usiptr = (unsigned short int *)malloc(1*sizeof(unsigned short int));

  FOUND = FALSE;
  DAILY = TRUE;
  for(idx=0;idx<FILEFORMATS;idx++) {
    if(!FOUND) {
      rewind(fin[0]);
      /* year */
      fread(usiptr,1,sizeof(unsigned short int),fin[0]);
      syear = (int)usiptr[0];
      /* month */
      fread(cptr,1,sizeof(char),fin[0]);
      smonth = (int)cptr[0];
      /* day */
      fread(cptr,1,sizeof(char),fin[0]);
      sday = (int)cptr[0];
      if(!DAILY) {
        /* hour */
        fread(cptr,1,sizeof(char),fin[0]);
        shour = (int)cptr[0];
      }
      /* skip data */
      fread(tmpdata,Nskip[idx],sizeof(char),fin[0]);
      /* year */
      fread(usiptr,1,sizeof(unsigned short int),fin[0]);
      eyear = (int)usiptr[0];
      /* month */
      fread(cptr,1,sizeof(char),fin[0]);
      emonth = (int)cptr[0];
      /* day */
      fread(cptr,1,sizeof(char),fin[0]);
      eday = (int)cptr[0];
      if(!DAILY) {
        /* hour */
        fread(cptr,1,sizeof(char),fin[0]);
        ehour = (int)cptr[0];
      }
  
      printf("first date: %i/%i/%i - %i\nsecond date: %i/%i/%i - %i\n",smonth,sday,syear,shour,emonth,eday,eyear,ehour);
  
      if(syear==eyear || (syear+1==eyear)) {
        if(DAILY) {
	  dt = 24;
	  FOUND = TRUE;
        }
        else {
          dt = ehour - shour;
          if ( dt < 0 ) dt += 24;
	  FOUND = TRUE;
        }
        if(idx>=2) *FROZEN=TRUE;
        else *FROZEN=FALSE;
      }
    }
    if(DAILY) DAILY=FALSE;
    else DAILY=TRUE;
  }

  if(!FOUND) {
    fprintf(stderr,"ERROR: Unable to determine LADS flux file time step or format.\n");
    exit(0);
  }
  return(dt);

}

#undef FILEFORMATS

void get_record(FILE **fin, int **date, float **data, int DAILY, int FROZEN) {

  char *cptr;
  short int *siptr;
  unsigned short int *usiptr;
  int   *iptr;
  float *fptr;
  int i, j;
  
  cptr = (char *)malloc(1*sizeof(char));
  siptr = (short int *)malloc(1*sizeof(short int));
  usiptr = (unsigned short int *)malloc(1*sizeof(unsigned short int));
  iptr = (int *)malloc(1*sizeof(int));
  fptr = (float *)malloc(1*sizeof(float));

  /* year */
  fread(usiptr,1,sizeof(unsigned short int),fin[0]);
  if(!feof(fin[0])) {
    date[0][0] = (int)usiptr[0];
    /* month */
    fread(cptr,1,sizeof(char),fin[0]);
    date[0][1] = (int)cptr[0];
    /* day */
    fread(cptr,1,sizeof(char),fin[0]);
    date[0][2] = (int)cptr[0];
    /* hour */
    if(!DAILY) {
      fread(cptr,1,sizeof(char),fin[0]);
      date[0][3] = (int)cptr[0];
    }
    /* prec */
    fread(usiptr,1,sizeof(unsigned short int),fin[0]);
    data[0][0] = (float)usiptr[0]/100;
    /* evap */
    fread(siptr,1,sizeof(short int),fin[0]);
    data[0][1] = (float)siptr[0]/100;
    /* runoff */
    fread(fptr,1,sizeof(float),fin[0]);
    data[0][2] = fptr[0];
    /* baseflow */
    fread(fptr,1,sizeof(float),fin[0]);
    data[0][3] = fptr[0];
    /* moist */
    for(i=0;i<NUM_LAYERS;i++) {
      fread(usiptr,1,sizeof(unsigned short int),fin[0]);
      data[0][4+i] = (float)usiptr[0]/10;
    }
    /* swq */
    fread(usiptr,1,sizeof(unsigned short int),fin[0]);
    data[0][4+NUM_LAYERS] = (float)usiptr[0]/100;
    /* net_short */
    fread(siptr,1,sizeof(short int),fin[0]);
    data[0][5+NUM_LAYERS] = (float)siptr[0]/10;
    /* in_long */
    fread(siptr,1,sizeof(short int),fin[0]);
    data[0][6+NUM_LAYERS] = (float)siptr[0]/10;
    /* r_net */
    fread(siptr,1,sizeof(short int),fin[0]);
    data[0][7+NUM_LAYERS] = (float)siptr[0]/10;
    /* latent */
    fread(siptr,1,sizeof(short int),fin[0]);
    data[0][8+NUM_LAYERS] = (float)siptr[0]/10;
    /* sensible */
    fread(siptr,1,sizeof(short int),fin[0]);
    data[0][9+NUM_LAYERS] = (float)siptr[0]/10;
    /* grnd_flux */
    fread(siptr,1,sizeof(short int),fin[0]);
    data[0][10+NUM_LAYERS] = (float)siptr[0]/10;
    /* albedo */
    fread(usiptr,1,sizeof(unsigned short int),fin[0]);
    data[0][11+NUM_LAYERS] = (float)usiptr[0]/10000;
    /* surf_temp */
    fread(siptr,1,sizeof(short int),fin[0]);
    data[0][12+NUM_LAYERS] = (float)siptr[0]/100;
    /* rel_humid */
    fread(usiptr,1,sizeof(unsigned short int),fin[0]);
    data[0][13+NUM_LAYERS] = (float)usiptr[0]/100;
    /* air_temp */
    fread(siptr,1,sizeof(unsigned short int),fin[0]);
    data[0][14+NUM_LAYERS] = (float)siptr[0]/100;
    if(FROZEN) {
      /* ice */
      for(i=0;i<NUM_LAYERS;i++) {
        fread(usiptr,1,sizeof(unsigned short int),fin[0]);
        data[0][15+NUM_LAYERS+i] = (float)usiptr[0]/10;
      }
      for(i=0;i<NUM_FRONTS;i++) {
        fread(usiptr,1,sizeof(unsigned short int),fin[0]);
        data[0][15+2*NUM_LAYERS+2*i] = (float)usiptr[0]/100;
        fread(usiptr,1,sizeof(unsigned short int),fin[0]);
        data[0][16+2*NUM_LAYERS+2*i] = (float)usiptr[0]/100;
      }
    }
  }  

  free((char *)cptr);
  free((char *)siptr);
  free((char *)usiptr);
  free((char *)iptr);
  free((char *)fptr);

}

