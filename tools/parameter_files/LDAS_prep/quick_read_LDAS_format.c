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
#define TRUE !FALSE

int get_file_time_step(FILE **, int *);
void get_record(FILE **, int **, float **,int,int);

main(int argc, char *argv[]) {
/*****************************************************************
  Quick read of the binary LDAS flux file format, output sent to
  stdout as ASCII columns.

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

*****************************************************************/

  FILE *fin;
  char *cptr;
  short int *siptr;
  unsigned short int *usiptr;
  int   *iptr, *date, *tmpdate;
  int    DAILY;
  int    FROZEN;
  int    Ncols;
  float *fptr, *data;
  int i, j, dt;
  
  cptr = (char *)malloc(1*sizeof(char));
  siptr = (short int *)malloc(1*sizeof(short int));
  usiptr = (unsigned short int *)malloc(1*sizeof(unsigned short int));
  iptr = (int *)malloc(1*sizeof(int));
  fptr = (float *)malloc(1*sizeof(float));
  date = (int *)malloc(4*sizeof(int));
  tmpdate = (int *)malloc(4*sizeof(int));
 
  if(argc!=2) {
    fprintf(stderr,"Usage: %s <flux file name>\n",argv[0]);
    exit(0);
  }

  if((fin=fopen(argv[1],"rb"))==NULL) {
    fprintf(stderr,"ERROR: Unable to open input file %s.\n",argv[1]);
    exit(0);
  }

  /** Check output file time step **/
  dt = get_file_time_step(&fin,&FROZEN);
  if(dt==24) {
    DAILY = TRUE;
    date[3] = 0;
    fprintf(stderr,"Model time step is daily.\n");
  }
  else {
    DAILY = FALSE;
    fprintf(stderr,"Model time step = %i hours.\n", dt);
  }
  if(FROZEN) {
    fprintf(stderr,"Frozen soils activated.\n");
    Ncols = NUM_DATA_VALS+2*NUM_LAYERS+2*NUM_FRONTS;
  }
  else Ncols = NUM_DATA_VALS+NUM_LAYERS;
  data = (float *)malloc((Ncols)*sizeof(float));
  rewind(fin);

  get_record(&fin,&date,&data,DAILY,FROZEN);
 
  while(!feof(fin)) {

    if(DAILY) for(i=0;i<3;i++) fprintf(stdout,"%i\t",date[i]);
    else for(i=0;i<4;i++) fprintf(stdout,"%i\t",date[i]);
    for(i=0;i<Ncols-1;i++) fprintf(stdout,"%f\t",data[i]);
    fprintf(stdout,"%f\n",data[Ncols-1]);

    get_record(&fin,&date,&data,DAILY,FROZEN);

  }

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

  fprintf(stderr,"Bytes per file format: %i, %i, %i, %i\n",Nskip[0],Nskip[1],Nskip[2],Nskip[3]);

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
  
      fprintf(stderr,"first date: %i/%i/%i - %i\nsecond date: %i/%i/%i - %i\n",smonth,sday,syear,shour,emonth,eday,eyear,ehour);
  
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
        data[0][15+2*NUM_LAYERS+2*i+1] = (float)usiptr[0]/100;
      }
    }
  }  

  free((char *)cptr);
  free((char *)siptr);
  free((char *)usiptr);
  free((char *)iptr);
  free((char *)fptr);

}

