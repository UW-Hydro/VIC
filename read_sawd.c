#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
 
void read_sawd(atmos_data_struct *temp,
               FILE              *sawdf,
               int               *nrecs,
               int                dt,
               int                prec)
/**********************************************************************
	read_sawd	Keith Cherkauer		July 25, 1996

  This routine reads in atmospheric data values from surface
  airways gridded hourly data files.  If time step is less than
  hourly, extra data is skipped.

				Input	Output
				Units	Units
	air_temp	-	C	C
	rel_humid	-	%	%
	tskc		-	%	fraction
	wind		-	m/s	m/s
	pressure	-	Pa	kPa
	prec		-	mm	mm

**********************************************************************/
{
  extern option_struct options;
  extern debug_struct debug;
  extern param_set_struct param_set;
 
  int    i, j, n, rec, maxline = 210;
  int    fixcnt, headcnt;
  char   str[210],jnkstr[MAXSTRING];

  /** Count Records **/
  n = 0;
  while (fgets(str,maxline,sawdf) != '\0') n++;
  printf("nrecs = %d\n",n);
  if(*nrecs*dt>n) {
    fprintf(stderr,"WARNING: SAWD file does not have as many records as defined in the global parameter file, truncating run to %i records.\n",n);
    *nrecs=n/dt;
  }
  if(*nrecs*dt<n) {
    fprintf(stderr,"WARNING: SAWD file has more records then were defined in the global parameter file, run will stop after %i records.\n",*nrecs);
    n = *nrecs*dt;
  }

  rewind(sawdf);

  /** Check for Header, and Skip **/
  fgets(jnkstr,MAXSTRING,sawdf);
  headcnt=0;
  if((jnkstr[0]<48 || jnkstr[0]>57) && jnkstr[0]!=46) {
    fgets(jnkstr,MAXSTRING,sawdf);
    fprintf(stderr,"SAWD ... skipping header\n");
    n--;
    headcnt++;
    if(prec) *nrecs = n;
  }
  rewind(sawdf);
  for(i=0;i<headcnt;i++) fgets(jnkstr,MAXSTRING,sawdf);

  fixcnt = 0;
  rec = 0;
  for (i = 0; i < n; i++) {
    if(i == rec*dt) {
      fscanf(sawdf,"%*s");
      fscanf(sawdf,"%s",str);
      temp[rec].air_temp = atof(str);
      fscanf(sawdf,"%s",str);
      temp[rec].pressure = atof(str) / 1000.;
      if(temp[rec].pressure <= 0. && rec>0) {
	temp[rec].pressure = temp[rec-1].pressure;
	fixcnt++;
      }
      fscanf(sawdf,"%s",str);
      temp[rec].wind = atof(str);
      fscanf(sawdf,"%s",str);
      temp[rec].rel_humid = atof(str);
      if(temp[rec].rel_humid == 0. && rec>0) {
	temp[rec].rel_humid = temp[rec-1].rel_humid;
      fixcnt++;
      }
      fscanf(sawdf,"%s",str);
      temp[rec].tskc = atof(str) / 100.;
      if(prec) {
	fscanf(sawdf,"%s",str);
	temp[rec].prec = atof(str);
      }
      rec++;
    }
    else fgets(jnkstr,MAXSTRING,sawdf);

  }

  if(fixcnt>0) {
    fprintf(stderr,"WARNING: Had to fix %i values in sawd file.\n",fixcnt);
  }

  param_set.WIND = param_set.AIR_TEMP = param_set.PRESSURE = TRUE;
  param_set.TSKC = param_set.REL_HUMID = TRUE;
  if(prec) param_set.PREC = TRUE;

}
