#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
 
#ifndef _LEAPYR
#define LEAPYR(y) (!((y)%400) || (!((y)%4) && ((y)%100)))
#endif

static char vcid[] = "$Id$";

dmy_struct *make_dmy(global_param_struct *global)
/**********************************************************************
	make_dmy	Dag Lohmann		January 1996

  This subroutine creates an array of structures that contain 
  information about the day, month and year of each time step.

  modifications:
  7-25-96  Added hour count, so that model can run on less than
           a daily time step.					KAC
  5-17-99  Modified routine to use LEAPYR function, make use of 
           simulation ending dates, and to skip over initial 
	  forcing data records so that the model can be run on
	  subsets of more complete data records.            KAC

**********************************************************************/
{
  dmy_struct *temp;
  int    hr, year, day, month, jday, ii, daymax;
  int    days[12]={31,28,31,30,31,30,31,31,30,31,30,31};
  int    endmonth, endday, endyear;
  int    tmpmonth, tmpday, tmpyear, tmphr, tmpjday, step;
  char   DONE;

  hr = global->starthour;
  year = global->startyear;
  day = global->startday;
  month = global->startmonth;
  
  /** Check if user defined end date instead of number of records **/
  if(global->nrecs<0) {
    if(global->endyear<0 || global->endmonth<0 || global->endday<0) {
      nrerror("The model global file MUST define EITHER the number of records to simulate (NRECS), or the year (ENDYEAR), month (ENDMONTH), and day (ENDDAY) of the last full simulation day");
    }
    endday   = global->endday;
    endmonth = global->endmonth;
    endyear  = global->endyear;
    if(( (endyear%4 == 0) && ( (endyear%100 != 0) || (endyear%400 ==0) ))) days[1] = 29;
    if(endday < days[global->endmonth-1]-1) endday++;
    else {
      endday = 1;
      endmonth++;
      if(endmonth > 12) {
	endmonth = 1;
	endyear++;
      }
    }

    DONE = FALSE;
    ii   = 0;

    tmpyear = year;
    tmpmonth = month;
    tmpday = day;
    tmphr = hr;
    while(!DONE) {
      get_next_time_step(&tmpyear,&tmpmonth,&tmpday,&tmphr,&tmpjday,global->dt);
      ii++;
      if(tmpyear == endyear)
	if(tmpmonth == endmonth)
	  if(tmpday == endday)
	    DONE = TRUE;
    }
    global->nrecs = ii;

  }

  temp = (dmy_struct*) calloc(global->nrecs, sizeof(dmy_struct));

  /** Create Date Structure for each Modeled Time Step **/
  jday = day;
  for (ii=0;ii<month-1;ii++) 
    jday += days[ii];
  if ((month > 2 ) && LEAPYR(year))
    jday += 1;
  
  DONE = FALSE;
  ii   = 0;
  
  while(!DONE) {
    temp[ii].hour = hr;
    temp[ii].day   = day;
    temp[ii].month = month;
    temp[ii].year  = year;
    temp[ii].day_in_year = jday;
    temp[ii].day_count = ii;

    get_next_time_step(&year,&month,&day,&hr,&jday,global->dt);

    ii++;
    if(global->nrecs>0 && ii==global->nrecs) DONE=TRUE;
    else if(year == endyear)
      if(month == endmonth)
	if(day == endday)
	  DONE = TRUE;

  }

  if(global->nrecs < 0) global->nrecs = ii;

  /** Determine number of forcing records to skip before model start time **/
  if(global->forceyear > 0) {
    tmpyear  = global->forceyear;
    tmpmonth = global->forcemonth;
    tmpday   = global->forceday;
    tmphr    = global->forcehour;
    tmpjday = tmpday;
    for (ii=0;ii<tmpmonth-1;ii++) 
      tmpjday += days[ii];
    if ((tmpmonth > 2 ) && LEAPYR(tmpyear))
      tmpjday += 1;

    step     = (int)(1./((float)global->dt/24.));
    while(tmpyear < temp[0].year || 
	  (tmpyear == temp[0].year && tmpjday < temp[0].day_in_year)) {
      
      get_next_time_step(&tmpyear,&tmpmonth,&tmpday,&tmphr,
			 &tmpjday,global->dt);
      
      global->forceskip ++;
    }
  }

  return temp;
}

void get_next_time_step(int *year, 
			int *month, 
			int *day, 
			int *hr, 
			int *jday, 
			int dt) {
  
  int    days[12]={31,28,31,30,31,30,31,31,30,31,30,31};
  int daymax;
  
  *hr += dt;
  if(*hr >= 24) {
    *hr=0;
    *day += 1;
    *jday += 1;
    
    if(LEAPYR(*year)) days[1] = 29;
    else days[1] = 28;
    
    if(*day > days[*month-1]) {
      *day = 1;
      *month += 1;
      if(*month == 13){
	*month = 1;
	*jday  = 1;
	*year += 1;
      }
    } 
  }
  
}
