#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
 
dmy_struct *make_dmy(global_param_struct global)
/**********************************************************************
	make_dmy	Dag Lohmann		January 1996

  This subroutine creates an array of structures that contain 
  information about the day, month and year of each time step.

  modifications:
  7-25-96  Added hour count, so that model can run on less than
           a daily time step.					KAC

NOTE need to modify code so that number of records will be total times
model is run, while number of days is used to reference daily data.
number of records can then be used to track hourly time steps.
**********************************************************************/
{
  dmy_struct *temp;
  int    hr, i, j, k, l, ii, daymax;
  int    days[12]={31,28,31,30,31,30,31,31,30,31,30,31};

  temp = (dmy_struct*) calloc(global.nrecs, sizeof(dmy_struct));

  hr = global.starthour;
  i = global.startyear;
  j = global.startday;
  k = global.startmonth;

  l = j;
  for (ii=0;ii<k-1;ii++) 
    l += days[ii];
  if ((k > 2 ) && ( (i%4 == 0) && ( (i%100 != 0) || (i%400 ==0) )))
    l += 1;

  for(ii=0;ii<global.nrecs;ii++) {
    temp[ii].hour = hr;
    temp[ii].day   = j;
    temp[ii].month = k;
    temp[ii].year  = i;
    temp[ii].day_in_year = l;
    temp[ii].day_count = ii;

    hr += global.dt;
    if(hr >= 24) {
      hr=0;
      j++;
      l++;

      daymax = days[k-1];
      if ( (k == 2) && ( (i%4 == 0) && ( (i%100 != 0) || (i%400 ==0))))
          daymax = 29;
  
      if(j > daymax) {
        j = 1;
        k++;
        if(k == 13){
          k = 1;
	  l = 1;
          i++;
        }
      } 
    }
  }
  return temp;
}

