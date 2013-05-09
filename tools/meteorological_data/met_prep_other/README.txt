The following is the list of programs included in the met_prep_other.tgz tarball.  

ascend.c

	 This program reads a mask file and generates a new maskfile with unique integers in each valid cell.
 
format_wind.scr

	cleans the wind data file and station file before simple gridding program is used

get_prism_airtemp.aml

	 Script used to calculate airtemp PRISM monthly means in each gridcell output from this script will be used to scale airtemp data.


implement_wind_vic.c

	This program reads monthly winddata from Surface Airways CDs. The rawdata have to be formatted with format_wind.scr
	before this program can read the windfields. The VIC inputdata (Prcp, Tmax,Tmin) is read and
  	a daily timeseries of wind is added as a fourth  coloumn to the metdata.


mk_monthly_airtemp.c

	This program reads the output file from the regridding program. It then generates monthly values for jan,feb,mar,etc
	that will be used together with PRISM data to scale the airtemperature timeseries.  

rescale_airtemp.c

	This program rescales the gridded airtemperature data from the regrid program 
	with Tmax = Tmax + TmaxPrism -TmaxMonth. This is done to get the long term mean airtemperatures to be the same as PRISM 

run_mk_monthly_airtemp.c

	runs mk_monthly_airtemp.c

run_rescale_airtemp.c

	runs rescale_airtemp.c
