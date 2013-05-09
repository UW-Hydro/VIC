README.txt

The following files are contained in the tarball temporal_agg.tgz:

make_file_list.script

	To simplify the process of making a file list when you have a large number of grid cells

extract_daily_means_from_LDAS_binary.c

	determines the time step of the LDAS binary file and whether or not it
contains frozen soil information.  Then it reads through the file and
aggregates variables from sub-daily to daily values. Precipitation, runoff,
evaporation, and baseflow are summed to give daily totals.  All other
variables are averaged to produce daily average values for each day.  

agg_time.pl

	perl script that will aggregate any ascii column-formatted data from
one time interval to a larger one, ranging from hourly through daily, monthly,
and annual.  “agg_time.pl -h” will give a comprehensive usage message.  (you
might want to paste the us
