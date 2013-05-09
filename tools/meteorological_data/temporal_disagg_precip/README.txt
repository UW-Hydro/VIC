README.txt

The following files are contained in the tarball temporal_disagg_precip.tgz:

disagg_prec.c

	disaggregates daily VIC precip rfom a 4 column binary input data file

make_lookup_table.c

	creates a set of look up tables (CDFs) for precipitation duration and
start-time/occurence time on an hourly basis, by month

reformat_prcp.scr

	reformats hourly precipitation retrieved from the EarthInfo CDs in
ASCII format. Produces a station file and a data file with hourly data

run_disagg_prec.scr

	runs disagg_prec.c

run_make_lookup_table.scr

	runs make_lookup_table.c
