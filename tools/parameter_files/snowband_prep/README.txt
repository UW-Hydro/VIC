README.txt

The following directories and files are contained in the tarball bands_prep.tgz:

Please note that snowband.c and elevband.c are two alternate programs for
creating the same file, created by two different people.  The user must decide
which to use.  


elevband.c
	
	prepares the snow elevation band input file based on the ASCII column
	format soil input file and an ASCII raster elevation file (with Arc/Info
	header) at a resolution finer than the resolution of the VIC model run

snowband.c

	Purpose: Calculates snow band elevations and fractions for grid cells
	used in VIC, based on a finer resolution DEM, and PRISM precip scalings for 
	precip fractions in each snow band.


run_snowband.scr

	compiles snowband.c and runs it to make a snowband parameter file

