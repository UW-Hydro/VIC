README.txt

The following directories and files are contained in the tarball veg_param_prep.tgz:

make_veg.aml.gz

	generates one ASCII grid for each vegetation class.  In each gridcell
the areal cover as a fraction is given.  If you use 20 vegetation classes
output from thr program is 20 ASCII grids..

make_vegetation.c.gz

	the output file, output_veg_paramfile, is the vegetation parameter
file that the VIC model needs

LDAS_2_vegparam.c

	the raw LDAS vegetation data can be converted to a VIC-style with the
program LDAS_2_vegparam.c

aggregate_LDAS_veg_param.c

	for the vegetation parameter file format including LAI, aggregates the
file
