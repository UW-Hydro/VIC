README.txt

The following directories and files are contained in the tarball soil_param_prep.tgz:

standard/
	Directory containing files for preparation of standard soil parameter files

	source.usda.triangle.zip

		calls a function that classifies soil in the USDA textual
		triangle using sand and clay %

	run_awk_nv.scr

		creaets hydraulic parameters according to the "cosby.1984.csv"
		look up table

	texture2prop.awk

	cosby.1984.csv

		assign hydraulic parameter to soil texture (grouped with
		run_awk_nv.scr and texture2prop.awk)

	make_parametersfiles.scr
		
		makes Ks.txt, porosity.txt, fieldCap.txt and WiltPint.txt and
		b.txt, each with 3 columns for each soil layer, based on the created
		para_hydro.txt

	clipSoilVariable.c

	clipSoilVariableClass.c

	make_clipsoil.scr

	create_soilf.c

		to read soil parameters from various soil file directories and
		put them together as the final soil file (VIC version 4.0.1)

	fix_resid.pl

		quality control on soil file

	set_moist.pl

		sets the initial soil moisture based on soil characteristics

	convert_ARCINFO_to_ASCII_soil_file.c

		reads a runfile mask and an Arc/Info soil parameter control
		file and creates an ASCII column soil data file for the grid cells activated
		in the masks file

	Scalesrf.c
		
		aggregate 5 minute .srf files produced by the soil data
		program to other resolutions (multiple of 5 minutes)

	aggr_soil_2vic.c

		aggregate soil data from arc30 resolution into 1/8 VIC grid
		cell

	convert_srf_to_text.scr

		script to remove the carriage return in the srf files

	get_van_GN.c

		to get the brooks-corey or van Genuchten N from the sand and
		clay percentages; determine lambda using eqns on pg.5.15 in handbook of
		hydrologyl; use the lambda to get N using eqns on pg. 5.6

		modified by N. Voisin to handle the clipped sandclay files directly 
		(output from the usda.triangle that has been clipped out )
		note that the porosity must be in fraction format, NOT percent!

	make_5min_II_data.scr

		this script is a copy of  make_arc30_ll_data.scr  but I want
		to keep the 5 minutes resolution as this is my final
		resolution

	make_arc30_II_data.scr

		# This script is used to make lon lat data file in arc30
		# resoution for *.srf file, which dissegregate the 5-minutes
		# *.srf file into 30second resolution. This output will be
		# used to later aggregation into 1/8 degree VIC cell.   

	make_my_ll_data.scr

	run_soil_aggr.scr

	





