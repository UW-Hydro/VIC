SHELL = /bin/csh

CC = gcc -I. -g
#CC = cc -I.
#LIBRARY = -lm
LIBRARY = -lm /nfs/hydro4/usr5/cherkaue/source/electric_fence_meter/lib/libefence.a
#LIBRARY = -lm -lc
UTILDIR = ./
DISTDIR = ./
FULLDIR = ./
FROZDIR = ./

INCFILE = vicNl.h vicNl_def.h global.h rad_and_vpd.h snow.h user_def.h

OBJS = vicNl.o cmd_proc.o check_files.o read_soilparam.o open_file.o \
       make_in_and_outfiles.o read_atmosdata.o read_vegparam.o make_dmy.o \
       get_global_param.o initialize_soil.o close_files.o rad_and_vpd.o \
       calc_trans.o calc_netshort.o aurad.o priestley.o svp.o \
       long_shortwave.o nrerror.o vicerror.o \
       make_veg_var.o make_cell_data.o initialize_veg.o penman.o arno_evap.o \
       put_data.o polint.o initialize_global.o read_snowband.o \
       read_snowmodel.o read_sawd.o read_sawd_binary.o \
       initialize_atmos.o correct_precip.o compress_files.o \
       write_soilparam.o write_vegparam.o write_atmosdata.o write_data.o \
       write_vegvar.o write_layer.o make_dist_prcp.o free_dist_prcp.o \
       initialize_energy.o frozen_soil.o \
       soil_conduction.o make_energy_bal.o canopy_evap.o runoff.o \
       modify_Ksat.o full_energy.o func_surf_energy_bal.o read_rosemount.o \
       initialize_snow.o snow_melt.o root_brent.o calc_water_density.o \
       SnowPackEnergyBalance.o make_snow_data.o StabilityCorrection.o \
       snow_utility.o func_snow_ground_flux.o soil_thermal_eqn.o \
       calc_surf_energy_bal.o write_debug.o dist_prec.o open_debug.o \
       initialize_new_storm.o redistribute_during_storm.o \
       snow_intercept.o massrelease.o CalcAerodynamic.o \
       calc_veg_params.o read_veglib.o calc_rainonly.o \
       calc_long_shortwave.o calc_air_temperature.o read_PILPS2c.o \
       calc_snow_ground_flux.o read_forcing_data.o read_soilparam_arc.o \
       read_arcinfo_ascii.o solve_snow.o


all:
	make clean
	make model

default:
	make model

clean:
	/bin/rm -f *.o core log

model: $(OBJS)
	$(CC) -o vicNl $(OBJS) $(LIBRARY)

vicNl.o:  vicNl.c     $(INCFILE)
	$(CC) -c vicNl.c
cmd_proc.o:     	 cmd_proc.c $(INCFILE)
	$(CC) -c cmd_proc.c
runoff.o:   	  runoff.c $(INCFILE)
	$(CC) -c runoff.c
initialize_global.o:   	  initialize_global.c $(INCFILE)
	$(CC) -c initialize_global.c

# DATA HANDLING CODE
make_veg_var.o:   	make_veg_var.c $(INCFILE)
	$(CC) -c make_veg_var.c
make_cell_data.o:   	make_cell_data.c $(INCFILE)
	$(CC) -c make_cell_data.c
make_dmy.o:   	  make_dmy.c $(INCFILE)
	$(CC) -c make_dmy.c
initialize_soil.o:   	  initialize_soil.c $(INCFILE)
	$(CC) -c initialize_soil.c
initialize_veg.o:   	  initialize_veg.c $(INCFILE)
	$(CC) -c initialize_veg.c
initialize_atmos.o:   	  initialize_atmos.c $(INCFILE)
	$(CC) -c initialize_atmos.c
write_vegvar.o:	write_vegvar.c $(INCFILE)
	$(CC) -c write_vegvar.c
write_layer.o:	write_layer.c $(INCFILE)
	$(CC) -c write_layer.c

# INPUT/OUTPUT CODE
read_forcing_data.o:   	read_forcing_data.c $(INCFILE)
	$(CC) -c read_forcing_data.c
read_soilparam.o:   	read_soilparam.c $(INCFILE)
	$(CC) -c read_soilparam.c
read_soilparam_arc.o:   	read_soilparam_arc.c $(INCFILE)
	$(CC) -c read_soilparam_arc.c
read_arcinfo_ascii.o:   	read_arcinfo_ascii.c $(INCFILE)
	$(CC) -c read_arcinfo_ascii.c
write_soilparam.o:   	write_soilparam.c $(INCFILE)
	$(CC) -c write_soilparam.c
open_file.o:		open_file.c $(INCFILE)
	$(CC) -c open_file.c
open_debug.o:		open_debug.c $(INCFILE)
	$(CC) -c open_debug.c
read_vegparam.o:   	read_vegparam.c $(INCFILE)
	$(CC) -c read_vegparam.c
read_veglib.o:   	read_veglib.c $(INCFILE)
	$(CC) -c read_veglib.c
write_vegparam.o:   	write_vegparam.c $(INCFILE)
	$(CC) -c write_vegparam.c
make_in_and_outfiles.o:	make_in_and_outfiles.c $(INCFILE)
	$(CC) -c make_in_and_outfiles.c
read_atmosdata.o:   	read_atmosdata.c $(INCFILE)
	$(CC) -c read_atmosdata.c
read_snowmodel.o:   	read_snowmodel.c $(INCFILE)
	$(CC) -c read_snowmodel.c
read_snowband.o:   	read_snowband.c $(INCFILE)
	$(CC) -c read_snowband.c
read_sawd.o:   	read_sawd.c $(INCFILE)
	$(CC) -c read_sawd.c
read_sawd_binary.o:   	read_sawd_binary.c $(INCFILE)
	$(CC) -c read_sawd_binary.c
read_PILPS2c.o:   	read_PILPS2c.c $(INCFILE)
	$(CC) -c read_PILPS2c.c
write_atmosdata.o:   	write_atmosdata.c $(INCFILE)
	$(CC) -c write_atmosdata.c
get_global_param.o:   	  get_global_param.c $(INCFILE)
	$(CC) -c get_global_param.c
close_files.o:   	  close_files.c $(INCFILE)
	$(CC) -c close_files.c
put_data.o:   	  put_data.c  $(INCFILE)
	$(CC) -c put_data.c
write_data.o:	  write_data.c  $(INCFILE)
	$(CC) -c write_data.c
write_debug.o:	  write_debug.c  $(INCFILE)
	$(CC) -c write_debug.c
check_files.o:   	check_files.c $(INCFILE)
	$(CC) -c check_files.c

# RADIATION AND EVAPORATION CODE
rad_and_vpd.o:   	  rad_and_vpd.c $(INCFILE)
	$(CC) -c rad_and_vpd.c
calc_trans.o:   	  calc_trans.c $(INCFILE)
	$(CC) -c calc_trans.c
calc_netshort.o:   	  calc_netshort.c $(INCFILE)
	$(CC) -c calc_netshort.c
aurad.o:   	  aurad.c $(INCFILE)
	$(CC) -c aurad.c
priestley.o:   	  priestley.c $(INCFILE)
	$(CC) -c priestley.c
svp.o:   	  svp.c $(INCFILE)
	$(CC) -c svp.c
long_shortwave.o: long_shortwave.c $(INCFILE)
	$(CC) -c long_shortwave.c
canopy_evap.o:    canopy_evap.c $(INCFILE)
	$(CC) -c canopy_evap.c
penman.o:   	  penman.c $(INCFILE)
	$(CC) -c penman.c
arno_evap.o:   	  arno_evap.c $(INCFILE)
	$(CC) -c arno_evap.c

# DITRIBUTED PRECIPITATION CODE
polint.o:	polint.c $(INCFILE)
	$(CC) -c polint.c
nrerror.o:	nrerror.c $(INCFILE)
	$(CC) -c nrerror.c
vicerror.o:	vicerror.c $(INCFILE)
	$(CC) -c vicerror.c
make_prcp_var.o:   	  make_prcp_var.c $(INCFILE)
	$(CC) -c make_prcp_var.c
make_dist_prcp.o:	  make_dist_prcp.c  $(INCFILE)
	$(CC) -c make_dist_prcp.c
free_dist_prcp.o:	  free_dist_prcp.c  $(INCFILE)
	$(CC) -c free_dist_prcp.c
initialize_new_storm.o:	  initialize_new_storm.c  $(INCFILE)
	$(CC) -c initialize_new_storm.c
redistribute_during_storm.o:	  redistribute_during_storm.c  $(INCFILE)
	$(CC) -c redistribute_during_storm.c

# SNOW MELT CODE
make_snow_data.o:	make_snow_data.c $(INCFILE)
	$(CC) -c make_snow_data.c
initialize_snow.o:	initialize_snow.c $(INCFILE)
	$(CC) -c initialize_snow.c
solve_snow.o:	solve_snow.c $(INCFILE)
	$(CC) -c solve_snow.c
snow_melt.o:	snow_melt.c $(INCFILE)
	$(CC) -c snow_melt.c
root_brent.o:	root_brent.c $(INCFILE)
	$(CC) -c root_brent.c
calc_water_density.o:	calc_water_density.c $(INCFILE)
	$(CC) -c calc_water_density.c
SnowPackEnergyBalance.o:	SnowPackEnergyBalance.c $(INCFILE)
	$(CC) -c SnowPackEnergyBalance.c
StabilityCorrection.o:	StabilityCorrection.c $(INCFILE)
	$(CC) -c StabilityCorrection.c
func_snow_ground_flux.o:	func_snow_ground_flux.c $(INCFILE)
	$(CC) -c func_snow_ground_flux.c
calc_snow_ground_flux.o:	calc_snow_ground_flux.c $(INCFILE)
	$(CC) -c calc_snow_ground_flux.c
snow_utility.o:	snow_utility.c $(INCFILE)
	$(CC) -c snow_utility.c
soil_thermal_eqn.o:	soil_thermal_eqn.c $(INCFILE)
	$(CC) -c soil_thermal_eqn.c
snow_intercept.o:	snow_intercept.c $(INCFILE)
	$(CC) -c snow_intercept.c
massrelease.o:	massrelease.c $(INCFILE)
	$(CC) -c massrelease.c

# FROZEN SOIL CODE
frozen_soil.o:	frozen_soil.c $(INCFILE)
	$(CC) -c frozen_soil.c
initialize_energy.o:	initialize_energy.c $(INCFILE)
	$(CC) -c initialize_energy.c
soil_conduction.o:	soil_conduction.c $(INCFILE)
	$(CC) -c soil_conduction.c
make_energy_bal.o:	make_energy_bal.c $(INCFILE)
	$(CC) -c make_energy_bal.c
modify_Ksat.o:	modify_Ksat.c $(INCFILE)
	$(CC) -c modify_Ksat.c

# FULL ENERGY BALANCE CODE
dist_prec.o: dist_prec.c $(INCFILE)
	$(CC) -c dist_prec.c
full_energy.o: full_energy.c $(INCFILE)
	$(CC) -c full_energy.c
func_surf_energy_bal.o: func_surf_energy_bal.c $(INCFILE)
	$(CC) -c func_surf_energy_bal.c
calc_surf_energy_bal.o:	calc_surf_energy_bal.c $(INCFILE)
	$(CC) -c calc_surf_energy_bal.c
read_rosemount.o: read_rosemount.c $(INCFILE)
	$(CC) -c read_rosemount.c
correct_precip.o:	correct_precip.c $(INCFILE)
	$(CC) -c correct_precip.c
compress_files.o:	compress_files.c $(INCFILE)
	$(CC) -c compress_files.c
CalcAerodynamic.o:	CalcAerodynamic.c $(INCFILE)
	$(CC) -c CalcAerodynamic.c
calc_veg_params.o: calc_veg_params.c $(INCFILE)
	$(CC) -c calc_veg_params.c
calc_rainonly.o: calc_rainonly.c $(INCFILE)
	$(CC) -c calc_rainonly.c
calc_air_temperature.o: calc_air_temperature.c $(INCFILE)
	$(CC) -c calc_air_temperature.c
calc_long_shortwave.o: calc_long_shortwave.c $(INCFILE)
	$(CC) -c calc_long_shortwave.c
