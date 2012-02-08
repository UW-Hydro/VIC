#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

static char vcid[] = "$Id$";

void write_snow_data(snow_data_struct snow,
		     int              veg,
		     int              band)
/**********************************************************************
  write_snow_data		Keith Cherkauer		October 26, 2000

  This routine writes all values stored in the snow data structure.

  Modifications:

  2012-Feb-08 Renamed depth_full_snow_cover to max_snow_distrib_slope
	      and clarified the descriptions of the SPATIAL_SNOW
	      option.							TJB
**********************************************************************/
{
  printf("Snowpack Data Variables: veg type %i, snow band %i\n", veg, band);
  printf("snow              = %x\n", snow.snow);
  printf("Qnet              = %f W/m^2\n", snow.Qnet);
  printf("albedo            = %f\n", snow.albedo);
  printf("canopy_vapor_flux = %f m\n", snow.canopy_vapor_flux);
  printf("coldcontent       = %f\n", snow.coldcontent);
  printf("coverage          = %f\n", snow.coverage);
  printf("density           = %f kg/m^3\n", snow.density);
  printf("depth             = %f m\n", snow.depth);
  printf("mass_error        = %f W/m^2\n", snow.mass_error);
  printf("max_snow_depth    = %f m\n", snow.max_snow_depth);
  printf("pack_temp         = %f C\n", snow.pack_temp);
  printf("pack_water        = %f m\n", snow.pack_water);
  printf("snow_canopy       = %f m\n", snow.snow_canopy);
  printf("surf_temp         = %f C\n", snow.surf_temp);
  printf("surf_water        = %f m\n", snow.surf_water);
  printf("swq               = %f m\n", snow.swq);
  printf("tmp_int_storage   = %f m\n", snow.tmp_int_storage);
  printf("vapor_flux        = %f m\n", snow.vapor_flux);
  printf("last_snow         = %i\n", snow.last_snow);

}

