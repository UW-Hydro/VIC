# VIC Forcing File

The VIC Image Driver requires a NetCDF file with gridded subdaily forcings. Forcing timestep must be the same as snow model timestep, which is specified by the `SNOW_STEPS_PER_DAY` parameter in the [Global Parameter File](GlobalParam.md). The required forcing variables and units are listed below and must also be specified in the [Global Parameter File](GlobalParam.md):

#### Meteorological Forcings, Required in all simulations:

| Variable   | Description                         | Units           |   
|------------|-------------------------------------|---------------- |
| AIR_TEMP   | Average air temperature             | C               |   
| PREC       | Total precipitation (rain and snow) | mm              |   
| PRESSURE   | Atmospheric pressure                | kPa             |   
| SWDOWN     | Incoming shortwave radiation        | W/m<sup>2</sup> |
| LWDOWN     | Incoming longwave radiation         | W/m<sup>2</sup> |
| VP         | Vapor pressure                      | kPa             |   
| WIND       | Wind speed                          | m/s             |   

The forcing data must be chunked by calendar year, with each NetCDF file named by the year, e.g. `prefix.$year.nc`.

Example output from `ncdump -h Stehekin_image_test.forcings_10days.1949` should look like this:

```
netcdf Stehekin_image_test.forcings_10days.1949 {
dimensions:
    time = 240 ;
    lon = 5 ;
    lat = 4 ;
variables:
    int time(time) ;
        time:units = "hours since 1949-01-01" ;
        time:calendar = "proleptic_gregorian" ;
    double lon(lon) ;
        lon:standard_name = "longitude" ;
        lon:long_name = "longitude of grid cell center" ;
        lon:units = "degrees_east" ;
        lon:axis = "X" ;
    double lat(lat) ;
        lat:standard_name = "latitude" ;
        lat:long_name = "latitude of grid cell center" ;
        lat:units = "degrees_north" ;
        lat:axis = "Y" ;
    float prcp(time, lat, lon) ;
        prcp:_FillValue = 9.96921e+36f ;
        prcp:long_name = "PREC" ;
        prcp:column = 0 ;
        prcp:units = "mm/step" ;
        prcp:description = "PREC" ;
    float tas(time, lat, lon) ;
        tas:_FillValue = 9.96921e+36f ;
        tas:long_name = "AIR_TEMP" ;
        tas:column = 1 ;
        tas:units = "C" ;
        tas:description = "AIR_TEMP" ;
    float dswrf(time, lat, lon) ;
        dswrf:_FillValue = 9.96921e+36f ;
        dswrf:long_name = "SWDOWN" ;
        dswrf:column = 2 ;
        dswrf:units = "W/m2" ;
        dswrf:description = "SWDOWN" ;
    float dlwrf(time, lat, lon) ;
        dlwrf:_FillValue = 9.96921e+36f ;
        dlwrf:long_name = "LWDOWN" ;
        dlwrf:column = 3 ;
        dlwrf:units = "W/m2" ;
        dlwrf:description = "LWDOWN" ;
    float pres(time, lat, lon) ;
        pres:_FillValue = 9.96921e+36f ;
        pres:long_name = "PRESSURE" ;
        pres:column = 5 ;
        pres:units = "kPa" ;
        pres:description = "PRESSURE" ;
    float vp(time, lat, lon) ;
        vp:_FillValue = 9.96921e+36f ;
        vp:long_name = "VP" ;
        vp:column = 6 ;
        vp:units = "kPa" ;
        vp:description = "VP" ;
    float wind(time, lat, lon) ;
        wind:_FillValue = 9.96921e+36f ;
        wind:long_name = "WIND" ;
        wind:column = 7 ;
        wind:units = "m/s" ;
        wind:description = "WIND" ;

// global attributes:
        :title = "Stehekin 10-day example (19490101-19490110)" ;
        :history = "Created: Fri May  6 09:34:06 2016 by ymao" ;
        :institution = "University of Washington" ;
        :source = "/home/ymao/.conda/envs/vic5/bin/vic_utils" ;
        :references = "Primary Historical Reference for VIC: Liang,X., D. P. Lettenmaier, E. F. Wood, and S. J. Burges,1994: A Simple hydrologically Based Model of LandSurface Water and Energy Fluxes for GSMs, J. Geophys.Res., 99(D7), 14,415-14,428." ;
        :comment = "Output from the Variable Infiltration Capacity(VIC) Macroscale Hydrologic Model" ;
        :conventions = "CF-1.6" ;
        :hostname = "compute-0-2.local" ;
        :username = "ymao" ;
        :version = "VIC.4.2.b" ;
        :grid = "BPA Grid" ;
}

```

#### Non-meteorological Forcings

Documentation for vegetation, lake and carbon cycle forcings will be added at a later date.
