# VIC Domain file

The Image Driver uses the [NetCDF](http://www.unidata.ucar.edu/software/netcdf/) file format to define model running domain.

Below is a list of variables in the domain netCDF file. The dimensions of the netCDF file are `lat` and `lon`. Note that here only the type of variables (i.e., MASK, AREA, FRAC, LAT and LON) is listed; corresponding variable names in the input netCDF file are specified by user in the [Global Parameter File](GlobalParam.md). All the listed variables are required.

| Variable   | Dimension   | Units    | Type   | Description |
|------------|-------------|----------|--------|-------------|
| LAT        | [lat]       | degree   | double | Latitudes   |
| LON        | [lon]       | degree   | double | Longitues   |
| MASK       | [lat, lon]  | N/A      | integer | Mask of domain. 1 for grid cells inside considered domain; 0 for grid cells outside of domain. Cells outside of domain will not be run. Use run_cell variable in the parameter file to turn on/off active cells inside domain. |
| AREA       | [lat, lon]  | m2       | double | Area of grid cells.   |
| FRAC       | [lat, lon]  | N/A      | double | Fraction of grid cells that is land. |

# Example netCDF format VIC 5 image driver domain file

```shell
ncdump -h /ArkRed.domain.nc
netcdf ArkRed.domain {                                                                       
dimensions:                                                                                  
        lat = 66 ;                                                                           
        lon = 125 ;                                                                          
variables:                                                                                   
        int mask(lat, lon) ;
                mask:comment = "0 indicates grid cell outside of domain" ;
                mask:long_name = "domain mask" ;
        double lon(lon) ;
                lon:long_name = "longitude coordinate" ;
                lon:units = "degrees_east" ;
        double lat(lat) ;
                lat:long_name = "latitude coordinate" ;
                lat:units = "degrees_north" ;
        double frac(lat, lon) ;
                frac:long_name = "fraction of grid cell that is active" ;
                frac:units = "1" ;
        double area(lat, lon) ;
                area:standard_name = "area" ;
                area:long_name = "area of grid cell" ;
                area:units = "m2" ;

// global attributes:
                :title = "VIC domain data" ;
                :Conventions = "CF-1.6" ;
                :history = "Wed Oct 12 15:48:42 2016: ncap2 -s mask=int(mask) ArkRed.domain.nc.float_mask ArkRed.domain.nc\n",
                        "created by ymao, 2016-09-23 18:17:58.761256" ;
                :user_comment = "VIC domain data" ;
                :source = "generated from VIC CONUS 1.8 deg model parameters, see Maurer et al. (2002) for more information" ;
                :nco_openmp_thread_number = 1 ;
}
```

