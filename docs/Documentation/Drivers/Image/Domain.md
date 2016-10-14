# VIC Domain file

The Image Driver uses the [NetCDF](http://www.unidata.ucar.edu/software/netcdf/) file format to define model running domain.

Below is a list of variables in the domain netCDF file. The dimensions of the netCDF file are `lat` and `lon`. Note that here only the type of variable (i.e., MASK, AREA, FRAC, LAT and LON) are listed; corresponding variable names in the input netCDF file is specified by user in the [Global Parameter File](GlobalParam.md). All the listed variables are required.

| Variable   | Dimension   | Units    | Type   | Description |
|------------|-------------|----------|--------|-------------|
| LAT        | [lat]       | degree   | double | Latitudes   |
| LON        | [lon]       | degree   | double | Longitues   |
| MASK       | [lat, lon]  | N/A      | integer | Mask of VIC run. 1 for active cells for VIC run; 0 indicates inactive grid cells. VIC will not run at grid cells with MASK = 0 or missing. |
| AREA       | [lat, lon]  | m2       | double | Area of grid cell.   |
| FRAC       | [lat, lon]  | N/A      | double | Fraction of grid cell that is land. |

