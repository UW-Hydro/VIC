FORCE_TYPE defines which forcing variables are in the forcing file, and the order in which they appear. These lines must follow either FORCING1 or FORCING2, and N_TYPES so that the model knows which file to associate the columns with and how any columns are present. If the forcing file is in binary format, FORCE_TYPE must also define whether each value is a SIGNED or UNSIGNED short int and the factor by which it needs to be multiplied.

Possible forcing file data types are:

* AIR_TEMP - sub-daily air temperature (C)
* ALBEDO - surface albedo (fraction)
* DENSITY - atmospheric density (kg/m^3)
* PREC - precipitation (mm)
* PRESSURE - atmospheric pressure (kPa)
* SHORTWAVE - shortwave radiation (W/m<sup>2</sup>)
* TMAX - daily maximum temperature (C)
* TMIN - daily minimum temperature (C)
* VP - atmospheric vapor pressure (kPa)
* WIND - wind speed (m/s)
* SKIP - used to indicate a data column which is not read into the model</menu>

_Examples._ a standard four column daily forcing data file will be defined as:

## ASCII File
```
FORCING1  FORCING_DATA/LDAS_ONE_DEGREE/data_
N_TYPES    4
FORCE_TYPE  PREC
FORCE_TYPE  TMAX
FORCE_TYPE  TMIN
FORCE_TYPE  WIND
FORCE_FORMAT  ASCII
FORCE_DT  24
```

## Binary File
```
FORCING1  FORCING_DATA/LDAS_ONE_DEGREE/data_
N_TYPES    4
FORCE_TYPE  PREC  UNSIGNED  40
FORCE_TYPE  TMAX  SIGNED    100
FORCE_TYPE  TMIN  SIGNED    100
FORCE_TYPE  WIND  SIGNED    100
FORCE_FORMAT  BINARY
FORCE_ENDIAN  LITTLE
FORCE_DT  24
```
