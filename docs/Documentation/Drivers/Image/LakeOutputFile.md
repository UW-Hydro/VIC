# Lake Output File

The lake file contains information about the lake fraction of the grid cell. Lake files are only written when LAKES is set equal to a valid lake parameter file, in the [global parameter file](GlobalParam.md). Lake file name is `lake.**.nc`, where `**` is determined by AGGFREQ variable specified in the [global parameter file](GlobalParam.md).

The default flux output file is in netCDF format, with the same data structure described in the [image driver output format](OutputFormatting.md). the list of default flux variables are:

| Variable   Name     | Units    | Description                                               |
|---------------------|----------|-----------------------------------------------------------|
| OUT_LAKE_ICE_TEMP   | C        | Temperature   of lake ice                                 |
| OUT_LAKE_ICE_HEIGHT | cm       | Thickness   of lake ice                                   |
| OUT_LAKE_ICE_FRACT  | fraction | Fractional   coverage of lake ice                         |
| OUT_LAKE_DEPTH      | m        | Lake depth   (distance between surface and deepest point) |
| OUT_LAKE_SURF_AREA  | m2       | Lake   surface area                                       |
| OUT_LAKE_VOLUME     | m3       | Lake   volume                                             |
| OUT_LAKE_SURF_TEMP  | C        | Lake   surface temperature                                |
| OUT_EVAP_LAKE       | mm       | Net   evaporation from lake surface                       |
