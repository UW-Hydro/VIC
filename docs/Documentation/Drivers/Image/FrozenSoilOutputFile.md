# Frozen Soil Output File - Image Driver

When the model is run with the frozen soil algorithms, a third output file is produced which contains soil thermal output parameters. Since the frozen soils algorithm only works in with the full energy balance, there is only one format for the output file.

The default flux output file is in netCDF format, with the same data structure described in the [image driver output format](OutputFormatting.md). the list of default flux variables are:

| Variable   Name     | Dimension                | Units    | Description                                                 |
|---------------------|--------------------------|----------|-------------------------------------------------------------|
| OUT_FDEPTH          | [time, front, lat, lon]  | cm       | Depth of   the freezing front                               |
| OUT_TDEPTH          | [time, front, lat, lon]  | cm       | Depth of   the thawing front                                |
| OUT_SOIL_MOIST      | [time, nlayer, lat, lon] | mm       | Total   soil moisture in each layer (liquid water plus ice) |
| OUT_SOIL_MOIST_FRAC | [time, lat, lon]         | fraction | Fraction   of total soil moisture in each layer             |
