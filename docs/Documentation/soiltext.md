# VIC Model - Soil Texture

## Sample Index of Soil Hydraulic Properties\*

| USDA Class    | Soil Type     | % Sand    | % Clay    | Bulk Density  | Field Capacity    | Wilting Point     | Porosity  | Saturated Hydraulic Conductivity  | Slope of Retention Curve (in log space) `b` |
|------------   |-----------    |--------   |--------   |-------------- |----------------   |---------------    |---------- |---------------------------------- |---------------------------------------------  |
|               |               |           |           | g/cm<sup>3</sup>| cm<sup>3</sup></sup>  | cm<sup>3</sup></sup>  | fraction  | cm/hr                   |                                               |
| 1             | s             | 94.83     | 2.27      | 1.49          | 0.08              | 0.03              | 0.43      | 38.41                             | 4.1                                           |
| 2             | ls            | 85.23     | 6.53      | 1.52          | 0.15              | 0.06              | 0.42      | 10.87                             | 3.99                                          |
| 3             | sl            | 69.28     | 12.48     | 1.57          | 0.21              | 0.09              | 0.4       | 5.24                              | 4.84                                          |
| 4             | sil           | 19.28     | 17.11     | 1.42          | 0.32              | 0.12              | 0.46      | 3.96                              | 3.79                                          |
| 5             | si            | 4.5       | 8.3       | 1.28          | 0.28              | 0.08              | 0.52      | 8.59                              | 3.05                                          |
| 6             | l             | 41        | 20.69     | 1.49          | 0.29              | 0.14              | 0.43      | 1.97                              | 5.3                                           |
| 7             | scl           | 60.97     | 26.33     | 1.6           | 0.27              | 0.17              | 0.39      | 2.4                               | 8.66                                          |
| 8             | sicl          | 9.04      | 33.05     | 1.38          | 0.36              | 0.21              | 0.48      | 4.57                              | 7.48                                          |
| 9             | cl            | 30.08     | 33.46     | 1.43          | 0.34              | 0.21              | 0.46      | 1.77                              | 8.02                                          |
| 10            | sc            | 50.32     | 39.3      | 1.57          | 0.31              | 0.23              | 0.41      | 1.19                              | 13                                            |
| 11            | sic           | 8.18      | 44.58     | 1.35          | 0.37              | 0.25              | 0.49      | 2.95                              | 9.76                                          |
| 12            | c             | 24.71     | 52.46     | 1.39          | 0.36              | 0.27              | 0.47      | 3.18                              | 12.28                                         |

\* Source is "Average hydraulic properties of ARS soil texture classes," draft dated February, 2000 by J. Schaake. This expanded the work of others and included a total of 2128 soil samples. Wilting point is the fractional water content at 15 bar tension; field capacity is the fractional water content at 1/3 bar tension.

`b` is as used in Campbell's equation. ref: Cosby et al., A Statistical exploration of the relationships of soil moisture characteristics to the physical properties of soils, _Water Resources Research_ 20(6): 682-690, 1984.

**Note that units in this table differ from those used in the VIC model. In particular:**

1.  Bulk Density for the VIC model should be in kg/m<sup>3</sup> (g/cm<sup>3</sup> * 1000 = kg/m<sup>3</sup>
2.  Field Capacity and Wilting Point for the VIC model are described as a fraction of the maximum moisture, where the maximum moisture for each soil layer is the depth times the porosity. This can be obtained from the above data in cm<sup>3</sup>/cm<sup>3</sup> by dividing the values by the fractional porosity in the table.
3.  Saturated Hydraulic Conductivity for the VIC model is in mm/day (cm/hr * 240 = mm/day)
4.  By default, VIC uses the Brooks-Corey relationship for unsaturated flow. The Brooks-Corey exponent, n, can be estimated from the slope of the Retention Curve (in log space), b, in the table above by: `n = 3 + 2b`
