# VIC Model State File

VIC can save the hydrologic state from any point in the simulation (usually the final state) to a file for the purpose of re-starting the simulation later (as an initial state file). This is useful for simulations that require lengthy spin-up periods or ensemble methods. The initial state file is not required; if it is not specified, VIC will use a default starting condition, and will take initial soil moisture contents from the values specified in the [Soil Parameter File](SoilParam.md). An Initial State File may be prepared simply by running VIC with the necessary state file options in the [global parameter file](GlobalParam.md#DefineStateFiles).

The model state file contains all information needed by the VIC model to "warm"-start a simulation (i.e. start from "realistic" conditions, or re-start a simulation exactly where the model previously stopped). To read an initial state file, or to save a "final" state file, the appropriate options should be set in the [global parameter file](GlobalParam.md#DefineStateFiles).

The timestamp of the state file represents the instantaneous time for which the values of the state variables are valid. This corresponds to the end of the time interval after which the state file is written out and the beginning of the time interval for which the model is started. For example, if the MODEL_STEPS_PER_DAY is 24 (hourly) and the last time step for which the model is run is 1999-09-20 23:00:00, the the state file will be stamped 1999-09-21 00:00:00. This state file can then be used to restart a model run whose starting time will be 1999-09-21 00:00:00.

The state file has two header lines used by the VIC model to verify that the model is set-up correctly. Following the header there are repeating blocks of lines which define all variables for each grid cell. Each block starts with a line defining variables held constant across the grid cell and is then followed by lines for all vegetation types (including bare soil) and all snow bands [ `(number of vegetation types) * (number of snow bands) = (number of lines)` ].

*   [File Header](#FileHeader)
*   [Grid Cell Information](#GridInfo)
*   [Vegetation and Snow Band Information](#VSBInfo)
*   [Lake Information](#LakeInfo)</menu>

* * *

## File Header

The file header appears once at the top of the state file and is used by the VIC model to verify the parameters initialized with the model control file match those for the data stored in the state file. An error is reported if the model state file cannot be used to initialize the current state.

Define date of state file, always saved at hour 0 of the defined day.

| Column    | Name          | Type  | Description                       |
|--------   |------------   |------ |-------------------------------    |
| 1         | STATEYEAR     | int   | Year of the model state file      |
| 2         | STATEMONTH    | int   | Month of the model state file     |
| 3         | STATEDAY      | int   | Day of the model state file       |

Define global soil limits.

| Column    | Name      | Type  | Description                       |
|--------   |--------   |------ |--------------------------------   |
| 1         | Nlayer    | int   | Number of soil moisture layers    |
| 2         | Nnodes    | int   | Number of soil thermal nodes      |

* * *

## Grid Cell Information

This line is used to define the parameters that are constant for the entire grid cell. It first appears after the file header lines and must be defined for each grid cell.

| Column                    | Name          | Type      | Description                                                                                                               |
|-------------------------  |------------   |--------   |-------------------------------------------------------------------------------------------------------------------------  |
| 1                         | cellnum       | int       | Current cell number                                                                                                       |
| 2                         | Nveg          | int       | Number of vegetation types in the grid cell                                                                               |
| 3                         | Nbands        | int       | Number of snow elevation bands in the grid cell                                                                           |
| 4:(3+Nnodes)              | dz_node       | double    | Distances between soil thermal nodes [m]                                                                                  |
| (4+Nnodes):(3+2\*Nnodes)  | node_depth    | double    | Depth from surface of each soil thermal node (first node should have a depth of 0m indicating it is at the surface) [m]   |

* * *

## Vegetation and Snow Band Information

Each grid cell information line must be followed with lines for all defined vegetation types (plus bare soil) and snow bands ( `[Nveg+1]*Nbands = number of lines required` ). These lines contain information about the storage of moisture within each fractional coverage type. If the model is being run with distributed precipitation, the wet and dry fractions are averaged before the model state is stored and the model is always initialized with a mu value of 1.

| Column                                                                        | Name              | Type      | Description                                                                           |
|------------------------------------------------------------------------------ |---------------    |--------   |-------------------------------------------------------------------------------------  |
| 1                                                                             | veg               | int       | Current vegetation number as defined in the vegetation parameter file                 |
| 2                                                                             | band              | int       | Current snow band number                                                              |
| 3:(2+Nlayer)                                                                  | moist             | double    | Soil layer total moisture contents including ice [mm]                                 |
| (3+Nlayer):(2+2\*Nlayer)                                                      | ice               | double    | Soil layer ice contents [mm]                                                          |
| (3+2\*Nlayer)                                                                 | Wdew              | double    | Amount of dew stored in the vegetation [mm]. Not defined for bare soil (veg = Nveg)   |
| (The next 5 terms appear only if CARBON=TRUE)                                 |                   |           |                                                                                       |
| (4+2\*Nlayer)                                                                 | AnnualNPP         | double    | running total annual NPP [gC/m<sup>2</sup>]                                           |
| (5+2\*Nlayer)                                                                 | AnnualNPPPrev     | double    | total annual NPP from previous year [gC/m<sup>2</sup>]                                |
| (6+2\*Nlayer)                                                                 | CLitter           | double    | carbon storage in litter pool [gC/m<sup>2</sup>]                                      |
| (7+2\*Nlayer)                                                                 | CInter            | double    | carbon storage in intermediate pool [gC/m<sup>2</sup>]                                |
| (8+2\*Nlayer)                                                                 | CSlow             | double    | carbon storage in slow pool [gC/m<sup>2</sup>]                                        |
| (The following terms always appear; if CARBON is TRUE add 5 to column index)  |                   |           |                                                                                       |
| (4+2\*Nlayer)                                                                 | last_snow         | int       | Number of model time steps since the last new snow                                    |
| (5+2\*Nlayer)                                                                 | MELTING           | char      | Flag to indicate whether snowpack is in accumulation or melting phase                 |
| (6+2\*Nlayer)                                                                 | coverage          | double    | Fraction of grid cell area covered by snow                                            |
| (7+2\*Nlayer)                                                                 | swq               | double    | Snow water equivalence of the pack [m]                                                |
| (8+2\*Nlayer)                                                                 | surf_temp         | double    | Temperature of the surface layer in the snow algorithm [C]                            |
| (9+2\*Nlayer)                                                                 | surf_water        | double    | Liquid water content of the surface layer [m]                                         |
| (10+2\*Nlayer)                                                                | pack_temp         | double    | Temperature of the pack layer in the snow algorithm [C]                               |
| (11+2\*Nlayer)                                                                | pack_water        | double    | Liquid water content of the pack layer [m]                                         |
| (12+2\*Nlayer)                                                                | density           | double    | Density of the snowpack [kg/m<sup>3</sup>]                                            |
| (13+2\*Nlayer)                                                                | coldcontent       | double    | Cold content of snow pack [J/m<sup>2</sup>]                                           |
| (14+2\*Nlayer)                                                                | snow_canopy       | double    | Snow water equivalence stored in the canopy [m]                                       |
| (15+2\*Nlayer):(15+2\*Nlayer+Nnodes)                                          | node_T            | double    | Soil temperature at each of the defined soil thermal nodes [C]                        |

* * *

## Lake Information (only when LAKES are turned on in the [global parameter file](GlobalParam.md#DefineStateFiles))

| Column                                                                        | Name          | Type      | Description                                                                           |
|------------------------------------------------------------------------------ |-------------- |--------   |------------------------------------------------------------------------------------   |
| 1:(0+Nlayer)                                                                  | moist         | double    | Soil layer total moisture contents including ice [mm]                                 |
| (1+Nlayer):(0+2\*Nlayer)                                                      | ice           | double    | Soil layer ice contents [mm]                                                          |
| (The next 3 terms appear only if CARBON=TRUE)                                 |               |           |                                                                                       |
| (1+2\*Nlayer)                                                                 | CLitter       | double    | carbon storage in litter pool [gC/m2]                                                 |
| (2+2\*Nlayer)                                                                 | CInter        | double    | carbon storage in intermediate pool [gC/m2]                                           |
| (3+2\*Nlayer)                                                                 | CSlow         | double    | carbon storage in slow pool [gC/m2]                                                   |
| (The following terms always appear; if CARBON is TRUE add 3 to column index)  |               |           |                                                                                       |
| (1+2\*Nlayer)                                                                 | last_snow     | int       | Number of model time steps since the last new snow                                    |
| (2+2\*Nlayer)                                                                 | MELTING       | char      | Flag to indicate whether snowpack is in accumulation or melting phase                 |
| (3+2\*Nlayer)                                                                 | coverage      | double    | Fraction of grid cell area covered by snow                                            |
| (4+2\*Nlayer)                                                                 | swq           | double    | Snow water equivalence of the pack [m]                                                |
| (5+2\*Nlayer)                                                                 | surf_temp     | double    | Temperature of the surface layer in the snow algorithm [C]                            |
| (6+2\*Nlayer)                                                                 | surf_water    | double    | Liquid water content of the surface layer [m]                                         |
| (7+2\*Nlayer)                                                                 | pack_temp     | double    | Temperature of the pack layer in the snow algorithm [C]                               |
| (8+2\*Nlayer)                                                                 | pack_water    | double    | Liquid water content of the surface layer [m]                                         |
| (9+2\*Nlayer)                                                                 | density       | double    | Density of the snowpack [kg/m<sup>3</sup>]                                            |
| (10+2\*Nlayer)                                                                | coldcontent   | double    | Cold content of snow pack [J/m<sup>2</sup>]                                           |
| (11+2\*Nlayer)                                                                | snow_canopy   | double    | Snow water equivalence stored in the canopy [m]                                       |
| (12+2\*Nlayer):(12+2\*Nlayer+Nnodes)                                          | node_T        | double    | Soil temperature at each of the defined soil thermal nodes [C]                        |
| (13+2\*Nlayer+Nnodes)                                                         | activenod     | int       | number of active lake nodes                                                           |
| (14+2\*Nlayer+Nnodes)                                                         | dz            | double    | Vertical thickness of all horizontal water layers below the surface layer (m)         |
| (15+2\*Nlayer+Nnodes)                                                         | surfdz        | double    | Vertical thickness of surface (top) water layer (m)                                   |
| (16+2\*Nlayer+Nnodes)                                                         | ldepth        | double    | Current depth of liquid water in lake (distance from surface to deepest point) (m)    |
| (17+2\*Nlayer+Nnodes):(16+2\*Nlayer+Nnodes+numnod)                            | surface       | double    | Area of horizontal cross-section of liquid water in lake at each node (m<sup>2</sup>) |
| (17+2\*Nlayer+Nnodes+numnod)                                                  | sarea         | double    | Current surface area of ice+liquid water on lake surface (m<sup>2</sup>)              |
| (18+2\*Nlayer+Nnodes+numnod)                                                  | volume        | double    | Current lake water volume, including liquid water equivalent of lake ice (m<sup>3</sup>) |
| (19+2\*Nlayer+Nnodes+numnod):(18+2\*Nlayer+Nnodes+2\*numnod)                  | temp          | double    | Lake water temperature at each node (C)                                               |
| (19+2\*Nlayer+Nnodes+2\*numnod)                                               | tempavg       | double    | Average liquid water temperature of entire lake (C)                                   |
| (20+2\*Nlayer+Nnodes+2\*numnod)                                               | areai         | double    | Area of ice coverage (at beginning of time step) (m<sup>2</sup>)                      |
| (21+2\*Nlayer+Nnodes+2\*numnod)                                               | new_ice_area  | double    | Area of ice coverage (at end of time step) (m<sup>2</sup>)                            |
| (22+2\*Nlayer+Nnodes+2\*numnod)                                               | ice_water_eq  | double    | Liquid water equivalent volume of lake ice (m<sup>3</sup>)                            |
| (23+2\*Nlayer+Nnodes+2\*numnod)                                               | hice          | double    | Height of lake ice at thickest point (m)                                              |
| (24+2\*Nlayer+Nnodes+2\*numnod)                                               | tempi         | double    | Lake ice temperature (C)                                                              |
| (25+2\*Nlayer+Nnodes+2\*numnod)                                               | swe           | double    | Water equivalence of lake snow cover - end of step (m<sup>3</sup>)                    |
| (26+2\*Nlayer+Nnodes+2\*numnod)                                               | surf_temp     | double    | Temperature of surface snow layer (C)                                                 |
| (27+2\*Nlayer+Nnodes+2\*numnod)                                               | pack_temp     | double    | Temperature of pack snow layer (C)                                                    |
| (28+2\*Nlayer+Nnodes+2\*numnod)                                               | coldcontent   | double    | cold content of snow pack                                                             |
| (29+2\*Nlayer+Nnodes+2\*numnod)                                               | surf_water    | double    | Water content of surface snow layer (m)                                   |
| (30+2\*Nlayer+Nnodes+2\*numnod)                                               | pack_water    | double    | Water content of pack snow layer (m)                                      |
| (31+2\*Nlayer+Nnodes+2\*numnod)                                               | SAlbedo       | double    | Albedo of lake snow (fraction)                                                        |
| (32+2\*Nlayer+Nnodes+2\*numnod)                                               | sdepth        | double    | Depth of snow on top of ice (m)                                           |

## State File Example

From the Stehekin basin, using 3 soil layers and 10 thermal nodes. Note: indented text indicates the continuation of the previous line. Only the first grid cell is included.

    1948 12 31
    3 10
    86340 5 5 0.100000  0.100000  0.100000  0.576923  0.576923  0.576923  0.576923  0.576923  0.576923  0.576923  0.000000  0.100000  0.200000  0.538462  1.115385  1.692308  2.269231  2.846154  3.423077  4.000000
    1.000000 88 -999
    0 0 17.061740 56.710901 154.076105 0.000000 0.000000 0.000000 0.000000 49 0 1.000000 0.282294 -26.007470 0.000000 -10.385944 0.000000 287.137869 -6826960.769076 0.006871 -22.711803 -0.876051 -0.763897 -0.457698 -0.124445  0.062717 0.167832 0.226867 0.260023 0.302500
    0 1 16.681905 55.980988 293.638277 0.000000 0.000000 0.000000 0.000000 49 0 1.000000 0.415990 -27.784704 0.000000 -11.030049 0.000000 293.191662 -7293484.811361 0.009397 -25.822053 -0.876051 -0.763897 -0.457698 -0.124445  0.062717 0.167832 0.226867 0.260023 0.302500
    0 2 16.610667 55.843397 348.196273 0.000000 0.000000 0.000000 0.000000 49 0 1.000000 0.446568 -28.489163 0.000000 -11.626803 0.000000 292.111407 -7478405.215830 0.009737 -27.036902 -0.876051 -0.763897 -0.457698 -0.124445  0.062717 0.167832 0.226867 0.260023 0.302500
    0 3 16.613013 55.867093 377.347556 0.000000 0.000000 0.000000 0.000000 49 0 1.000000 0.454764 -29.028473 0.000000 -12.121350 0.000000 289.896465 -7619974.228963 0.009780 -27.959252 -0.876051 -0.763897 -0.457698 -0.124445  0.062717 0.167832 0.226867 0.260023 0.302500
    0 4 16.607877 55.858645 395.803394 0.000000 0.000000 0.000000 0.000000 49 0 1.000000 0.468959 -29.680128 0.000000 -12.676810 0.000000 288.316397 -7791033.525378 0.009791 -29.007702 -0.876051 -0.763897 -0.457698 -0.124445  0.062717 0.167832 0.226867 0.260023 0.302500
    1.000000 -56 -999
    1 0 17.267470 57.916699 140.952406 0.000000 0.000000 0.000000 0.000000 24 0 1.000000 0.267034 -21.793677 0.000000 -10.021876 0.000000 286.635428 -5720840.159419 0.000000 -22.711803 -0.876051 -0.763897 -0.457698 -0.124445  0.062717 0.167832 0.226867 0.260023 0.302500
    1 1 16.687072 56.058484 222.699231 0.000000 0.000000 0.000000 0.000000 24 0 1.000000 0.414656 -24.283520 0.000000 -10.497186 0.000000 296.026782 -6374423.878210 0.000000 -25.822053 -0.876051 -0.763897 -0.457698 -0.124445  0.062717 0.167832 0.226867 0.260023 0.302500
    1 2 16.657359 56.045383 283.580321 0.000000 0.000000 0.000000 0.000000 24 0 1.000000 0.447827 -25.253470 0.000000 -11.192870 0.000000 295.065053 -6629036.003708 0.000000 -27.036902 -0.876051 -0.763897 -0.457698 -0.124445  0.062717 0.167832 0.226867 0.260023 0.302500
    1 3 16.640994 55.986534 317.416672 0.000000 0.000000 0.000000 0.000000 24 0 1.000000 0.459124 -26.010046 0.000000 -11.661158 0.000000 293.237947 -6827637.187789 0.000000 -27.959252 -0.876051 -0.763897 -0.457698 -0.124445  0.062717 0.167832 0.226867 0.260023 0.302500
    1 4 16.615881 55.901837 358.898564 0.000000 0.000000 0.000000 0.000000 24 0 1.000000 0.475501 -26.857030 0.000000 -12.244615 0.000000 291.835035 -7049970.387308 0.000000 -29.007702 -0.876051 -0.763897 -0.457698 -0.124445  0.062717 0.167832 0.226867 0.260023 0.302500
    1.000000 62 -999
    2 0 16.779661 55.074733 137.923388 0.000000 0.000000 0.000000 0.000000 49 0 1.000000 0.259040 -30.180932 0.000000 -10.236263 0.000000 279.660711 -7922494.701340 0.000000 -22.711803 -0.876051 -0.763897 -0.457698 -0.124445  0.062717 0.167832 0.226867 0.260023 0.302500
    2 1 16.631105 55.800206 210.712879 0.000000 0.000000 0.000000 0.000000 49 0 1.000000 0.400361 -31.141543 0.000000 -10.549196 0.000000 289.567469 -8174654.958190 0.000000 -25.822053 -0.876051 -0.763897 -0.457698 -0.124445  0.062717 0.167832 0.226867 0.260023 0.302500
    2 2 16.620325 55.890479 244.646706 0.000000 0.000000 0.000000 0.000000 49 0 1.000000 0.439078 -31.621836 0.000000 -11.055820 0.000000 290.450159 -8300732.030048 0.000000 -27.036902 -0.876051 -0.763897 -0.457698 -0.124445  0.062717 0.167832 0.226867 0.260023 0.302500
    2 3 16.616265 55.908883 263.550104 0.000000 0.000000 0.000000 0.000000 49 0 1.000000 0.453719 -32.009960 0.000000 -11.547218 0.000000 289.714849 -8402614.488123 0.000000 -27.959252 -0.876051 -0.763897 -0.457698 -0.124445  0.062717 0.167832 0.226867 0.260023 0.302500
    2 4 16.617371 55.931029 287.245168 0.000000 0.000000 0.000000 0.000000 49 0 1.000000 0.472990 -32.466258 0.000000 -12.127620 0.000000 289.322428 -8522392.634618 0.000000 -29.007702 -0.876051 -0.763897 -0.457698 -0.124445  0.062717 0.167832 0.226867 0.260023 0.302500
    1.000000 0 -999
    3 0 16.781958 56.452564 171.547740 0.000000 0.000000 0.000000 0.000000 49 0 1.000000 0.242929 -29.987999 0.000000 -10.293697 0.000000 275.675680 -7871849.840076 0.000000 -22.711803 -0.876051 -0.763897 -0.457698 -0.124445  0.062717 0.167832 0.226867 0.260023 0.302500
    3 1 16.631348 55.971010 247.826283 0.000000 0.000000 0.000000 0.000000 49 0 1.000000 0.392567 -31.009633 0.000000 -10.519804 0.000000 288.836171 -8140028.766553 0.000000 -25.822053 -0.876051 -0.763897 -0.457698 -0.124445  0.062717 0.167832 0.226867 0.260023 0.302500
    3 2 16.620545 55.932068 279.497319 0.000000 0.000000 0.000000 0.000000 49 0 1.000000 0.433544 -31.509586 0.000000 -10.985153 0.000000 290.232655 -8271266.356764 0.000000 -27.036902 -0.876051 -0.763897 -0.457698 -0.124445  0.062717 0.167832 0.226867 0.260023 0.302500
    3 3 16.616829 55.925264 296.552013 0.000000 0.000000 0.000000 0.000000 49 0 1.000000 0.450121 -31.908823 0.000000 -11.465947 0.000000 289.775964 -8376066.045837 0.000000 -27.959252 -0.876051 -0.763897 -0.457698 -0.124445  0.062717 0.167832 0.226867 0.260023 0.302500
    3 4 16.617693 55.939020 318.298371 0.000000 0.000000 0.000000 0.000000 49 0 1.000000 0.471458 -32.372915 0.000000 -12.049658 0.000000 289.642090 -8497890.284478 0.000000 -29.007702 -0.876051 -0.763897 -0.457698 -0.124445  0.062717 0.167832 0.226867 0.260023 0.302500
    1.000000 88 -999
    4 0 16.785370 55.829994 123.658906 0.000000 0.000000 0.000000 0.000000 49 0 1.000000 0.227140 -29.809920 0.000000 -10.454044 0.000000 271.511679 -7825104.061686 0.000000 -22.711803 -0.876051 -0.763897 -0.457698 -0.124445  0.062717 0.167832 0.226867 0.260023 0.302500
    4 1 16.629006 55.695153 189.592901 0.000000 0.000000 0.000000 0.000000 49 0 1.000000 0.384277 -30.886970 0.000000 -10.557114 0.000000 287.944539 -8107829.508228 0.000000 -25.822053 -0.876051 -0.763897 -0.457698 -0.124445  0.062717 0.167832 0.226867 0.260023 0.302500
    4 2 16.616612 55.824128 225.926605 0.000000 0.000000 0.000000 0.000000 49 0 1.000000 0.427881 -31.405603 0.000000 -10.956261 0.000000 289.827263 -8243970.868187 0.000000 -27.036902 -0.876051 -0.763897 -0.457698 -0.124445  0.062717 0.167832 0.226867 0.260023 0.302500
    4 3 16.614931 55.881099 246.627049 0.000000 0.000000 0.000000 0.000000 49 0 1.000000 0.446246 -31.814116 0.000000 -11.418063 0.000000 289.633796 -8351205.485392 0.000000 -27.959252 -0.876051 -0.763897 -0.457698 -0.124445  0.062717 0.167832 0.226867 0.260023 0.302500
    4 4 16.616750 55.921897 271.619558 0.000000 0.000000 0.000000 0.000000 49 0 1.000000 0.469795 -32.285656 0.000000 -12.000788 0.000000 289.840151 -8474984.674525 0.000000 -29.007702 -0.876051 -0.763897 -0.457698 -0.124445  0.062717 0.167832 0.226867 0.260023 0.302500
    1.000000 -56 -999
    5 0 15.442000 46.327000 154.425000 0.000000 0.000000 0.000000 0 0 0.000000  0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000   0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000    0.000000 0.000000
    5 1 15.442000 46.327000 154.425000 0.000000 0.000000 0.000000 0 0 0.000000  0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000   0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000    0.000000 0.000000
    5 2 15.442000 46.327000 154.425000 0.000000 0.000000 0.000000 0 0 0.000000  0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000   0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000    0.000000 0.000000
    5 3 15.442000 46.327000 154.425000 0.000000 0.000000 0.000000 0 0 0.000000  0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000   0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000    0.000000 0.000000
    5 4 15.442000 46.327000 154.425000 0.000000 0.000000 0.000000 0 0 0.000000  0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000   0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000    0.000000 0.000000
