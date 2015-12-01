# VIC 4.2.a.irrigation

- This is a VIC support branch developed by Ingjerd Haddeland and Ted Bohn to  model irrigation requirements ("potential irrigation"). Actual irrigation, i.e. with water limitation included, requires use of the routing (reservoir) model with dams included. This version was developed in Jan 2015 and was used by Tian Zhou to run simulations for the [ISI-MIP](https://www.pik-potsdam.de/research/climate-impacts-and-vulnerabilities/research/rd2-cross-cutting-activities/isi-mip) 2.1 global water sector .

### File changes made from VIC 4.2:
- full_energy.c: This is where irrigation water is added to crops if needed. Either potential irrigation (i.e. not limited by water availability) or actual irrigation where water availability is taken into account. If water availability is considered, VIC first looks for water originating from the same cell as the one in question (runoff and/or baseflow already produced in cell ), from the river network flowing through the cell in question (routing model must interact with VIC), or from more distant upstream reservoir(s) (also requires interaction with the routing model). Water from the river is input as forcing variable `irr_run`, and water from reservoirs is input as forcing variable `irr_with`.  
The subroutine also takes into account possible varying crop fraction area, i.e. the area that is to be irrigated may vary in time.
Irrigation does not happen when snow is present or when the temperature is below 7 degrees Celsius, even if your vegetation library defines month in question as a month when irrigation may happen.
- vicNl_def.h: Some additional irrigation parameters/constants/variables defined.
- make_in_and_outputfiles.c: Reading of `crop_frac` files enabled
- initialize_global.c: Global parameters, input forcings
  `IRR_FREE FALSE` is default, i.e. no irrigation
- read_forcing_data.c: Reading of `crop_frac` files enabled
- get_global_param.c: Reading of `crop_frac` files enabled
- get_force_type.c: Reading of `irr_with` and `irr_run` forcings enabled (not applicable without the routing scheme interacting with VIC)
- close_files.c: Close additional forcing files.
- display_current_settings.c: Minor changes
- output_list_utils.c: Additional output variables enabled
- initialize_soil.c: Added `cell[veg][band].irr_extract=0`
- put_data.c: Output irrigation terms calculated
- alloc_atmos.c: Allocate/free memory for `atmos->irr_run` and `irr_with`
- initialize_veg.c: Initializing `crop_frac`, `irrig` and `irr_apply` for every `snow_band`
- global.h: `double ref_crop_frac = { 0.0, 0.0, 1.0, 1.0 }` added
- read_vegparam.c: Minor  changes related to reading of `crop_frac` info.
- surface_fluxes.c: `ppt += veg_var->irrig` (irrig water added to ppt)
- read_veglib.c: Reading and initializing irrigation parameters/variables from veg_lib file.
- alloc_veg_hist.c: Allocate and free memory, `veg_hist[][].crop_frac`
- initialize_atmos.c: Initialize `irr_run` and `irr_with`. Assign crop area fraction from veg_lib or from (sub)daily crop_frac files if provided, normalize veg. Read information on available water for irrigation from external (forcing) files if provided.
- vicNl.h: Argument `all_vars_struct *` added
- initialize_model_state.c: Argument `all_vars_struct` in list of arguments. Crop structures initialized
- vicNl.c: Crop variables added

### Input files differences compared to regular VIC inputs

- The sample input and output files are located in **samples.VIC.4.2.a.irrigation** directory.

- Global file: Two extra lines are included.
`IRRIGATION      FALSE/TRUE`
`IRR_FREE        FALSE/TRUE`
`IRRIGATION` indicates whether you want irrigation included (`TRUE`); `IRR_FREE` indicates whether the water for irrigation is freely available (`TRUE`) or not (`FALSE`). If `IRR_FREE` set to `FALSE`, two extra columns are needed in forcing files (`IRR_RUN` and `IRR_WITH`) describing the runoff and water availability in the current grid cell.

- Vegetation parameter file: An example of a cell in the veg parameter file with irrigated vegetation. There is one extra flag (0 or 1) on the first line of each tile. In this case, the last tile (veg number 79, defined in the veg library file) has the flag set to 1, which means this tile will be irrigated and the varying crop fraction will be read from crop_frac forcing file. Note: In order to clearly see the effects of the irrigation, the fraction of the crop tile (0.3229) is made 10 times larger than the real value in this example.   

```
38597 6
7 0.1094	 0.30  0.60  0.70  0.40  0
0.5870  1.0750  1.3500  1.5870  1.7870  1.7120  1.7880  2.1000  2.0370  1.5500  1.0250  0.6620
8 0.2149	 0.30  0.70  0.70  0.30  0
0.2620  0.2870  0.4000  0.4750  0.3000  0.2500  0.4870  0.9250  0.5620  0.2750  0.2250  0.2620
9 0.2242	 0.30  0.70  0.70  0.30  0
0.5500  0.6120  0.7620  0.4500  0.2250  0.3120  0.7250  0.9380  0.4250  0.2000  0.2000  0.3630
10 0.1233	 0.30  0.80  0.70  0.20  0
0.2370  0.2750  0.3880  0.5870  0.8500  1.0500  0.5870  0.4130  0.3130  0.2870  0.2500  0.2250
11 0.0003	 0.30  0.50  0.70  0.50  0
0.4300  0.8080  1.3730  1.7830  1.2110  0.3660  0.6290  0.9010  0.5010  0.2400  0.2250  0.2960
79 0.3229	 0.30  0.50  0.70  0.50  1
1.4460  4.2880  4.5380  2.4280  1.2810  0.3810  1.2810  1.2810  1.2810  1.2810  0.1330  0.6290
```

- Crop fraction file: A crop fraction file is needed as an additional forcing file when the crop fraction changes over time. Similar to the climate forcing file, the crop fraction file is one file per grid cell, describing the crop area fraction (0-1) as it changes over time in the crop tile for each grid cell. Note that the fraction only represents the crop fraction in the crop tile, not in the entire grid cell.

- Vegetation library: There are 12 more columns (JAN-IRR to DEC-IRR) set to 0 or 1 for each veg type, indicating if this veg type will be irrigated when `IRRIGATION` is set to `TRUE` in global file.

- Soil file: Same as regular soil file.
