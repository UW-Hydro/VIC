# VIC 4.2.a.irrigation

- Developed by Ingjerd Haddeland and Ted Bohn. This stand alone version only possible to model irrigation requirements ("potential irrigation"). Actual irrigation, i.e. with water limitation included, requires use of the routing (reservoir) model with dams included. This version was developed in Jan 2015 and was used to simulate ISI-MIP 2.1 global water sector by Tian Zhou.

### Input files differences compared to regular VIC inputs

- The sample input and output files are located in **samples.VIC.4.2.a.irrigation** directory.

- Global file: Two extra lines are included.
IRRIGATION      FALSE/TRUE
IRR_FREE        FALSE/TRUE
IRRIGATION means if you want irrigation included; IRR_FREE means if the water for irrigation is freely available or not. If IRR_FREE set to FALSE, two extra columns are needed in forcing files (IRR_RUN and IRR_WITH) describing the runoff and water availability in current grid cell.

- Vegetation parameter file: An example of a cell in the veg parameter file with irrigated vegetation. There is one extra flag (0 or 1) on the first line of each tile. In this case, the last tile (veg number 79, defined in the veg library file) has the flag set to 1, which means this tile will be irrigated and the varying crop fraction will be read from crop_frac forcing file.

>38597 6
<br>7 0.1094	 0.30  0.60  0.70  0.40  0
<br>0.5870  1.0750  1.3500  1.5870  1.7870  1.7120  1.7880  2.1000  2.0370  1.5500  1.0250  0.6620
<br>8 0.2149	 0.30  0.70  0.70  0.30  0
<br>0.2620  0.2870  0.4000  0.4750  0.3000  0.2500  0.4870  0.9250  0.5620  0.2750  0.2250  0.2620
<br>9 0.4242	 0.30  0.70  0.70  0.30  0
<br>0.5500  0.6120  0.7620  0.4500  0.2250  0.3120  0.7250  0.9380  0.4250  0.2000  0.2000  0.3630
<br>10 0.2233	 0.30  0.80  0.70  0.20  0
<br>0.2370  0.2750  0.3880  0.5870  0.8500  1.0500  0.5870  0.4130  0.3130  0.2870  0.2500  0.2250
<br>11 0.0003	 0.30  0.50  0.70  0.50  0
<br>0.4300  0.8080  1.3730  1.7830  1.2110  0.3660  0.6290  0.9010  0.5010  0.2400  0.2250  0.2960
<br>79 0.0229	 0.30  0.50  0.70  0.50  1
<br>1.4460  4.2880  4.5380  2.4280  1.2810  0.3810  1.2810  1.2810  1.2810  1.2810  0.1330  0.6290

- Crop fraction file: Crop fraction file is needed as an additional forcing file when changing crop fraction is simulated. Similar to the climate forcing file, the crop fraction file is one file per grid cell, describing the fraction change (0-1) over time in the crop tile for each grid cell. Note that the fraction only represents the crop fraction in the crop tile, not in the entire grid cell.

- Vegetation library: There are 12 more columns (JAN-IRR to DEC-IRR) set to 0 or 1 for each veg type, indicating if this veg type will be irrigated when IRRGATION is set to TRUE in global file.

- Soil file: Same as regular soil file.
