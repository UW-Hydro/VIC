VIC 4.2.a with irrigation. 

Developed by Ingjerd Haddeland and Ted Bohn. This stand alone version only possible to model irrigation requirements ("potential irrigation"). Actual irrigation, i.e. with water limitation included, requires use of the routing (reservoir) model with dams included.

####################################
Input files differences compared to regular VIC inputs

Global file: Two extra lines are included.
IRRIGATION      FALSE/TRUE 
IRR_FREE        FALSE/TRUE 
IRRIGATION means if you want irrigation included; IRR_FREE means if the water for irrigation is freely avaiable or not. If IRR_FREE set to FALSE, two extra columns are needed in forcing files (IRR_RUN and IRR_WITH) describing the runoff and water avaiability in current grid cell.

Vegetation parameter file: 

An example of a cell in the veg param file cell with irrigated vegetation can be seen below. There is one extra flag (0 or 1) on the first line of each tile. In this case, the last tile (veg number 79, defined in the veg library file) has the flag set to 1, which means this tile will be irrigated and the varying crop fraction will be read from crop_frac forcing file.
 
38597 6
 7 0.1094	 0.30  0.60  0.70  0.40  0
 0.5870  1.0750  1.3500  1.5870  1.7870  1.7120  1.7880  2.1000  2.0370  1.5500  1.0250  0.6620 
 8 0.2149	 0.30  0.70  0.70  0.30  0
 0.2620  0.2870  0.4000  0.4750  0.3000  0.2500  0.4870  0.9250  0.5620  0.2750  0.2250  0.2620 
 9 0.4242	 0.30  0.70  0.70  0.30  0
 0.5500  0.6120  0.7620  0.4500  0.2250  0.3120  0.7250  0.9380  0.4250  0.2000  0.2000  0.3630 
 10 0.2233	 0.30  0.80  0.70  0.20  0
 0.2370  0.2750  0.3880  0.5870  0.8500  1.0500  0.5870  0.4130  0.3130  0.2870  0.2500  0.2250 
 11 0.0003	 0.30  0.50  0.70  0.50  0
 0.4300  0.8080  1.3730  1.7830  1.2110  0.3660  0.6290  0.9010  0.5010  0.2400  0.2250  0.2960 
 79 0.0229	 0.30  0.50  0.70  0.50  1
 1.4460  4.2880  4.5380  2.4280  1.2810  0.3810  1.2810  1.2810  1.2810  1.2810  0.1330  0.6290

Vegetation library: There are 12 more columns (JAN-IRR to DEC-IRR) set to 0 or 1 for each veg type, indicating if this veg type will be irrigated when IRRGATION is set to TRUE in global file. 

Soil file: Same as regular soil file.


