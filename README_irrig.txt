VIC 4.0.6 with irrigation. 

VIC stand-alone version: Only possible to model irrigation requirements ("potential irrigation"). Actual irrigation, i.e. with water limitation included, requires use of the routing model with dams included in addition.

Limitations when irrigation is included: Only tested in daily water balance mode, one elevation band, and without distributed precipitation. Also, you can get water balance errors at the beginning of the second simulation month (has to do with initialization), so always include a spin up period! Does not allow for the use of INIT_STATE. You can only have one irrigated vegetation type within each cell.

PotEvap does not funtction properly!
####################################
Input files

Global file: Two extra lines are included. For now, set both these to TRUE if you want irrigation included (IRRIGATION = TRUE and IRR_FREE = FALSE means irrigation with water limitations, which cannot be done with the VIC irrigation stadn-alone version). 
IRRIGATION      FALSE/TRUE 
IRR_FREE        FALSE/TRUE 

Vegetation parameter file: 
A cell without irrigated vegetation should look like this (i.e. one extra 0 on the first line):
20000 5 0
	1 0.1608	0.30 0.30 0.70 0.70 
	5.2130 5.2130 5.2130 5.2130 5.2130 5.2130 5.2130 5.2130 5.2130 5.2130 5.2130 5.2130 
	5 0.4148	0.30 0.30 0.70 0.70 
	0.1000 0.1000 0.1620 0.4620 3.6620 5.2750 5.5880 5.3630 3.0750 0.6500 0.2000 0.1500 
	6 0.2401	0.30 0.60 0.70 0.40 
	0.2870 0.3120 0.3120 0.3120 1.2500 3.4630 4.3380 3.5250 2.0880 0.7620 0.6130 0.4620 
	7 0.1284	0.30 0.60 0.70 0.40 
	0.3250 0.2500 0.2500 0.2500 1.0620 3.1500 3.8750 3.0870 1.1750 0.4250 0.4250 0.4250 
	11 0.0558	0.30 0.50 0.70 0.50 
	0.1500 0.1620 0.2880 0.2640 3.0110 4.6190 4.4940 4.2060 2.3230 0.6900 0.1900 0.1900 
An example of a cell in the veg param file cell with irrigated vegetation can be seen below. There should be one extra 1 on the first line, + information on percent irrigated area pr month on last line. Number of vegetation types is, as you can see, only 6, although there may seem to be information on 7 types listed. The vegetation type "110" is irrigated vegetation. The percentage tells the model how much of the area equipped for irrigation within that cell (0.0487+0.0487) that is actually irrigated that month (information taken from e.g. FAO). Do not go below 1 percent or above 99 percent on the last line!  
20376 6 1
	1 0.0410	0.30 0.30 0.70 0.70 
	4.8780 4.8780 4.8780 4.8780 4.8780 4.8780 4.8780 4.8780 4.8780 4.8780 4.8780 4.8780 
	6 0.0608	0.30 0.60 0.70 0.40 
	0.3380 0.3630 0.3630 0.3750 0.7880 2.1250 3.1500 2.6370 1.4620 0.6250 0.3870 0.3870 
	7 0.2061	0.30 0.60 0.70 0.40 
	0.4630 0.6500 0.9880 1.9380 3.3630 3.2370 2.6630 1.8120 1.9250 1.1380 0.6130 0.3380 
	10 0.0230	0.30 0.80 0.70 0.20 
	0.2250 0.3250 0.4750 0.8370 2.4000 3.1380 2.9380 2.1250 1.7000 1.0500 0.4500 0.2500 
	11 0.5719	0.30 0.50 0.70 0.50 
	0.1080 0.1560 0.2470 0.5410 2.1790 3.6170 3.1640 1.2390 0.9970 0.6120 0.3250 0.1750 
	110 0.0487	0.30 0.50 0.70 0.50 
	0.1000 0.1000 0.1000 0.1000 0.2230 0.9340 2.2950 2.4310 1.0530 0.1000 0.1000 0.1000 
	110 0.0487	0.30 0.50 0.70 0.50 
	0.1000 0.1000 0.1000 0.1000 0.2230 0.9340 2.2950 2.4310 1.0530 0.1000 0.1000 0.1000 
	10.0 10.0 10.0 10.0 80.0000 80.0 80.0 80.0 80.0 10.0 10.0 10.0 

Vegetation library: Same as 4.0.6 (but you may have to add lines, e.g. I've added information on the vegetation type "110". 

Soil file: Three extra columns added, although they have nothing to do with the irrigation scheme, - I just needed them for other purposes). Given three soil layers, the added cols are: Col 54: Not used (although it may seem so if you read read_soilparam.c). Col 55: options.BASEFLOW (overrules info in global file). Col 56: options.ROOT_ZONES (overrules info in global file). 
