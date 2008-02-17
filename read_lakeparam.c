#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
#include <string.h>

static char vcid[] = "$Id$";

lake_con_struct read_lakeparam(FILE            *lakeparam, 
			       soil_con_struct  soil_con, 
			       veg_con_struct  *veg_con)
/**********************************************************************
	read_lakeparam		Laura Bowling		2000

  This routine reads in lake parameters for the current grid cell.  It 
  will either calculate the lake area v. depth profile from a parabolic 
  curve or read in constant values depending on the LAKE_PROFILE flag.

  Parameters Read from File:
  TYPE   NAME                    UNITS   DESCRIPTION
  double  maxdepth                 m      Maximum lake depth   
  int     numnod                   -      Number of lake profile nodes.
  double *surface                  m^2    Area of lake at each node. 
  double  b                        -      Exponent controlling lake depth
                                          profile (y=Ax^b)
  float rpercent;                  -      Fraction of the grid cell runoff 
                                          routed through the lake.
  float bpercent;                  -      Fraction of the grid cell baseflow
                                          routed through the lake. 

  Parameters Computed from those in the File:
  TYPE   NAME                    UNITS   DESCRIPTION
  double  dz                       m      Thickness of each solution layer.
  double *surface                  m^2    Area of the lake at each node.

  Modifications:
  03-11-01 Modified Cv_sum so that it includes the lake fraction,
	   thus 1 - Cv_sum is once again the bare soil fraction.  KAC
  04-Oct-04 Merged with Laura Bowling's updated lake model code.		TJB
  02-Nov-04 Modified the adjustment of Cv_sum so that veg fractions
	    share in area reduction along with the lake fraction.		TJB
  22-Feb-05 Merged with Laura Bowling's second update to lake model code.	TJB
  2005-03-08 Added EQUAL_AREA option.  When TRUE, res variable is interpreted
	     to be the grid cell area in km^2.  When FALSE, res is interpreted
	     to be the length of a grid cell side, in degrees (as before).	TJB
  2005-03-17 Laura Bowling's update had included what appeared to be temporary
	     code that expected to read the lake node depths from the lake param
	     file.  This has code has been removed.				TJB
  2006-Nov-07 Removed LAKE_MODEL option.					TJB
  2007-Oct-24 Modified to handle grid cells with empty or very shallow starting
	      lake depths.							KAC via TJB
  2007-Nov-06 Updated to read new set of lake parameters (broad-crested wier)	LCB via TJB
  2008-Feb-16 Added !feof(lakeparam) condition to loop over lakeparam file.	TJB
**********************************************************************/

{
  extern option_struct   options;
#if LINK_DEBUG
  extern debug_struct    debug;
#endif

  int    i;
  int    lakecel;
  int    junk, flag;
  double tempdz;
  double radius, A, x, y;
  char   instr[251];
  char   tmpstr[MAXSTRING+1];
  int    ErrFlag;

  lake_con_struct temp;
  
#if !NO_REWIND
  rewind(lakeparam);
#endif // NO_REWIND
    
  /*******************************************************************/
  /* Read in general lake parameters.                           */
  /******************************************************************/

  // Locate current grid cell
  fscanf(lakeparam, "%d", &lakecel);
  while ( lakecel != soil_con.gridcel && !feof(lakeparam) ) {
    fgets(tmpstr, MAXSTRING, lakeparam);
    fscanf(lakeparam, "%d", &lakecel);
  }

  // cell number not found
  if ( feof(lakeparam) ) {
    sprintf(tmpstr, "Unable to find cell %i in the lake parameter file, check the file or set NO_REWIND to FALSE", soil_con.gridcel);
    nrerror(tmpstr);
  }

  // read lake parameters from file
  fscanf(lakeparam, "%d", &temp.numnod);
  fscanf(lakeparam, "%lf %lf", &temp.mindepth, &temp.wfrac);
  fscanf(lakeparam, "%lf", &temp.depth_in);
  fscanf(lakeparam, "%f", &temp.rpercent);
  temp.wetland_veg_class = 0;
  temp.bpercent=0.0;

  if(temp.numnod > MAX_LAKE_NODES) {
    nrerror("Number of lake nodes exceeds the maximum allowable, edit user_def.h.");
  }

  /*******************************************************************/
  /* Find lake basin area with depth.                           */
  /******************************************************************/

  /* Read in parameters to calculate lake profile. */
  if(!options.LAKE_PROFILE) { 

    fprintf(stderr, "LAKE PROFILE being computed. \n");

    fscanf(lakeparam, "%lf %lf", &temp.z[0], &temp.Cl[0]);
    temp.maxdepth = temp.z[0];
    if(temp.Cl[0] < 0.0 || temp.Cl[0] > 1.0)
      nrerror("Lake area must be a fraction between 0 and 1, check the lake parameter file.");
    
    fgets(tmpstr, MAXSTRING, lakeparam);
    	
    temp.basin[0] = temp.Cl[0] * soil_con.cell_area;
	
    /**********************************************
    Compute depth area relationship.
    **********************************************/
  
    radius = sqrt(temp.basin[0] / PI);

    temp.maxvolume = 0.0;
    for(i=1; i<= temp.numnod; i++) {
      temp.z[i] = (temp.numnod - i) * tempdz;             
      if(temp.z[i] < 0.0) temp.z[i] = 0.0;
      x = pow(temp.z[i]/temp.maxdepth,BETA)*radius;
      temp.basin[i] = PI * x * x;	
      temp.maxvolume += (temp.basin[i] + temp.basin[i-1]) * tempdz/2.;
    }

  }
  
  /* Read in basin area for each layer depth. */
  /* Assumes that the lake bottom area is zero. */

  else{       

    fprintf(stderr, "Reading in the specified lake profile.\n");
    temp.maxvolume=0.0;
    temp.Cl[0] = 0; // initialize to 0 in case no lake is defined
    for ( i = 0; i < temp.numnod; i++ ) {
      fscanf(lakeparam, "%lf %lf", &temp.z[i], &temp.Cl[i]);
      temp.basin[i] = temp.Cl[i] * soil_con.cell_area;
      
      if(i==0)
        temp.maxdepth = temp.z[i];

      if(temp.Cl[i] < 0.0 || temp.Cl[i] > 1.0)
	nrerror("Lake area must be a fraction between 0 and 1, check the lake parameter file.");
    }
    temp.z[temp.numnod] = 0.0;
    temp.basin[temp.numnod] = 0.0;
    temp.Cl[temp.numnod] = 0.0;

    for ( i = 1; i <= temp.numnod; i++ ) {
      temp.maxvolume += (temp.basin[i] + temp.basin[i-1]) * (temp.z[i-1] - temp.z[i]) / 2.;
    }
    fprintf(stderr, "maxvolume = %e km3\n",temp.maxvolume/(1000.*1000.*1000.));

  }

  if(temp.depth_in > temp.maxdepth) {
    nrerror("Initial depth exceeds the specified maximum lake depth.");
  }

  ErrFlag = get_volume(temp, temp.mindepth, &(temp.minvolume));
  if (ErrFlag == ERROR) {
    fprintf(stderr, "ERROR: problem in get_volume(): depth %f volume %f rec %d\n", temp.mindepth, temp.minvolume, 0);
    exit(0);
  }

  /* Compute water layer thickness */
  tempdz = (temp.maxdepth) / ((float) temp.numnod);

  // Add lake fraction to Cv_sum
  (veg_con[0].Cv_sum) += temp.Cl[0];
  if ( veg_con[0].Cv_sum > 1.000 ) {
    // Adjust Cv_sum so that it is equal to 1
    fprintf(stderr,"WARNING: Total Cv of veg fractions and lake fraction exceeds 1.0 at grid cell %d, fractions being adjusted to equal 1\n", soil_con.gridcel);
    for( i = 0; i < veg_con[0].vegetat_type_num; i++ )
      veg_con[i].Cv = veg_con[i].Cv / veg_con[0].Cv_sum;
    temp.Cl[0] = temp.Cl[0] / veg_con[0].Cv_sum;
    veg_con[0].Cv_sum = 1;
  }

  fprintf(stderr, "Lake area = %e km2\n",temp.basin[0]/(1000.*1000.));
  return temp;
}

