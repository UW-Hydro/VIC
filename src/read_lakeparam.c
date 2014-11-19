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
  2008-Mar-01 Moved assignment of tempdz so that it is always assigned a value.	TJB
  2008-Jun-16 Added a second fgets to loop over grid cells, to correctly parse
	      file.								LCB via TJB
  2008-Sep-09 Deleted the fprintf statement for maxiumum lake volume.		LCB via TJB
  2009-Jul-31 Added new parameter: index of veg tile (from veg param file) that
	      contains the lake/wetland.					TJB
  2010-Sep-24 Clarified the fprintf statement describing lake basin area as
	      lake plus wetland area.						TJB
  2010-Oct-10 Corrected the check on initial depth to allow initial depth <
	      mindepth.								TJB
  2012-Jan-02 Modified to turn off lakes in a grid cell if lake_idx is < 0.
	      Added validation of parameter values.				TJB
  2012-Jan-16 Removed LINK_DEBUG code						BN
  2013-Jul-25 Fixed bug in parsing lakeparam file in case of no lake
	      in the cell.							TJB
  2013-Dec-28 Removed NO_REWIND option.					TJB
**********************************************************************/

{
  extern option_struct   options;
  int    i;
  int    lakecel;
  int    junk, flag;
  double tempdz;
  double radius, A, x, y;
  char   instr[251];
  char   tmpstr[MAXSTRING+1];
  int    ErrFlag;
  double tmp_mindepth, tmp_maxdepth;

  lake_con_struct temp;

  /*******************************************************************/
  /* Read in general lake parameters.                           */
  /******************************************************************/

  fscanf(lakeparam, "%d %d", &lakecel, &temp.lake_idx);
  while ( lakecel != soil_con.gridcel && !feof(lakeparam) ) {
    fgets(tmpstr, MAXSTRING, lakeparam); // grid cell number, etc.
    if (temp.lake_idx >= 0)
      fgets(tmpstr, MAXSTRING, lakeparam); // lake depth-area relationship
    fscanf(lakeparam, "%d %d", &lakecel, &temp.lake_idx);
  }

  // cell number not found
  if ( feof(lakeparam) ) {
    sprintf(tmpstr, "Unable to find cell %i in the lake parameter file", soil_con.gridcel);
    nrerror(tmpstr);
  }

  // read lake parameters from file
  if (temp.lake_idx >= 0) {
    veg_con[temp.lake_idx].LAKE = 1;
    fscanf(lakeparam, "%d", &temp.numnod);
    if (temp.numnod < 1) {
      sprintf(tmpstr, "Number of vertical lake nodes (%d) for cell %d specified in the lake parameter file is < 1; increase this number to at least 1.", temp.numnod, soil_con.gridcel);
      nrerror(tmpstr);
    }
    if(temp.numnod > MAX_LAKE_NODES) {
      sprintf(tmpstr, "Number of lake nodes (%d) in cell %d specified in the lake parameter file exceeds the maximum allowable (%d), edit MAX_LAKE_NODES in user_def.h.", temp.numnod, soil_con.gridcel, MAX_LAKE_NODES);
      nrerror(tmpstr);
    }
    fscanf(lakeparam, "%lf", &temp.mindepth);
    if (temp.mindepth < 0) {
      sprintf(tmpstr, "Minimum lake depth (%f) for cell %d specified in the lake parameter file is < 0; increase this number to at least 0.", temp.mindepth, soil_con.gridcel);
      nrerror(tmpstr);
    }
    fscanf(lakeparam, "%lf", &temp.wfrac);
    if (temp.wfrac < 0 || temp.wfrac > 1) {
      sprintf(tmpstr, "Lake outlet width fraction (%f) for cell %d specified in the lake parameter file falls outside the range 0 to 1.  Change wfrac to be between 0 and 1.", temp.wfrac, soil_con.gridcel);
      nrerror(tmpstr);
    }
    fscanf(lakeparam, "%lf", &temp.depth_in);
    if (temp.depth_in < 0) {
      sprintf(tmpstr, "Initial lake depth (%f) for cell %d specified in the lake parameter file is < 0; increase this number to at least 1.", temp.depth_in, soil_con.gridcel);
      nrerror(tmpstr);
    }
    fscanf(lakeparam, "%f", &temp.rpercent);
    if (temp.rpercent < 0 || temp.rpercent > 1) {
      sprintf(tmpstr, "Fraction of runoff entering lake catchment (%f) for cell %d specified in the lake parameter file falls outside the range 0 to 1.  Change rpercent to be between 0 and 1.", temp.rpercent, soil_con.gridcel);
      nrerror(tmpstr);
    }
  }
  else { // no lake exists anywhere in this grid cell
    temp.numnod = 0;
    temp.mindepth = 0;
    temp.maxdepth = 0;
    temp.Cl[0] = 0;
    temp.basin[0] = 0;
    temp.z[0] = 0;
    temp.minvolume = 0;
    temp.maxvolume = 0;
    temp.wfrac = 0;
    temp.depth_in = 0;
    temp.rpercent = 0;
    temp.bpercent = 0;
    fgets(tmpstr, MAXSTRING, lakeparam); // skip to end of line
    return temp;
  }
  temp.bpercent = temp.rpercent;

  /*******************************************************************/
  /* Find lake basin area with depth.                           */
  /******************************************************************/

  /* Read in parameters to calculate lake profile. */
  if(!options.LAKE_PROFILE) { 

    fprintf(stderr, "LAKE PROFILE being computed. \n");

    fscanf(lakeparam, "%lf %lf", &temp.z[0], &temp.Cl[0]);
    temp.maxdepth = temp.z[0];
    tempdz = (temp.maxdepth) / ((float) temp.numnod);
    if(temp.Cl[0] < 0.0 || temp.Cl[0] > 1.0) {
      sprintf(tmpstr, "Lake area fraction (%f) for cell (%d) specified in the lake parameter file must be a fraction between 0 and 1.", temp.Cl[0], soil_con.gridcel);
      nrerror(tmpstr);
    }
    
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
      
      if(i==0) {
        temp.maxdepth = temp.z[i];
        tempdz = (temp.maxdepth) / ((float) temp.numnod);
      }

      if(temp.Cl[0] < 0.0 || temp.Cl[0] > 1.0) {
        sprintf(tmpstr, "Lake area fraction (%f) for cell (%d) specified in the lake parameter file must be a fraction between 0 and 1.", temp.Cl[0], soil_con.gridcel);
        nrerror(tmpstr);
      }
    }
    temp.z[temp.numnod] = 0.0;
    temp.basin[temp.numnod] = 0.0;
    temp.Cl[temp.numnod] = 0.0;

    for ( i = 1; i <= temp.numnod; i++ ) {
      temp.maxvolume += (temp.basin[i] + temp.basin[i-1]) * (temp.z[i-1] - temp.z[i]) / 2.;
    }

  }

  // Compute volume corresponding to mindepth
  ErrFlag = get_volume(temp, temp.mindepth, &(temp.minvolume));
  if (ErrFlag == ERROR) {
    sprintf(tmpstr, "ERROR: problem in get_volume(): depth %f volume %f rec %d\n", temp.mindepth, temp.minvolume, 0);
    nrerror(tmpstr);
  }

  // Make sure min < max
  if(temp.mindepth > temp.maxdepth) {
    sprintf(tmpstr, "Adjusted minimum lake depth %f exceeds the adjusted maximum lake depth %f.", temp.mindepth, temp.maxdepth);
    nrerror(tmpstr);
  }

  // Validate initial conditions
  if(temp.depth_in > temp.maxdepth) {
    fprintf(stderr, "WARNING: Initial lake depth %f exceeds the maximum lake depth %f; setting initial lake depth equal to max lake depth.\n", temp.depth_in, temp.maxdepth);
    temp.depth_in = temp.maxdepth;
  }
  else if(temp.depth_in < 0) {
    fprintf(stderr, "WARNING: Initial lake depth %f < 0; setting to 0.\n", temp.depth_in);
    temp.depth_in = 0;
  }

  fprintf(stderr, "Lake plus wetland area = %e km2\n",temp.basin[0]/(1000.*1000.));
  return temp;
}

