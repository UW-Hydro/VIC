#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
#include <string.h>

static char vcid[] = "$Id$";

#if LAKE_MODEL

lake_con_struct read_lakeparam(FILE            *lakeparam, 
			       soil_con_struct  soil_con, 
			       float            res,
			       double          *Cv_sum)
/**********************************************************************
	read_lakeparam		Laura Bowling		2000

  This routine reads in lake parameters for the current grid cell.  It 
  will either calculate the lake area v. depth profile from a parabolic 
  curve or read in constant values depending on the LAKE_PROFILE flag.

  Parameters Read from File:
  TYPE   NAME                    UNITS   DESCRIPTION
  double eta_a                   kPa/m (?) Decline of solar rad. w/ depth        
  double  maxdepth                 m      Maximum lake depth   
  int    numnod                    -      Number of lake profile nodes.
  double *surface                  m^2    Area of lake at each node. 
  double    b                      -      Exponent controlling lake depth
                                          profile (y=Ax^b)
  float rpercent;                  -      Fraction of the grid cell runoff 
                                          routed through the lake.
  float bpercent;                  -      Fraction of the grid cell baseflow
                                          routed through the lake. 

  Parameters Computed from those in the File:
  TYPE   NAME                    UNITS   DESCRIPTION
  double  dz                       m      Thickness of each solution layer.
  double *surface                 m^2    Area of the lake at each node.
  double cell_area                m^2    Area of the model grid cell.


  Modifications:
  03-11-01 Modified Cv_sum so that it includes the lake fraction,
	   thus 1 - Cv_sum is once again the bare soil fraction.  KAC
  04-Oct-04 Merged with Laura Bowling's updated lake model code.	TJB
  
**********************************************************************/

{
  extern option_struct   options;
#if LINK_DEBUG
  extern debug_struct    debug;
#endif

  int    i;
  int    junk, flag;
  double lat, lng;
  double start_lat, right_lng, left_lng, delta;
  double dist;
  double tempdz;
  double radius, A, x, y;
  char   instr[251];
  char   tmpstr[MAXSTRING+1];

  lake_con_struct temp;
  
#if !NO_REWIND
  rewind(lakeparam);
#endif // NO_REWIND
    
  /*******************************************************************/
  /* Calculate grid cell area. */
  /******************************************************************/

  lat = fabs(soil_con.lat);
  lng = fabs(soil_con.lng);

  start_lat = lat - res / 2.;
  right_lng = lng + res / 2;
  left_lng  = lng - res / 2;

  delta = get_dist(lat,lng,lat+res/10.,lng);

  dist = 0.;

  for ( i = 0; i < 10; i++ ) {
    dist += get_dist(start_lat,left_lng,start_lat,right_lng) * delta;
    start_lat += res/10;
  }
  temp.cell_area = dist * 1000. * 1000.; /* Grid cell area in m^2. */

  /*******************************************************************/
  /* Read in general lake parameters.                           */
  /******************************************************************/

  /* VIC format files: */
  
 /*   fgets(instr, 2500, lakeparam); */

  // Locate current grid cell
  fscanf(lakeparam, "%d", &temp.gridcel);
  while ( temp.gridcel != soil_con.gridcel ) {
    fgets(tmpstr, MAXSTRING, lakeparam);
    fscanf(lakeparam, "%d", &temp.gridcel);
  }

  // cell number not found
  if ( feof(lakeparam) ) {
    sprintf(tmpstr, "Unable to find cell %i in the lake parameter file, check the file or set NO_REWIND to FALSE", soil_con.gridcel);
    nrerror(tmpstr);
  }

  // read lake parameters from file
  fscanf(lakeparam, "%lf %d", &temp.maxdepth, &temp.numnod);
  fscanf(lakeparam, "%lf %lf", &temp.mindepth, &temp.maxrate);

  // Comment out the following if depthfrac and ratefrac aren't used. 
 // fscanf(lakeparam, "%lf %lf", &temp.ratefrac, &temp.depthfrac);

  fscanf(lakeparam, "%lf", &temp.depth_in);
  fscanf(lakeparam, "%f", &temp.rpercent);
  temp.wetland_veg_class = 0;
  temp.bpercent=0.0;
  
  if(temp.numnod >= MAX_LAKE_NODES) {
    nrerror("Number of lake nodes exceeds the maximum allowable, edit user_def.h.");
  }

  if(temp.depth_in > temp.maxdepth) {
    nrerror("Initial depth exceeds the specified maximum lake depth.");
  }

  /**********************************************
      Compute water layer thickness
  **********************************************/
    tempdz = (temp.maxdepth) / ((float) temp.numnod); 

  /*******************************************************************/
  /* Find lake basin area with depth.                           */
  /******************************************************************/

  /* Read in parameters to calculate lake profile. */
  if(!options.LAKE_PROFILE) { 

    fprintf(stderr, "LAKE PROFILE being computed. \n");

    fscanf(lakeparam, "%lf", &temp.Cl[0]);
 if(temp.Cl[0] < 0.0 || temp.Cl[0] > 1.0)
	nrerror("Lake area must be a fraction between 0 and 1, check the lake parameter file.");
    
    fgets(tmpstr, MAXSTRING, lakeparam);
    	
    temp.basin[0] = temp.Cl[0] * temp.cell_area;
	
    /**********************************************
    Compute depth area relationship.
    **********************************************/
  
    radius = sqrt(temp.basin[0] / PI);
    temp.z[0] = temp.maxdepth;
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
    for ( i = 0; i < temp.numnod; i++ ) {
      fscanf(lakeparam, "%lf %lf", &temp.z[i], &temp.Cl[i]);
      temp.basin[i] = temp.Cl[i] * temp.cell_area;
      
//      if(i==0)
//	temp.z[i] = temp.maxdepth;
//      else
//	temp.z[i] = (temp.numnod - i) * tempdz;

      if(temp.Cl[i] < 0.0 || temp.Cl[i] > 1.0)
	nrerror("Lake area must be a fraction between 0 and 1, check the lake parameter file.");
    }
    temp.z[temp.numnod] = 0.0;
    temp.basin[temp.numnod] = 0.0;

    for ( i = 1; i <= temp.numnod; i++ ) {
         temp.maxvolume += (temp.basin[i] + temp.basin[i-1]) * tempdz/2.;
	    }
   fprintf(stderr, "maxvolume = %e km3\n",temp.maxvolume/(1000.*1000.*1000.));
	 }

  // Add lake fraction to Cv_sum
  (*Cv_sum) += temp.Cl[0];
  if ( *Cv_sum > 0.999 ) {
    // Adjust Cv_sum so that it is equal to 1 
    temp.Cl[0] += 1. - *Cv_sum;
    *Cv_sum = 1;
  }
  fprintf(stderr, "Lake area = %e km2\n",temp.basin[0]/(1000.*1000.));
  return temp;
}

/*******************************************************************************
  Function: double distance(double lat1, double long1, double lat2, double long2)
  Returns : distance between two locations
********************************************************************************/

double get_dist(double lat1, double long1, double lat2, double long2)
{
  double theta1;
  double phi1;
  double theta2;
  double phi2;
  double dtor;
  double term1;
  double term2;
  double term3;
  double temp;
  double distance;

  dtor = 2.0*PI/360.0;
  theta1 = dtor*long1;
  phi1 = dtor*lat1;
  theta2 = dtor*long2;
  phi2 = dtor*lat2;
  term1 = cos(phi1)*cos(theta1)*cos(phi2)*cos(theta2);
  term2 = cos(phi1)*sin(theta1)*cos(phi2)*sin(theta2);
  term3 = sin(phi1)*sin(phi2);
  temp = term1+term2+term3;
  temp = (double) (1.0 < temp) ? 1.0 : temp;
  distance = RADIUS*acos(temp);

  return distance;
}  

#endif // LAKE_MODEL
