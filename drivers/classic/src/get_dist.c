#include <stdlib.h>
#include <stdio.h>
#include <vicNl.h>

static char vcid[] = "$Id$";

double get_dist(double lat1, double long1, double lat2, double long2)
/*******************************************************************************
  Function: double get_dist(double lat1, double long1, double lat2, double long2)
  Returns : distance between two locations

  Modifications:
  2007-Nov-06 Moved to separate file from read_lakeparam.c.		TJB
********************************************************************************/
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

