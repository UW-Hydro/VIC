#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
#include <math.h>

static char vcid[] = "$Id: gl_flow.c,v 5.12.2.20 2013/10/7 03:45:12 vicadmin Exp $";


int gl_flow(soil_con_struct     *soil_con,
	        veg_con_struct      *veg_con)

/**********************************************************************
	gl_flow      Bibi S. Naz	October 07, 2013

  this routine solve the glacier flow algorithm as part of VIC glacier flow model. see more detail for ...... 

 						
**********************************************************************/
{
  extern veg_lib_struct *veg_lib;
  extern option_struct   options;
#if LINK_DEBUG
  extern debug_struct    debug;
#endif // LINK_DEBUG
  double kterminus;
  double  Hhi;
  double gradhi;
  double areahi;
  double  Hlo;
  double gradlo;
  double arealo;
  double Qnet;
  double delH;
  double divQ;
  double gamma;
  double delswe, deliwe;
  int rhoi=910;
  int rhow=1000;
  int Lvert=100;    // m - length scale for surface slopes/ice flux
  int                    iveg;
  int                    Nveg;
  int                    veg_class;
  int                    band;
  int                    Nbands, Glbands;
  
  double                 Cv;
  snow_data_struct     **snow;
  gamma=1.0e-8;                   // tune this value for the ice dynamics
 
  for(iveg = 0; iveg <= Nveg; iveg++){
    if (veg_con[iveg].Cv > 0.0) {
      Cv = veg_con[iveg].Cv;
      Nbands = options.SNOW_BAND;
      // find the highest ice
      for (band = 0; band <=Nbands; band++) {	
	printf("gl_area=%f\n", snow[iveg][band].iwq);
	if ((soil_con->GlAreaFract[band] > 0) && (soil_con->GlAreaFract[band+1]==0))
	  {
	    Glbands=band;
	    printf("glbands=%d\n",Glbands);
	  }
	else
	  Glbands=Nbands;
      }
	// look for the first occurence and flag the terminus
      //for (band = Glbands; band >=0; band--) {
      for (band = 0; band <=Glbands; band++) {		
	if( soil_con->AreaFract[band] > 0 ) {
	  printf("band= %d\n", band);
	  //printf("swq= %f\n", snow[iveg][band].iwq);
	  delswe = snow[iveg][band].swq - snow[iveg][band].swqold;  //change in swe in previous month
	  snow[iveg][band].swqold = snow[iveg][band].swq;
	  deliwe = snow[iveg][band].iwq - snow[iveg][band].iwqold; // change in iwe in previous month
	  snow[iveg][band].iwqold = snow[iveg][band].iwq;  
	  snow[iveg][band].bn = delswe + deliwe; //glacier mass balance 
	  
	  if (band > 0)
	    {
	      if ((snow[iveg][band].iwq == 0) && (snow[iveg][band-1].iwq > 0))
		kterminus=band-1;
	    }
	  
	  // calculate the flux into the cell
	  if (band < Glbands)   // below the top cell
	    {
	      if(snow[iveg][band+1].iwq > 0) //ice above
		{
		  Hhi = (snow[iveg][band].iwq + snow[iveg][band+1].iwq)/2.;
		  gradhi = (soil_con->BandElev[band+1]-soil_con->BandElev[band])/Lvert;
		  areahi = (soil_con->GlAreaFract[band]+soil_con->GlAreaFract[band+1])/2;
		}
	      else
		{
		  Hhi=0;
		  gradhi=0;
		  areahi=0;
		}
	    }
	  else // top cell
	    {
	      if(snow[iveg][band].iwq > 0) //ice here
		{
		  Hhi = snow[iveg][band].iwq;
		  gradhi = 0; // shuts off incoming ice
		  areahi = 0;
		}
	      else
		{
		  Hhi=0;
		  gradhi=0;
		  areahi=0;
		}
	    }
	  // calculate the flux out of the cell
	  if (band > 0)     // above the valley bottom
	    {
	      gradlo=(soil_con->BandElev[band]-soil_con->BandElev[band-1])/Lvert;
	      Hlo = (snow[iveg][band].iwq + snow[iveg][band-1].iwq)/2;
	      //arealo=(hyps_future(k)+hyps_future(k-1))/2;
	      arealo=(soil_con->GlAreaFract[band]+soil_con->GlAreaFract[band-1])/2;
	    }
	  else
	    {
	      Hlo=snow[iveg][band].iwq;
	      gradlo=1;
	      arealo=soil_con->GlAreaFract[band];
	    }
	  
	  // assign fluxes at the interface
	  snow[iveg][band].Qin = gamma*areahi * pow(Hhi,5)* pow(gradhi,3);
	  snow[iveg][band].Qout = gamma*arealo * pow(Hlo,5) * pow(gradlo,3);
	  Qnet = snow[iveg][band].Qin - snow[iveg][band].Qout;
	  if ((Qnet > 0) && (soil_con->GlAreaFract[band]==0))
	    soil_con->GlAreaFract[band] = soil_con->GlAreaFract[band+1];   // new incoming ice
	  
	  if (abs(Qnet) > 0)
	    {
	      if ((Qnet/soil_con->GlAreaFract[band])<5)
		divQ = Qnet/soil_con->GlAreaFract[band];
	      else
		divQ = 5;
	      
	    }
	  
	  else
	    divQ=0;
	  
	  delH = snow[iveg][band].bn + divQ;
	  if((snow[iveg][band].iwq + delH) > 0)
	    snow[iveg][band].iwq = snow[iveg][band].iwq + delH;
	  else
	    snow[iveg][band].iwq = 0;
	  
	  soil_con->BandElev[band] = soil_con->BandElev[band] + delH;  
	  // don't melt into the rock
	  if (snow[iveg][band].iwq < 0){
	    snow[iveg][band].iwq= 0;
	    soil_con->GlAreaFract[band]= 0;
	  }
	} 
      }// End loop through elevation bands  
    } /** end current vegetation type **/
  } /** end of vegetation loop **/
  return(0);
}
  
  
