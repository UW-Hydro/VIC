/* to read soil parameters from various soil file directories and put them together as the final soil file - VIC version 4.0.1... 

   NOTE:- the final soil file has 53 columns.... 

   all the properties from the soil program are for two layers of soil..the first 30 cm and the next 100cm. In our model, we have the top two layers as 10 cm and 100cm. but the properties correspond to the two horizons in the soil program...

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define MAX 1000
#define MXR 65812  /* maximum number of grid cells */

#define NOT_REQD  -9999.0  /* no data value */

int main (int argc, char *argv[]) {
  
  int cellid[MXR],ll_cnt,RUN_FLAG,FS_ACT,i,n;

  float lng[MXR],lat[MXR],elev[MXR],CRITICAL,depth1,depth2,depth3_00,SOIL_ROUGH,SNOW_ROUGH,AVG_T,DP,res_moist30,res_moist100,lngarno[MXR],latarno[MXR],d1arno[MXR],d2arno[MXR],d3arno[MXR],d4arno[MXR],latwat[MXR],lngwat[MXR],watn30[MXR],watn100[MXR],dumwat1[MXR],dumwat2[MXR],latks30[MXR],lngks30[MXR],ks30[MXR],latks100[MXR],lngks100[MXR],ks100[MXR],latbd30[MXR],lngbd30[MXR],bd30[MXR],latbd100[MXR],lngbd100[MXR],bd100[MXR],latts30[MXR],lngts30[MXR],ts30[MXR],latts100[MXR],lngts100[MXR],ts100[MXR],latfc[MXR],lngfc[MXR],fc30[MXR],fc100[MXR],latwp[MXR],lngwp[MXR],wp30[MXR],wp100[MXR],latprcp[MXR],lngprcp[MXR],prcp[MXR],ts3,PHI,depth3,fc3,initmo1,initmo2,initmo3,wp3,soildens1,soildens2,soildens3,lngsd30[MXR],latsd30[MXR],lngsd100[MXR],latsd100[MXR],sd30[MXR],sd100[MXR],lngcl30[MXR],latcl30[MXR],lngcl100[MXR],latcl100[MXR],cl30[MXR],cl100[MXR],hb30,hb100;
  float lngcal[MXR],latcal[MXR],calbi[MXR],caln1[MXR],caln2[MXR],calks1[MXR],calks2[MXR],cald2[MXR];
  float BD;
  float INFIL;  /* infiltration parameter */
  float CALIB_CONST_FC_WP; /* FC and WP calibration factor */
  float CALIB_CONST_vG;    /* van Genuchten calibration factor */
  float CALIB_CONST_Ks;    /* Ksat calibration factor */

  char string[MAX+1];

  FILE *fll,*farno,*fwat,*fks30,*fks100,*fbd30,*fbd100,*fts30,*fts100,*ffc,*fwp,*fprcp,*fout,*fcal,*fsd30,*fsd100,*fcl30,*fcl100;

	
  /* usage */
  if(argc!=4) {
    printf("Usage: %s   lng lat file  soil output file prcpfile\n",argv[0]);
    exit(0);
  }
  
  /*open the output file */
  fout=fopen(argv[2],"w");
  
  /* open and read the lng lat file for the mekong */
  if((fll=fopen(argv[1],"r"))==NULL)   {
    printf("ERROR: Unable to open %s\n",argv[1]);
    exit(0);
  }
  
  ll_cnt=0;
  while (!feof(fll)) {
    fgets(string,MAX,fll);
    sscanf(string,"%d %f %f %f",&cellid[ll_cnt],&lng[ll_cnt],&lat[ll_cnt],&elev[ll_cnt]);
    ll_cnt++;
  }
  fclose(fll);

  /* fixed parameters */

  CRITICAL   = 0.8; /*this is set to 0.8 instead of 0.7 in order to keep wcr
		      greater than wwp*/

  SOIL_ROUGH = 0.001;  /*set for arctic - might change for global...was 0.005*/
  SNOW_ROUGH = 0.0005;   /*set for arctic - might change for global...was 0.01*/

  AVG_T  = NOT_REQD;
  PHI    = NOT_REQD;

  FS_ACT = 1;
  

  
  
  

  /* ARNO paramters */
  if((farno=fopen("ARNO/arno.global.param.30min","r"))==NULL) {
    printf("ERROR: Unable to open %s\n","ARNO/arno.global.param.30min");
    exit(0);
  }
  
  i=0;
  
  while (!feof(farno)) {
    fgets(string,MAX,farno);
    sscanf(string,"%f %f %f %f %f %f",&lngarno[i],&latarno[i],&d1arno[i],&d2arno[i],&d3arno[i],&d4arno[i]);
    i++;
  }
  fclose(farno);




  /* Van Genuchten N  - this is actually Brooks-Corey/Campbell*/
  if ((fwat=fopen("HYDR/Lambda_values","r"))==NULL)   {
    printf("ERROR: Unable to open %s\n","HYDR/Lambda_values");
    exit(0);
  }
  
  i=0;
  while (!feof(fwat)) {
    fgets(string,MAX,fwat);
    sscanf(string,"%f %f %f %f %f %f",&lngwat[i],&latwat[i],&dumwat1[i],&watn30[i],&dumwat2[i],&watn100[i]);
    i++;
  }
  fclose(fwat);



  
  /* Hydraulic conductivity */

  if ((fks30=fopen("soil_props/Ks30.masked","r"))==NULL)   {
    printf("ERROR: Unable to open %s\n","soil_props/Ks30.masked");
    exit(0);
  }
  
  i=0;
  while (!feof(fks30)) {
    fgets(string,MAX,fks30);
    sscanf(string,"%f %f %f",&lngks30[i],&latks30[i],&ks30[i]);
    i++;
  }
  fclose(fks30);
  
  if ((fks100=fopen("./soil_props/Ks100.masked","r"))==NULL)   {
    printf("ERROR: Unable to open %s\n","./soil_props/Ks100.masked");
    exit(0);
  }
  
  i=0;
  while (!feof(fks100)) { 
    fgets(string,MAX,fks100);
    sscanf(string,"%f %f %f",&lngks100[i],&latks100[i],&ks100[i]);
    i++;
  }
  fclose(fks100);
  



  /* Bulk density */
  if((fbd30=fopen("./soil_props/BulkDens30.masked","r"))==NULL)   {
    printf("ERROR: Unable to open %s\n","./soil_props/BulkDens30.masked");
    exit(0);
  }

  i=0;
  while (!feof(fbd30)) { 
    fgets(string,MAX,fbd30);
    sscanf(string,"%f %f %f",&lngbd30[i],&latbd30[i],&bd30[i]);
    i++;
  }
  fclose(fbd30);

  if((fbd100=fopen("./soil_props/BulkDens100.masked","r"))==NULL)   {
    printf("ERROR: Unable to open %s\n","./soil_props/BulkDens100.masked");
    exit(0);
  }

  i=0;
  while (!feof(fbd100)) { 
    fgets(string,MAX,fbd100);
    sscanf(string,"%f %f %f",&lngbd100[i],&latbd100[i],&bd100[i]);
    i++;
  }
  fclose(fbd100);



  /* Porosity...ThetaS*/
  if((fts30=fopen("./soil_props/ThetaS30.masked","r"))==NULL)   {
    printf("ERROR: Unable to open %s\n","./soil_props/ThetaS30.masked");
    exit(0);
  }

  i=0;
  while (!feof(fts30)) { 
    fgets(string,MAX,fts30);
    sscanf(string,"%f %f %f",&lngts30[i],&latts30[i],&ts30[i]);
    i++;
  }
  fclose(fts30);
  
  if((fts100=fopen("./soil_props/ThetaS100.masked","r"))==NULL)   {
    printf("ERROR: Unable to open %s\n","./soil_props/ThetaS100.masked");
    exit(0);
  }

  i=0;
  while (!feof(fts100)) { 
    fgets(string,MAX,fts100);
    sscanf(string,"%f %f %f",&lngts100[i],&latts100[i],&ts100[i]);
    i++;
  }
  fclose(fts100);
  

  /* quartz to be represented by sand content...*/
  if ((fsd30=fopen("soil_props/Sand30.masked","r"))==NULL)   {
    printf("ERROR: Unable to open %s\n","soil_props/Sand30.masked");
    exit(0);
  }
  
  i=0;
  while (!feof(fsd30)) { 
    fgets(string,MAX,fsd30);
    sscanf(string,"%f %f %f",&lngsd30[i],&latsd30[i],&sd30[i]);
    i++;
  }
  fclose(fsd30);

  if ((fsd100=fopen("soil_props/Sand100.masked","r"))==NULL)   {
    printf("ERROR: Unable to open %s\n","soil_props/Sand100.masked");
    exit(0);
  }
  
  i=0;
  while (!feof(fsd100)) { 
    fgets(string,MAX,fsd100);
    sscanf(string,"%f %f %f",&lngsd100[i],&latsd100[i],&sd100[i]);
    i++;
  }
  fclose(fsd100);



  /* Percent Clay...*/
  if ((fcl30=fopen("soil_props/Clay30.masked","r"))==NULL)   {
    printf("ERROR: Unable to open %s\n","soil_props/Clay30.masked");
    exit(0);
  }
  
  i=0;
  while (!feof(fcl30)) { 
    fgets(string,MAX,fcl30);
    sscanf(string,"%f %f %f",&lngcl30[i],&latcl30[i],&cl30[i]);
    i++;
  }
  fclose(fcl30);

  if ((fcl100=fopen("soil_props/Clay100.masked","r"))==NULL)   {
    printf("ERROR: Unable to open %s\n","soil_props/Clay100.masked");
    exit(0);
  }
  
  i=0;
  while (!feof(fcl100)) { 
    fgets(string,MAX,fcl100);
    sscanf(string,"%f %f %f",&lngcl100[i],&latcl100[i],&cl100[i]);
    i++;
  }
  fclose(fcl100);


  /* fieldcapacity ..*/
  if((ffc=fopen("HYDR/FieldCap.masked","r"))==NULL)   {
    printf("ERROR: Unable to open %s\n","HYDR/FieldCap.masked");
    exit(0);}

  i=0;
  while (!feof(ffc)) { 
    fgets(string,MAX,ffc);
    sscanf(string,"%f %f %f %f",&lngfc[i],&latfc[i],&fc30[i],&fc100[i]);
    i++;
  }
  fclose(ffc);





  /* wilting point...*/
  if ((fwp=fopen("HYDR/WiltPont.masked","r"))==NULL)   {
    printf("ERROR: Unable to open %s\n","HYDR/WiltPont.masked");
    exit(0);
  }
  
  i=0;
  while (!feof(fwp)) { 
    fgets(string,MAX,fwp);
    sscanf(string,"%f %f %f %f",&lngwp[i],&latwp[i],&wp30[i],&wp100[i]);
    i++;
  }
  fclose(fwp);



  
  /* mean annual precip */
  if((fprcp=fopen(argv[3],"r"))==NULL)   {
    printf("ERROR: Unable to open %s\n",argv[3]);
    exit(0);
  }
  
  i=0;
  while (!feof(fprcp)) {
    fgets(string,MAX,fprcp);
    sscanf(string,"%f %f %f",&lngprcp[i],&latprcp[i],&prcp[i]);
    i++;
  }
  fclose(fprcp);


  /* calibrated parameters aggregated to half degree */
  if((fcal=fopen("bart_CAL/cal_param.30min","r"))==NULL)   {
    printf("ERROR: Unable to open %s\n","bart_CAL/cal_param.30min");
    exit(0);
  }
  
  i=0;
  while (!feof(fcal)) {
    fgets(string,MAX,fcal);
    sscanf(string,"%f %f %f %f %f %f %f %f",&lngcal[i],&latcal[i],&calbi[i],&caln1[i],&caln2[i],&calks1[i],&calks2[i],&cald2[i]);
    i++;
  }
  fclose(fcal);



  /* write each soil parameter to the output file */
  for (n=0;n<ll_cnt-1;n++) 
    {
      
      RUN_FLAG = 1;

      /* for a region not belonging to any of the below categories */

      INFIL = 0.25;
      
      CALIB_CONST_FC_WP = 1.0;
      CALIB_CONST_vG = 1.0;
      CALIB_CONST_Ks = 1.0;

      depth1     = 0.3; /* depth of first soil layer in m */
      depth2     = 0.7; /* depth of second soil layer in m */
      depth3_00  = 100;  /* depth of 3rd soil layer in mm of moisture - 
			    will be converted to m of depth using porosity */

      /*to keep soil densities reasonable - fix porosity values that are too high- from figure 6-4 in Dingman, max is
	about 55%*/
      if(ts30[n] > 0.55)
	ts30[n] = 0.55;
      if(ts100[n] > 0.55)
	ts100[n] = 0.55;

      /* convert from cm/day to mm/day */
      ks30[n]=ks30[n]*10;
      ks100[n]=ks100[n]*10;
      
      /* convert from g/cm3 to kg/m3 */
      bd30[n]=bd30[n]*1000;   
      bd100[n]=bd100[n]*1000;
      
      /* 3rd layer = 2nd layer...since no info on 3rd layer */
      ts3=ts100[n];                
      /* convert from mm of mositure to m of soil depth and divide by porosity */
      depth3=depth3_00/(ts3*1000); 
      
      /* is a percentage and hence divided by 100 - not totally sure why need to
	 divide by porosity here....*/
      fc30[n]   = fc30[n]*CRITICAL/100/ts30[n];   
      fc100[n]  = fc100[n]*CRITICAL/100/ts100[n];
      /* 3rd layer = 2nd layer, since no info on 3rd layer */
      fc3       = fc100[n];               
      
      /* is a percentage and hence divide by 100 */ 
      wp30[n]  = wp30[n]/100/ts30[n];  
      wp100[n] = wp100[n]/100/ts100[n];
      /* 3rd layer = 2nd layer, since no info on 3rd layer */
      wp3      = wp100[n];     
      

      /* density = bulk.dens/(1-porosity) */
      soildens1=bd30[n]/(1-ts30[n]);   
      soildens2=bd100[n]/(1-ts100[n]);
      /* 3rd layer = 2nd layer, since no info on 3rd layer */
      soildens3=soildens2;             


      /* adopted from the legendary bart nijssen's code...probably not worth doing all these calculations....since its just an initial moisture value */
      initmo1=(fc30[n])*ts30[n]*depth1*1000; 
      initmo2=(fc100[n])*ts100[n]*depth2*1000;
      initmo3=(fc3)*ts3*depth3*1000;

      /* bubbling pressure - needed for Brooks and Corey */
      hb30 = exp(5.3396738 + 0.1845038*cl30[n] 
		 - 2.48394546*ts30[n] 
		 - 0.00213853*pow(cl30[n],2.)
		 - 0.04356349*sd30[n]*ts30[n]
		 - 0.61745089*cl30[n]*ts30[n]
		 + 0.00143598*pow(sd30[n],2.)
		 * pow(ts30[n],2.)
		 - 0.00855375*pow(cl30[n],2.)
		 * pow(ts30[n],2.)
		 - 0.00001282*pow(sd30[n],2.)*cl30[n]
		 + 0.00895359*pow(cl30[n],2.)*ts30[n]
		 - 0.00072472*pow(sd30[n],2.)*ts30[n]
		 + 0.00000540*pow(cl30[n],2.)*sd30[n]
		 + 0.50028060*pow(ts30[n],2.)*cl30[n]);

      hb100 = exp(5.3396738 + 0.1845038*cl100[n] 
		 - 2.48394546*ts100[n] 
		 - 0.00213853*pow(cl100[n],2.)
		 - 0.04356349*sd100[n]*ts100[n]
		 - 0.61745089*cl100[n]*ts100[n]
		 + 0.00143598*pow(sd100[n],2.)
		 * pow(ts100[n],2.)
		 - 0.00855375*pow(cl100[n],2.)
		 * pow(ts100[n],2.)
		 - 0.00001282*pow(sd100[n],2.)*cl100[n]
		 + 0.00895359*pow(cl100[n],2.)*ts100[n]
		 - 0.00072472*pow(sd100[n],2.)*ts100[n]
		 + 0.00000540*pow(cl100[n],2.)*sd100[n]
		 + 0.50028060*pow(ts100[n],2.)*cl100[n]);

      /* Brooks-Corey residual water content - volume fraction) */
      res_moist30 = -0.0182482+0.00087269*sd30[n]+0.00513488*cl30[n]+0.02939286*ts30[n]-0.00015395*pow(cl30[n],2.)-0.0010827*sd30[n]*ts30[n]-0.00018233*pow(cl30[n],2.)*pow(ts30[n],2.)+0.00030703*pow(cl30[n],2.)*ts30[n]-0.0023584*pow(ts30[n],2.)*cl30[n];

      res_moist100 = -0.0182482+0.00087269*sd100[n]+0.00513488*cl100[n]+0.02939286*ts100[n]-0.00015395*pow(cl100[n],2.)-0.0010827*sd100[n]*ts100[n]-0.00018233*pow(cl100[n],2.)*pow(ts100[n],2.)+0.00030703*pow(cl100[n],2.)*ts100[n]-0.0023584*pow(ts100[n],2.)*cl100[n];


      /*damping depth calculation*/
      BD = bd100[n]/1000.0; //express in t/m^3
      //DP = 1.0+2.5*BD/(BD+exp(6.53-5.63*BD));  
      //I have no idea if this representation is ok or not...constant better?
      DP = 4.0;  

      /* print all the 53 columns */
      
      fprintf(fout,"%d \t",RUN_FLAG);                   /* 1 */
      fprintf(fout,"%d \t",cellid[n]);                  /* 2 */
      fprintf(fout,"%.3f \t",lat[n]);                   /* 3 */
      fprintf(fout,"%.3f \t",lng[n]);                   /* 4 */
      fprintf(fout,"%.2f \t",calbi[n]);                    /* 5 */
      fprintf(fout,"%.4f \t",d1arno[n]);                /* 6 */
      fprintf(fout,"%.1f \t",d2arno[n]);                /* 7 */
      fprintf(fout,"%.4f \t",d3arno[n]);                /* 8 */
      fprintf(fout,"%.3f \t",d4arno[n]);                /* 9 */
      fprintf(fout,"%.2f \t",watn30[n]);                /* 10 */
      fprintf(fout,"%.2f \t",watn100[n]);              /* 11 */
      fprintf(fout,"%.2f \t",watn100[n]);               /* 12 */
      fprintf(fout,"%.2f \t",ks30[n]);                  /* 13 */
      fprintf(fout,"%.2f \t",ks100[n]);                  /* 14 */
      fprintf(fout,"%.2f \t",ks100[n]);                 /* 15 */
      fprintf(fout,"%.0f \t",PHI);                      /* 16 */
      fprintf(fout,"%.0f \t",PHI);                      /* 17 */
      fprintf(fout,"%.0f \t",PHI);                      /* 18 */
      fprintf(fout,"%.2f \t",initmo1);                  /* 19 */
      fprintf(fout,"%.2f \t",initmo2);                  /* 20 */
      fprintf(fout,"%.2f \t",initmo3);                  /* 21 */
      fprintf(fout,"%.1f \t",elev[n]);                  /* 22 */
      fprintf(fout,"%.6f \t",depth1);                   /* 23 */
      fprintf(fout,"%.6f \t",cald2[n]);                   /* 24 */
      fprintf(fout,"%.6f \t",depth3);                   /* 25 */
      fprintf(fout,"%.0f \t",AVG_T);                    /* 26 */
      fprintf(fout,"%.2f \t",DP);                       /* 27 */
      fprintf(fout,"%.2f \t",hb30);                   /* 28 */
      fprintf(fout,"%.2f \t",hb100);                   /* 29 */
      fprintf(fout,"%.2f \t",hb100);                   /* 30 */
      fprintf(fout,"%.3f \t",sd30[n]/100.0);                   /* 31 */
      fprintf(fout,"%.3f \t",sd100[n]/100.0);                   /* 32 */
      fprintf(fout,"%.3f \t",sd100[n]/100.0);                   /* 33 */
      fprintf(fout,"%.1f \t",bd30[n]);                  /* 34 */
      fprintf(fout,"%.1f \t",bd100[n]);                 /* 35 */
      fprintf(fout,"%.1f \t",bd100[n]);                 /* 36 */
      fprintf(fout,"%.2f \t",soildens1);                /* 37 */
      fprintf(fout,"%.2f \t",soildens2);                /* 38 */
      fprintf(fout,"%.2f \t",soildens3);                /* 39 */
      fprintf(fout,"%.1f \t",lng[n]/15);                /* 40 */
      fprintf(fout,"%.3f \t",fc30[n]);                  /* 41 */
      fprintf(fout,"%.3f \t",fc100[n]);/* 42 */
      fprintf(fout,"%.3f \t",fc3);                      /* 43 */
      fprintf(fout,"%.3f \t",wp30[n]);                  /* 44 */
      fprintf(fout,"%.3f \t",wp100[n]);/* 45 */
      fprintf(fout,"%.3f \t",wp3);                      /* 46 */
      fprintf(fout,"%.4f \t",SOIL_ROUGH);               /* 47 */
      fprintf(fout,"%.4f \t",SNOW_ROUGH);               /* 48 */
      fprintf(fout,"%.2f \t",prcp[n]);                  /* 49 */
      fprintf(fout,"%.1f \t",res_moist30);                      /* 50 */
      fprintf(fout,"%.1f \t",res_moist100);                      /* 51 */
      fprintf(fout,"%.1f \t",res_moist100);                      /* 52 */
      fprintf(fout,"%d \n",FS_ACT);                     /* 53 */
      /*fprintf(fout,"%.3f \t",ts30[n]);*/  /*temporary*/
      /*fprintf(fout,"%.3f \t",ts100[n]);*/   /*temporary*/
      /*fprintf(fout,"%.3f \n",ts3);*/    /*temporary*/
    }
      
  
  fclose(fout);
  
  
  return 0;
  

}
