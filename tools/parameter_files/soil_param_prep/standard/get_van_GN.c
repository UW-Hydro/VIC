/* 
   to get the brooks-corey or van Genuchten N from the sand and clay percentages; determine lambda using eqns on pg.5.15 in handbook of hydrologyl; use the lambda to get N using eqns on pg. 5.6

*/
// modified by N. Voisin to handle the clipped sandclay files directly 
//( output from the usda.triangle that has been clipped out )
// note that the porosity must be in fraction format, NOT percent!


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define MAX 300
#define MXR 1961

int main (int argc, char *argv[]) 
{

  int i,cellid[MXR],ll_cnt;
  float lng[MXR],lat[MXR],sand10[MXR],clay10[MXR],sand30[MXR],clay30[MXR],sand100[MXR],clay100[MXR],p10lat[MXR],p10lng[MXR],p30lat[MXR],p30lng[MXR],p100lat[MXR],p100lng[MXR],p10[MXR],p30[MXR],p100[MXR],lambda10[MXR],lambda30[MXR],lambda100[MXR],n10[MXR],n30[MXR],n100[MXR];
  char junk[50];
  float lon,lt;
  
  FILE *fout,*fll,*fsc10,*fsc30,*fsc100,*fp10,*fp30,*fp100;

  /* usage */
  if(argc!=6) {
    printf("Usage: %s  SaCl10  SaCl30   SaCl100   lng lat file   out file \n",argv[0]);
    exit(0);
  }
  
  fout = fopen(argv[5],"w");
  
  /* open and read the lng lat file */
  if((fll=fopen(argv[4],"r"))==NULL)   {
    printf("ERROR: Unable to open %s\n",argv[4]);
    exit(0);
  }
 
 for(i=0;i<MXR;i++) {
   if ( fscanf(fll,"%f %f",&lat[i], &lng[i]) != 2 ) {
     fprintf(stderr,"Error reading %s\n",argv[4]);
     exit(-1);
   }
 }
 fclose(fll);
 printf ("%.4f %.4f \n", lat[1],lng[1]); 

 /* open and read the SaCl10 file */
  if ((fsc10=fopen(argv[1],"r"))==NULL)   {
    printf("ERROR: Unable to open %s\n",argv[1]);
    exit(0);
  }
 for(i=0;i<MXR;i++) {
   //if ( fscanf(fsc10,"%f %f %[^\n]",&sand10[i], &clay10[i], &junk) != 3 ) {
   if ( fscanf(fsc10,"%f %f \n",&sand10[i], &clay10[i]) != 2 ) {
     fprintf(stderr,"Error reading %s, %d, %d\n",argv[1], i, MXR);
     exit(-1);
   }
 }
 fprintf(stdout,"10 %f %f \n",sand10[10],clay10[10]);
 fclose(fsc10); 
 
 /* open and read the SaCl30 file */
 if((fsc30=fopen(argv[2],"r"))==NULL)   {
    printf("ERROR: Unable to open %s\n",argv[2]);
    exit(0);
 }
 for(i=0;i<MXR;i++){
   //if ( fscanf(fsc30,"%f %f %[^\n]\n",&sand30[i], &clay30[i], &junk) != 3 ) {
   if ( fscanf(fsc30,"%f %f\n",&sand30[i], &clay30[i]) != 2) {
    fprintf(stderr,"Error reading %s\n",argv[2]);
    exit(-1);
   }
 }
  fprintf(stdout,"10 %f %f \n",sand30[10],clay30[10]);
  fclose(fsc30); 
  
 /* open and read the SaCl100 file */
  if((fsc100=fopen(argv[3],"r"))==NULL)   {
    printf("ERROR: Unable to open %s\n",argv[3]);
    exit(0);
  }
  for(i=0;i<MXR;i++) {
    //if ( fscanf(fsc100,"%f %f %[^\n]\n",&sand100[i], &clay100[i], &junk) != 3 ) {
    if ( fscanf(fsc100,"%f %f\n",&sand100[i], &clay100[i]) != 2 ) {
    fprintf(stderr,"Error reading %s\n",argv[3]);
    exit(-1);
   }
 }
  fprintf(stdout,"10 %f %f \n",sand100[10],clay100[10]);
  fclose(fsc100); 

 /* open and read the porosity file */ 
 // NV the cliiped files for all three layers
  if ((fp10=fopen("clipped_tmp_files/porosity.txt","r"))==NULL) {
    printf("ERROR: Unable to open %s\n","ThetaS10");
    exit(0);
  }
  
  for(i=0;i<MXR;i++) {
    if ( fscanf(fp10,"%f %f %f\n", &p10[i], &p30[i], &p100[i]) != 3 ) {
    fprintf(stderr,"Error reading clipped_tmp_files/porosity.txt\n");
    exit(-1);
   }
 }
  fprintf(stdout,"%f %f %f \n",p10[10],p30[10],p100[10]);
  for(i=0;i<MXR;i++) {
     p10[i]/=100;
     p30[i]/=100;
     p100[i]/=100;
  }
  fclose(fp10);

  for(i=0;i<MXR;i++) {
    lambda10[i] = exp((-0.7842831)+(0.0177544*sand10[i])-(1.062498*p10[i])-(0.00005304*sand10[i]*sand10[i])-(0.00273493*clay10[i]*clay10[i])+(1.11134946*p10[i]*p10[i])-(0.03088295*sand10[i]*p10[i])+(0.00026587*sand10[i]*sand10[i]*p10[i]*p10[i])-(0.00610522*clay10[i]*clay10[i]*p10[i]*p10[i])-(0.00000235*sand10[i]*sand10[i]*clay10[i])+(0.00798746*clay10[i]*clay10[i]*p10[i])-(0.00674491*p10[i]*p10[i]*clay10[i]));
    
    lambda30[i] = exp((-0.7842831)+(0.0177544*sand30[i])-(1.062498*p30[i])-(0.00005304*sand30[i]*sand30[i])-(0.00273493*clay30[i]*clay30[i])+(1.11134946*p30[i]*p30[i])-(0.03088295*sand30[i]*p30[i])+(0.00026587*sand30[i]*sand30[i]*p30[i]*p30[i])-(0.00610522*clay30[i]*clay30[i]*p30[i]*p30[i])-(0.00000235*sand30[i]*sand30[i]*clay30[i])+(0.00798746*clay30[i]*clay30[i]*p30[i])-(0.00674491*p30[i]*p30[i]*clay30[i]));
    
    lambda100[i] = exp((-0.7842831)+(0.0177544*sand100[i])-(1.062498*p100[i])-(0.00005304*sand100[i]*sand100[i])-(0.00273493*clay100[i]*clay100[i])+(1.11134946*p100[i]*p100[i])-(0.03088295*sand100[i]*p100[i])+(0.00026587*sand100[i]*sand100[i]*p100[i]*p100[i])-(0.00610522*clay100[i]*clay100[i]*p100[i]*p100[i])-(0.00000235*sand100[i]*sand100[i]*clay100[i])+(0.00798746*clay100[i]*clay100[i]*p100[i])-(0.00674491*p100[i]*p100[i]*clay100[i]));

    n10[i] = 3 + (2/lambda10[i]);
    n30[i] = 3 + (2/lambda30[i]);
    n100[i] = 3 + (2/lambda100[i]);

    fprintf(fout,"%d %.4f %.4f %.3f %.2f %.3f %.2f %.3f %.2f\n",i+1,lng[i],lat[i],lambda10[i],n10[i],lambda30[i],n30[i],lambda100[i],n100[i]);

  }


  fclose(fout);

 return 0; 

}
