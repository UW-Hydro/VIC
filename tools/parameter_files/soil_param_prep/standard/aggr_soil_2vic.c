/* This programming is used to aggregate soil data from arc30 resolution into 1/8 vic grid cell
 
 <VICcell>  column longitude latitude file for whole Mexico basin in vic cell resolution
 <arc30cell_vege_class> column longitude latitude soil_data for each arc30 cell

*/
   
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define MAXVICC 18641 //hardwared for the VIC grid cells

int main(int argc, void *argv[])
{  
int i,j,vcells,n,count,junk;
float data;
FILE *fpvic, *fpsoil,*fpout;
float VIClon[MAXVICC], VIClat[MAXVICC], vlon, vlat;
float aggrVIC[MAXVICC]; // aggregated soil data in each VIC cell
int cellnum[MAXVICC];   // numbers of each positive 1km cell in each VIC cell

if(argc!=4) {
   printf("Usage: %s <VICcell> <arc30cell_soil_data> <output_file> \n",argv[0]);
   exit(0);
 }

if((fpvic=fopen(argv[1],"r"))==NULL) {
   printf("ERROR: Unable to open %s\n",argv[1]);
   exit(0);    }

if((fpsoil=fopen(argv[2],"r"))==NULL) {
   printf("ERROR: Unable to open %s\n",argv[2]);
   exit(0);    }

// read in VIC cell information
for(vcells=0;vcells<MAXVICC;vcells++)
   fscanf(fpvic,"%f %f", 
       &VIClat[vcells], &VIClon[vcells]);
 fclose(fpvic);
 
// initialize...
 for(i=0;i<MAXVICC;i++) {
  aggrVIC[i] = 0;
  cellnum[i] = 0;   
  };
  
  n=0;
 count = 0;

 while(fscanf(fpsoil,"%f %f %f", &vlon, &vlat, &data)==3) {

   /*printf ("%f %f %f\n",vlon ,vlat,data);*/

 n++;
 if ( n == 100000 ) {count ++; printf( "%d \n", count ); n=0;};
     

//assign to VIC cell based on vlat, vlon
/* for(i=0;i<MAXVICC;i++) */
 for(i=0;i<MAXVICC;i++) {

     if((VIClon[i]-1/12/2) <= vlon && ( VIClon[i]+1/12/2) > vlon &&
        (VIClat[i]+1/12/2) > vlat && (VIClat[i]-1/12/2) <= vlat ) {
       
       
       if (data >= 0) {aggrVIC[i] += data;
                       cellnum[i]++;}
       }
   }  /* end finding cl_mod cell loop */

 }

 if( (fpout = fopen(argv[3], "w")) == NULL ){
       printf("Can't open %s\n",argv[3]);
       exit(0);
     }


/* write this stuff out in ASCII format */
/*for(i=0;i<MAXVICC;i++)*/
for(i=0;i<MAXVICC;i++) {

  printf("number of arc30 cells for vic cell%d: %d\n", i, cellnum[i]);
 
  fprintf(fpout,"%.4f\n", aggrVIC[i]/cellnum[i] );
  
}

fclose(fpout);

}
