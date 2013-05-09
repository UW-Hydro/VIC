#include <stdio.h>
#include <stdlib.h>
#include <math.h>
// Nathalie Voisin

#ifndef FALSE
#define FALSE 0
#define TRUE !FALSE
#endif

#define BOXCELL 4480  // number of cells in the box files from the soil program
#define CLIPCELL 1961 // number of cells in the clipped file

main( int argc, char** argv)
{
float *lat,*lon;
float llt, lln,value[3], junk;
int   i,j, istart;
FILE *fp,*fp2,*fp3,*fp4;
char  nomatch;
nomatch = TRUE;


if(argc != 3){
  printf("Wrong command line arguments: enter <inputfile> <outputfile>  .\n");
  exit (0);
} 

lat = (float*)calloc(CLIPCELL,sizeof(float));
lon = (float*)calloc(CLIPCELL,sizeof(float));
if ( lat == NULL || lon == NULL ) {
  fprintf(stderr,"Error allocaing\n");
  exit(-1);
}

// read clipped soil file
 fp3=fopen("zb_025.ll", "r");
 if (fp3==NULL){ 
    fprintf(stderr, "NULL   \n");
    exit(3);
 }
 for (j=0;j<CLIPCELL;j++){
   if (fscanf(fp3,"%f %f\n",&lat[j],&lon[j]) != 2 ) {
      fprintf(stderr, "error reading listcell\n");
      exit(-1);
   }
 }
 fclose(fp3);


//open input file and output file
 fp=fopen(argv[1], "r");
 fprintf(stderr,"opening inputfile\n");
 if (fp==NULL) {fprintf(stderr,"NULL  %s \n", argv[1]); exit(1);}

 fp2=fopen(argv[2], "w");
  if (fp2==NULL){ fprintf(stderr, "NULL   %s \n", argv[2]);exit(2);}

 fp4=fopen("zb_box_latlon","r");
 if (fp4==NULL){ fprintf(stderr, "NULL   lonlat \n");exit(4);}

//rewrite the file 
 for (j=0;j<BOXCELL;j++){
   if ( fscanf(fp4,"%f %f\n", &lln, &llt) != 2 ) {
     fprintf(stderr,"Error readind zb_box_latlon\n");
     exit(-1);
   }
   if (fscanf(fp,"%f %f %f\n", &value[0],&value[1],&value[2]) != 3 ) {
     fprintf(stderr, "error reading %s at %d\n", argv[1], j);
     exit(3);
   }
   nomatch = TRUE;  
   i=-1;
   while ( i<CLIPCELL &&  nomatch ){
     i++;
     if ( fabs(llt-lat[i])<0.001 && fabs(lln-lon[i])<0.001 ){
       fprintf(fp2," %.4f %.4f %.4f \n",value[0],value[1],value[2]);
       nomatch=FALSE;
     }
   }
 }

fclose(fp2);
fclose(fp);
fclose(fp4);

free(lat);
free(lon);
}   

