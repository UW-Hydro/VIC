/* This program rescales the gridded data from the program regrid with       
   the ratio of  precipvalue in cell * (prism.jan/monthly.jan). This is done 
   to achieve monthly precipitation in each grid cell equal to prism monthly 
   means. Inputfiles are the gridded data file, monthly values (mk_monthly)  
   and prism monthly values (get_prism).  Output is a new rescaled gridded   
   file that can be turned into the final VIC input files.                
   Last Modified 5-18-00:  I/O streamlined to remove tens of thousands of 
                           repetitive file reads                        */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define MAXLINE 80

int timesteps(int yy, int mt);
float *get_1d_mem_float(int dim1);
float **get_2d_mem_float(int dim1, int dim2);

void main(int argc, char **argv) 
{ 
  FILE *fpmask, *fpmonthly, *fpprism,*fpin_grd, *fpout_grd;        
  char maskfile[400];       // Total length of directory and filename < 400
  char monthly[400],old_monthly[400];  // Use files monthly.jan, .feb, etc. 
  char prism[400],old_prism[400];   // Using files prism.jan, prism.feb...
  char in_grd[400];
  char out_grd[400];
  char start_year[400],end_year[400],st_mon[400],end_mon[400];
  char *monthstr[]={"jan","feb","mar","apr","may","jun","jul",
                    "aug","sep","oct","nov","dec"};
  char str1[MAXLINE];
  float value, voidval;
  int year,endyy,st_mo,end_mo;
  int i, j, k, l, q;   
  int binflag,n_cells, active;
  float mm, pp, grd, mult, high, low;
  short int grd_int, value_int;
  int rows, cols;
  float *maskval;   // mask values (read once) 
  float **scalar;  // holds scalars for all cells and all months

  if (argc!=12) {    // Must be exactly 11 arguments behind the program name 
    printf("Not correct number of commandline arguments \n");
    exit(EXIT_FAILURE); 
  }
  strcpy(maskfile,argv[1]);      printf("%s \n",maskfile);
  strcpy(monthly,argv[2]);       printf("%s \n",monthly);
  strcpy(prism,argv[3]);         printf("%s \n",prism);
  strcpy(in_grd,argv[4]);        printf("%s \n",in_grd);
  strcpy(out_grd,argv[5]);       printf("%s \n",out_grd);
  strcpy(start_year,argv[6]);year=atoi(start_year);  printf("%d \n",year);
  // shift start and end months to be on 0-11 instead of 1-12 
  strcpy(st_mon,argv[7]); st_mo=atoi(st_mon)-1; printf("%d \n",st_mo+1);
  strcpy(end_year,argv[8]); endyy=atoi(end_year); printf("%d \n",endyy);
  strcpy(end_mon,argv[9]); end_mo=atoi(end_mon)-1; printf("%d \n",end_mo+1);
  binflag=atoi(argv[10]); printf("%d \n",binflag);
  mult=atof(argv[11]); printf("%.1f \n",mult);

  strcpy(old_monthly,monthly); 
  strcpy(old_prism,prism);  

  //open input and output files depending on binary status 
  if(binflag!=1) {
    if((fpin_grd = fopen(in_grd,"r"))==NULL){
      printf("Cannot open file %s \n",in_grd);exit(0);}
    if((fpout_grd = fopen(out_grd,"w"))==NULL){
      printf("Cannot open file %s \n",out_grd);exit(0);}
  } else {
    if((fpin_grd = fopen(in_grd,"rb"))==NULL){
      printf("Cannot open file %s \n",in_grd);exit(0);}
    if((fpout_grd = fopen(out_grd,"wb"))==NULL){
      printf("Cannot open file %s \n",out_grd);exit(0);}
  }

  // open mask file 
  if((fpmask = fopen(maskfile,"r"))==NULL){                  
    printf("Cannot open file %s \n",maskfile);
    exit(0);
  }
  // process mask file 
  fscanf(fpmask,"%*s %d",&cols); // read header
  fscanf(fpmask,"%*s %d",&rows);
  for(i=0;i<4;i++) 
    fscanf(fpmask,"%*s %s", &str1); //!ARC-style
  voidval = atof(str1);
  printf("void = %.2f\n",voidval);
  high = voidval + 0.001;
  low  = voidval - 0.001;
 
  // calloc mask array
  printf("rows=%d, cols=%d\n", rows, cols);
  maskval = get_1d_mem_float(rows*cols);

  n_cells = 0;   
  active = 0;
  while(fscanf(fpmask,"%s",&str1)!=EOF) {  // read body
    maskval[n_cells] = atof(str1);
    if(maskval[n_cells] < low || maskval[n_cells] > high) active++;
    n_cells++;
  }
  printf("total cells in mask = %d\n",n_cells);
  printf("active cells in mask = %d\n",active);
  fclose(fpmask);  
  if(n_cells == active) {     // check for incorrect void discrimination
    fprintf(stderr, "total mask cell number equals active cell number\n");
    fprintf(stderr, "this could indicate an incorrect void value\n");
  }

  // calloc scalar array
  scalar = get_2d_mem_float(rows*cols, 12);

  // get array of scalars for each month and grid cell
  for (k=0;k<12;k++) {  // 0=january, 1=february, ...
    // set filenames and open monthly.mnth and prism.mnth each month 
    strcpy(monthly,old_monthly);
    strcat(monthly,monthstr[k]);
    strcpy(prism,old_prism);
    strcat(prism,monthstr[k]);       
    if((fpmonthly = fopen(monthly,"r"))==NULL){
      printf("Cannot open file %s \n",monthly);exit(0);}
    if((fpprism = fopen(prism,"r"))==NULL){
      printf("Cannot open file %s \n",prism);exit(0);}

    for(i=0;i<6;i++){ //read over 6-line headers
      fgets(str1,MAXLINE,fpmonthly);
      fgets(str1,MAXLINE,fpprism);
    }
    for(q=0;q<n_cells;q++) {   //loop thru all mask cells 
      fscanf(fpmonthly,"%f",&mm); 
      fscanf(fpprism,"%f",&pp);
      if(maskval[q] < low || maskval[q] > high)  { //valid data
        if (pp<0.01)
          pp = mm; // reassign missing prism values 
        scalar[q][k] = pp/mm;
      } else {
        scalar[q][k] = -999.;  // void doesn't matter really
      }
    }
    fclose(fpmonthly);
    fclose(fpprism);         
  }

  // ------- START LOOPING THROUGH SPACE & TIME ---------------------
  for(l=year;l<=endyy;l++) {    /* Loop for each year */
    for (k=0;k<12;k++) {  // month number 0 is january and 1 is february 

      if(( l>year || k>=st_mo) && (l<endyy || k<=end_mo)) {
        printf("year %d month %d timesteps %d \n",l,k+1,timesteps(l,k));
   
	for(j=0;j<timesteps(l,k);j++) {// Timestep loop: days per mon --

          for(q=0;q<n_cells;q++) {   //loop thru all mask cells 
            if(maskval[q]<low || maskval[q]>high)  { //valid data
              if(binflag!=1) {       //reading ascii data
                fscanf(fpin_grd,"%f",&grd); //scan value from inputfile
                value = grd*scalar[q][k]; // rescale value 
                fprintf(fpout_grd,"%6.1f ",value);
              } else {
		fread(&grd_int,2,1,fpin_grd);
		// CONVERT BACK TO REAL PRECIP VALUE
		grd = ((float) grd_int)/mult;
		if(grd<0.0)printf("neg. prcp in grd file, cell %i yr %i mo %i\n",
				  q+1,l,k+1);
		if(scalar[q][k]<0)printf("negative scalar, cell %i mo %i\n",
					 q+1,k+1);
                value = grd*scalar[q][k];
		if(value<0.0)printf("value negative, cell %i mo %i\n",
					 q+1,k+1);
                value_int = (short int)(value*mult);
		fwrite(&value_int,sizeof(short int),1,fpout_grd);
	      }
            } // active cell condition end
          }   // mask loop done
          if(binflag!=1) fprintf(fpout_grd,"\n");
	}  //timestep loop done ----------------------------------------
      } // month condition end 
    } // month loop done
  } //year loop done
 
  fclose(fpout_grd);
  fclose(fpin_grd);
} /** END MAIN  ***/

// ------- SUBROUTINES ------------------------------------------------


int timesteps(int yy, int mt)
{
  int val;
  char *regular[]={"31","28","31","30","31","30","31","31",
                   "30","31","30","31"};
  char *leap[]=   {"31","29","31","30","31","30","31","31",
                   "30","31","30","31"};
  if (yy%4!=0) val = atoi(regular[mt%12]);
  if (yy%4==0) val = atoi(leap[mt%12]);
  return val;
} 


float *get_1d_mem_float(int dim1) 
{
  // allocate memory for a 2-d array, pass back pointer
  float *array_1d;

  if (!(array_1d = (float *) calloc(dim1, sizeof(float)))) {
    fprintf(stderr, "unable to allocate memory in get_1d_mem_float\n");
    fprintf(stderr, "dimension %d\n", dim1);
    exit(1);   
  }
  return array_1d;
}

float **get_2d_mem_float(int dim1, int dim2) 
{
  // allocate memory for a 2-d array, pass back pointer
  float **array_2d;
  int i;

  if (!(array_2d = (float **) calloc(dim1, sizeof(float *)))) {
    fprintf(stderr, "unable to allocate memory in get_2d_mem_float\n");
    fprintf(stderr, "dimension1= %d\n", dim1);
    exit(1);   
  }
  for (i = 0; i < dim1; i++) {
    if (!(array_2d[i] = (float *) calloc(dim2, sizeof(float)))) {
      fprintf(stderr, "unable to allocate memory in get_2d_mem_float\n");
      fprintf(stderr, "dimension1= %d, dimension2=%d\n", dim1,dim2);
      exit(1);
    }
  } 
  return array_2d;
}




