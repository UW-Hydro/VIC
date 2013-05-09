#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX_LAYERS 3
#define FALSE 0
#define TRUE !FALSE

typedef struct {
  char     FS_ACTIVE;                 /* if TRUE frozen soil algorithm is 
					 active in current grid cell */
  float   Ds;                        /* fraction of maximum subsurface flow 
					 rate */
  float   Dsmax;                     /* maximum subsurface flow rate 
					 (mm/day) */
  float   Ksat[MAX_LAYERS];          /* saturated hydraulic  conductivity 
					 (mm/day) */
  float   Wcr_FRACT[MAX_LAYERS];           /* critical moisture level for soil 
					 layer, evaporation is no longer 
					 affected moisture stress in the 
					 soil (mm) */
  float   Wpwp_FRACT[MAX_LAYERS];          /* soil moisture content at permanent 
					 wilting point (mm) */
  float   Ws;                        /* fraction of maximum soil moisture */
  float   annual_prec;               /* annual average precipitation (mm) */
  float   avg_temp;                  /* average soil temperature (C) */
  float   b_infilt;                  /* infiltration parameter */
  float   bubble[MAX_LAYERS];        /* bubbling pressure, HBH 5.15 (cm) */
  float   bulk_density[MAX_LAYERS];  /* soil bulk density (kg/m^3) */
  float   c;                         /* exponent */
  float   clay[MAX_LAYERS];
  float   depth[MAX_LAYERS];         /* thickness of each soil moisture 
					 layer (m) */
  float   dp;                        /* soil thermal damping depth (m) */
  float   expt[MAX_LAYERS];          /* pore-size distribution per layer, 
					 HBH 5.15 */
  float   init_moist[MAX_LAYERS];    /* initial layer moisture level (mm) */
  float   max_moist[MAX_LAYERS];     /* maximum moisture content (mm) per 
					 layer */
  float   off_gmt;
  float   phi_s[MAX_LAYERS];         /* soil moisture diffusion parameter 
					 (mm/mm) */
  float   porosity[MAX_LAYERS];      /* porosity (fraction) */
  float   quartz[MAX_LAYERS];        /* quartz content of soil (fraction) */
  float   resid_moist[MAX_LAYERS];   /* residual moisture content of soil 
					 layer */
  float   rough;                     /* soil surface roughness (m) */
  float   sand[MAX_LAYERS];
  float   snow_rough;                /* snow surface roughness (m) */
  float   soil_density[MAX_LAYERS];  /* soil partical density (kg/m^3) */
  float    elevation;                 /* grid cell elevation (m) */
  float    lat;                       /* grid cell central latitude */
  float    lng;                       /* grid cell central longitude */
  int      gridcel;                   /* grid cell number */
} soil_con_struct;

void read_arcinfo_mask(FILE *, int *, float **, float **, int **, int **);
void get_arc_info_values(FILE *, int, float *, float *, float *);

int main(int argc, char *argv[]) {
/****************************************************************
  convert_ARCINFO_to_ASCII_soil_file.c 
                          Keith Cherkauer          May 8, 2000

  This program was written to read a runfile mask and an ARC/INFO
  soil parameter control file and create an ASCII column soil
  data file for the grid cells activated in the mask file.
****************************************************************/

  FILE            *fmask, *farc, *fascii, *flist;
  char             arcdir[512], arcname[512], tmpstr[512];
  int              i, j;
  float           *values;
  float           *cell_lat;
  float           *cell_lng;
  int             *cell_row;
  int             *cell_col;
  int              Ncells;
  int              Nlayers;
  soil_con_struct *soil_con;
  
  if(argc!=6) {
    fprintf(stderr,"Usage: %s <basin mask> <ARC/INFO file list> <Arc/INFO file dir> <ASCII param file> <number of layers>\n",argv[0]);
    fprintf(stderr,"\tThis program was reads soil parameters from ARC/INFO grid cell parameter files and creates an ASCII column file.  A runfile mask is used to limit output cells to those active in the basin.  Activated cells have positive values in the run mask file.\n");
    exit(0);
  }
  
  if((fmask=fopen(argv[1],"r"))==NULL) {
    fprintf(stderr,"ERROR: Unable to open mask file %s\n",argv[1]);
    exit(0);
  }
  
  /** Read masked cells, need lat, long, row, col, total number **/
  read_arcinfo_mask(fmask, &Ncells, &cell_lat, &cell_lng, &cell_row, 
		    &cell_col);
  fclose(fmask);

  
  /** Allocate soil_con_struct for the number of grid cells present **/
  soil_con = (soil_con_struct *)calloc(Ncells, sizeof(soil_con_struct));
  values   = (float *)calloc(Ncells, sizeof(float));
  
  /** Open ARC/INFO file list **/
  if((flist=fopen(argv[2],"r"))==NULL) {
    fprintf(stderr,"ERROR: Unable to open ARC/INFO file list %s\n",argv[2]);
    exit(0);
  }
  
  /** Store ARC/INFO directory **/
  sprintf(arcdir,"%s",argv[3]);
  
  /** Open ACSII output file **/
  if((fascii=fopen(argv[4],"w"))==NULL) {
    fprintf(stderr,"ERROR: Unable to open %s\n",argv[4]);
    exit(0);
  }
  
  /** Set number of layers **/
  Nlayers = atoi(argv[5]);
  
  /** Read cell number **/
  fscanf(flist,"%s",tmpstr);
  sprintf(arcname, "%s/%s", arcdir, tmpstr);
  if((farc=fopen(arcname,"r"))==NULL) {
    fprintf(stderr,"ERROR: Unable to open %s\n",arcname);
    exit(0);
  }
  get_arc_info_values(farc, Ncells, cell_lat, cell_lng, values);
  for ( i = 0; i < Ncells; i++ ) soil_con[i].gridcel = (int)values[i];
  fclose(farc);
  
  /** Skip run mask **/
  fscanf(flist,"%s",tmpstr);

  /** Read elevation **/
  fscanf(flist,"%s",tmpstr);
  sprintf(arcname, "%s/%s", arcdir, tmpstr);
  if((farc=fopen(arcname,"r"))==NULL) {
    fprintf(stderr,"ERROR: Unable to open %s\n",arcname);
    exit(0);
  }
  get_arc_info_values(farc, Ncells, cell_lat, cell_lng, values);
  for ( i = 0; i < Ncells; i++ ) soil_con[i].elevation = values[i];
  fclose(farc);

  /** Read b_infilt **/
  fscanf(flist,"%s",tmpstr);
  sprintf(arcname, "%s/%s", arcdir, tmpstr);
  if((farc=fopen(arcname,"r"))==NULL) {
    fprintf(stderr,"ERROR: Unable to open %s\n",arcname);
    exit(0);
  }
  get_arc_info_values(farc, Ncells, cell_lat, cell_lng, values);
  for ( i = 0; i < Ncells; i++ ) soil_con[i].b_infilt = values[i];
  fclose(farc);

  /** Read Ds **/
  fscanf(flist,"%s",tmpstr);
  sprintf(arcname, "%s/%s", arcdir, tmpstr);
  if((farc=fopen(arcname,"r"))==NULL) {
    fprintf(stderr,"ERROR: Unable to open %s\n",arcname);
    exit(0);
  }
  get_arc_info_values(farc, Ncells, cell_lat, cell_lng, values);
  for ( i = 0; i < Ncells; i++ ) soil_con[i].Ds = values[i];
  fclose(farc);

  /** Read Dsmax **/
  fscanf(flist,"%s",tmpstr);
  sprintf(arcname, "%s/%s", arcdir, tmpstr);
  if((farc=fopen(arcname,"r"))==NULL) {
    fprintf(stderr,"ERROR: Unable to open %s\n",arcname);
    exit(0);
  }
  get_arc_info_values(farc, Ncells, cell_lat, cell_lng, values);
  for ( i = 0; i < Ncells; i++ ) soil_con[i].Dsmax = values[i];
  fclose(farc);

  /** Read Ws **/
  fscanf(flist,"%s",tmpstr);
  sprintf(arcname, "%s/%s", arcdir, tmpstr);
  if((farc=fopen(arcname,"r"))==NULL) {
    fprintf(stderr,"ERROR: Unable to open %s\n",arcname);
    exit(0);
  }
  get_arc_info_values(farc, Ncells, cell_lat, cell_lng, values);
  for ( i = 0; i < Ncells; i++ ) soil_con[i].Ws = values[i];
  fclose(farc);

  /** Read c **/
  fscanf(flist,"%s",tmpstr);
  sprintf(arcname, "%s/%s", arcdir, tmpstr);
  if((farc=fopen(arcname,"r"))==NULL) {
    fprintf(stderr,"ERROR: Unable to open %s\n",arcname);
    exit(0);
  }
  get_arc_info_values(farc, Ncells, cell_lat, cell_lng, values);
  for ( i = 0; i < Ncells; i++ ) soil_con[i].c = values[i];
  fclose(farc);

  /** Read avg_temp **/
  fscanf(flist,"%s",tmpstr);
  sprintf(arcname, "%s/%s", arcdir, tmpstr);
  if((farc=fopen(arcname,"r"))==NULL) {
    fprintf(stderr,"ERROR: Unable to open %s\n",arcname);
    exit(0);
  }
  get_arc_info_values(farc, Ncells, cell_lat, cell_lng, values);
  for ( i = 0; i < Ncells; i++ ) soil_con[i].avg_temp = values[i];
  fclose(farc);

  /** Read dp **/
  fscanf(flist,"%s",tmpstr);
  sprintf(arcname, "%s/%s", arcdir, tmpstr);
  if((farc=fopen(arcname,"r"))==NULL) {
    fprintf(stderr,"ERROR: Unable to open %s\n",arcname);
    exit(0);
  }
  get_arc_info_values(farc, Ncells, cell_lat, cell_lng, values);
  for ( i = 0; i < Ncells; i++ ) soil_con[i].dp = values[i];
  fclose(farc);

  /** Read off_gmt **/
  fscanf(flist,"%s",tmpstr);
  sprintf(arcname, "%s/%s", arcdir, tmpstr);
  if((farc=fopen(arcname,"r"))==NULL) {
    fprintf(stderr,"ERROR: Unable to open %s\n",arcname);
    exit(0);
  }
  get_arc_info_values(farc, Ncells, cell_lat, cell_lng, values);
  for ( i = 0; i < Ncells; i++ ) soil_con[i].off_gmt = values[i];
  fclose(farc);

  /** Read Wcr **/
  for ( j = 0; j < Nlayers ; j++ ) {
    fscanf(flist,"%s",tmpstr);
    sprintf(arcname, "%s/%s", arcdir, tmpstr);
    if((farc=fopen(arcname,"r"))==NULL) {
      fprintf(stderr,"ERROR: Unable to open %s\n",arcname);
      exit(0);
    }
    get_arc_info_values(farc, Ncells, cell_lat, cell_lng, values);
    for ( i = 0; i < Ncells; i++ ) soil_con[i].Wcr_FRACT[j] = values[i];
    fclose(farc);
  }

  /** Read Wpwp **/
  for ( j = 0; j < Nlayers ; j++ ) {
    fscanf(flist,"%s",tmpstr);
    sprintf(arcname, "%s/%s", arcdir, tmpstr);
    if((farc=fopen(arcname,"r"))==NULL) {
      fprintf(stderr,"ERROR: Unable to open %s\n",arcname);
      exit(0);
    }
    get_arc_info_values(farc, Ncells, cell_lat, cell_lng, values);
    for ( i = 0; i < Ncells; i++ ) soil_con[i].Wpwp_FRACT[j] = values[i];
    fclose(farc);
  }

  /** Read rough **/
  fscanf(flist,"%s",tmpstr);
  sprintf(arcname, "%s/%s", arcdir, tmpstr);
  if((farc=fopen(arcname,"r"))==NULL) {
    fprintf(stderr,"ERROR: Unable to open %s\n",arcname);
    exit(0);
  }
  get_arc_info_values(farc, Ncells, cell_lat, cell_lng, values);
  for ( i = 0; i < Ncells; i++ ) soil_con[i].rough = values[i];
  fclose(farc);

  /** Read snow_rough **/
  fscanf(flist,"%s",tmpstr);
  sprintf(arcname, "%s/%s", arcdir, tmpstr);
  if((farc=fopen(arcname,"r"))==NULL) {
    fprintf(stderr,"ERROR: Unable to open %s\n",arcname);
    exit(0);
  }
  get_arc_info_values(farc, Ncells, cell_lat, cell_lng, values);
  for ( i = 0; i < Ncells; i++ ) soil_con[i].snow_rough = values[i];
  fclose(farc);

  /** Read annual_prec **/
  fscanf(flist,"%s",tmpstr);
  sprintf(arcname, "%s/%s", arcdir, tmpstr);
  if((farc=fopen(arcname,"r"))==NULL) {
    fprintf(stderr,"ERROR: Unable to open %s\n",arcname);
    exit(0);
  }
  get_arc_info_values(farc, Ncells, cell_lat, cell_lng, values);
  for ( i = 0; i < Ncells; i++ ) soil_con[i].annual_prec = values[i];
  fclose(farc);

  /** Read sand **/
  for ( j = 0; j < Nlayers ; j++ ) {
    fscanf(flist,"%s",tmpstr);
    sprintf(arcname, "%s/%s", arcdir, tmpstr);
    if((farc=fopen(arcname,"r"))==NULL) {
      fprintf(stderr,"ERROR: Unable to open %s\n",arcname);
      exit(0);
    }
    get_arc_info_values(farc, Ncells, cell_lat, cell_lng, values);
    for ( i = 0; i < Ncells; i++ ) soil_con[i].sand[j] = values[i];
    fclose(farc);
  }

  /** Read clay **/
  for ( j = 0; j < Nlayers ; j++ ) {
    fscanf(flist,"%s",tmpstr);
    sprintf(arcname, "%s/%s", arcdir, tmpstr);
    if((farc=fopen(arcname,"r"))==NULL) {
      fprintf(stderr,"ERROR: Unable to open %s\n",arcname);
      exit(0);
    }
    get_arc_info_values(farc, Ncells, cell_lat, cell_lng, values);
    for ( i = 0; i < Ncells; i++ ) soil_con[i].clay[j] = values[i];
    fclose(farc);
  }

  /** Read Ksat **/
  for ( j = 0; j < Nlayers ; j++ ) {
    fscanf(flist,"%s",tmpstr);
    sprintf(arcname, "%s/%s", arcdir, tmpstr);
    if((farc=fopen(arcname,"r"))==NULL) {
      fprintf(stderr,"ERROR: Unable to open %s\n",arcname);
      exit(0);
    }
    get_arc_info_values(farc, Ncells, cell_lat, cell_lng, values);
    for ( i = 0; i < Ncells; i++ ) soil_con[i].Ksat[j] = values[i];
    fclose(farc);
  }

  /** Read phi_s **/
  for ( j = 0; j < Nlayers ; j++ ) {
    fscanf(flist,"%s",tmpstr);
    sprintf(arcname, "%s/%s", arcdir, tmpstr);
    if((farc=fopen(arcname,"r"))==NULL) {
      fprintf(stderr,"ERROR: Unable to open %s\n",arcname);
      exit(0);
    }
    get_arc_info_values(farc, Ncells, cell_lat, cell_lng, values);
    for ( i = 0; i < Ncells; i++ ) soil_con[i].phi_s[j] = values[i];
    fclose(farc);
  }

  /** Read init_moist **/
  for ( j = 0; j < Nlayers ; j++ ) {
    fscanf(flist,"%s",tmpstr);
    sprintf(arcname, "%s/%s", arcdir, tmpstr);
    if((farc=fopen(arcname,"r"))==NULL) {
      fprintf(stderr,"ERROR: Unable to open %s\n",arcname);
      exit(0);
    }
    get_arc_info_values(farc, Ncells, cell_lat, cell_lng, values);
    for ( i = 0; i < Ncells; i++ ) soil_con[i].init_moist[j] = values[i];
    fclose(farc);
  }

  /** Read depth **/
  for ( j = 0; j < Nlayers ; j++ ) {
    fscanf(flist,"%s",tmpstr);
    sprintf(arcname, "%s/%s", arcdir, tmpstr);
    if((farc=fopen(arcname,"r"))==NULL) {
      fprintf(stderr,"ERROR: Unable to open %s\n",arcname);
      exit(0);
    }
    get_arc_info_values(farc, Ncells, cell_lat, cell_lng, values);
    for ( i = 0; i < Ncells; i++ ) soil_con[i].depth[j] = values[i];
    fclose(farc);
  }

  /** Read bulk_density **/
  for ( j = 0; j < Nlayers ; j++ ) {
    fscanf(flist,"%s",tmpstr);
    sprintf(arcname, "%s/%s", arcdir, tmpstr);
    if((farc=fopen(arcname,"r"))==NULL) {
      fprintf(stderr,"ERROR: Unable to open %s\n",arcname);
      exit(0);
    }
    get_arc_info_values(farc, Ncells, cell_lat, cell_lng, values);
    for ( i = 0; i < Ncells; i++ ) soil_con[i].bulk_density[j] = values[i];
    fclose(farc);
  }

  /** Read porosity **/
  for ( j = 0; j < Nlayers ; j++ ) {
    fscanf(flist,"%s",tmpstr);
    sprintf(arcname, "%s/%s", arcdir, tmpstr);
    if((farc=fopen(arcname,"r"))==NULL) {
      fprintf(stderr,"ERROR: Unable to open %s\n",arcname);
      exit(0);
    }
    get_arc_info_values(farc, Ncells, cell_lat, cell_lng, values);
    for ( i = 0; i < Ncells; i++ ) soil_con[i].porosity[j] = values[i];
    fclose(farc);
  }

  /** Read FS_active **/
  fscanf(flist,"%s",tmpstr);
  sprintf(arcname, "%s/%s", arcdir, tmpstr);
  if((farc=fopen(arcname,"r"))==NULL) {
    fprintf(stderr,"ERROR: Unable to open %s\n",arcname);
    exit(0);
  }
  get_arc_info_values(farc, Ncells, cell_lat, cell_lng, values);
  for ( i = 0; i < Ncells; i++ ) 
    if ( values[i] == 1 ) soil_con[i].FS_ACTIVE = TRUE;
    else soil_con[i].FS_ACTIVE = FALSE;
  fclose(farc);

  /** Processes ARC/INFO data to get ASCII data **/
  for ( i = 0; i < Ncells; i++ ) {
    for(j = 0; j < Nlayers; j++) {
      soil_con[i].bulk_density[j] *= 1000.;
      soil_con[i].soil_density[j] = soil_con[i].bulk_density[j] 
	/ (1.0 - soil_con[i].porosity[j]);
      soil_con[i].quartz[j] = soil_con[i].sand[j] / 100.;
      soil_con[i].max_moist[j] = soil_con[i].depth[j] * soil_con[i].porosity[j] * 1000.;
      if(soil_con[i].init_moist[j] > soil_con[i].max_moist[j]) 
	soil_con[i].init_moist[j] = soil_con[i].max_moist[j];
      soil_con[i].bubble[j] = exp(5.3396738 + 0.1845038*soil_con[i].clay[j] 
				  - 2.48394546*soil_con[i].porosity[j] 
				  - 0.00213853*pow(soil_con[i].clay[j],2.)
				  - 0.04356349*soil_con[i].sand[j]*soil_con[i].porosity[j]
				  - 0.61745089*soil_con[i].clay[j]*soil_con[i].porosity[j]
				  + 0.00143598*pow(soil_con[i].sand[j],2.)
				  * pow(soil_con[i].porosity[j],2.)
				  - 0.00855375*pow(soil_con[i].clay[j],2.)
				  * pow(soil_con[i].porosity[j],2.)
				  - 0.00001282*pow(soil_con[i].sand[j],2.)*soil_con[i].clay[j]
				  + 0.00895359*pow(soil_con[i].clay[j],2.)*soil_con[i].porosity[j]
				  - 0.00072472*pow(soil_con[i].sand[j],2.)*soil_con[i].porosity[j]
				  + 0.00000540*pow(soil_con[i].clay[j],2.)*soil_con[i].sand[j]
				  + 0.50028060*pow(soil_con[i].porosity[j],2.)*soil_con[i].clay[j]);
      soil_con[i].expt[j] = exp(-0.7842831 + 0.0177544*soil_con[i].sand[j] 
				- 1.062498*soil_con[i].porosity[j] 
				- 0.00005304*pow(soil_con[i].sand[j],2.)
				- 0.00273493*pow(soil_con[i].clay[j],2.)
				+ 1.11134946*pow(soil_con[i].porosity[j],2.)
				- 0.03088295*soil_con[i].sand[j]*soil_con[i].porosity[j]
				+ 0.00026587*pow(soil_con[i].sand[j],2.)
				* pow(soil_con[i].porosity[j],2.)
				- 0.00610522*pow(soil_con[i].clay[j],2.)
				* pow(soil_con[i].porosity[j],2.)
				- 0.00000235*pow(soil_con[i].sand[j],2.)*soil_con[i].clay[j]
				+ 0.00798746*pow(soil_con[i].clay[j],2.)*soil_con[i].porosity[j]
				- 0.00674491*pow(soil_con[i].porosity[j],2.)*soil_con[i].clay[j]);
      soil_con[i].expt[j] = 2. / soil_con[i].expt[j] + 3.;
      soil_con[i].resid_moist[j] = - 0.0182482 + 0.00087269 * soil_con[i].sand[j]
	+ 0.00513488 * soil_con[i].clay[j] 
	+ 0.02939286 * soil_con[i].porosity[j] 
	- 0.00015395 * pow(soil_con[i].clay[j],2.) 
	- 0.00108270 * soil_con[i].sand[j] * soil_con[i].porosity[j] 
	- 0.00018233 * pow(soil_con[i].clay[j],2.) 
	* pow(soil_con[i].porosity[j],2.) 
	+ 0.00030703 * pow(soil_con[i].clay[j],2.0) 
	* soil_con[i].porosity[j] 
	- 0.00235840 * pow(soil_con[i].porosity[j],2.) 
	* soil_con[i].clay[j];
      
      /** Check for valid values of generated parameters **/
      if(soil_con[i].bubble[j]<1.36) {
	fprintf(stderr,"WARNING: estimated bubbling pressure too low (%f), resetting to minimum value (%f).\n",soil_con[i].bubble[j],1.36);
	soil_con[i].bubble[j] = 1.36;
      }
      if(soil_con[i].bubble[j]>187.2) {
	fprintf(stderr,"WARNING: estimated bubbling pressure too high (%f), resetting to maximum value (%f).\n",soil_con[i].bubble[j],187.2);
	soil_con[i].bubble[j] = 187.2;
      }
      if(soil_con[i].expt[j] < 2. / 1.090 + 3.) {
	fprintf(stderr,"WARNING: estimated exponential (expt) too low (%f), resetting to minimum value (%f).\n", soil_con[i].expt[j], 2. / 1.090 + 3.);
	soil_con[i].expt[j] = 2. / 1.090 + 3.;
      }
      if(soil_con[i].expt[j] > 2. / 0.037 + 3.) {
	fprintf(stderr,"WARNING: estimated exponential (expt) too high (%f), resetting to maximum value (%f).\n",soil_con[i].expt[j], 2. / 0.037 + 3.);
	soil_con[i].expt[j] = 2. / 0.037 + 3.;
      }
      if(soil_con[i].resid_moist[j] < -0.038) {
	fprintf(stderr,"WARNING: estimated residual soil moisture too low (%f), resetting to minimum value (%f).\n",soil_con[i].resid_moist[j],-0.038);
	soil_con[i].resid_moist[j] = -0.038;
      }
      if(soil_con[i].resid_moist[j] > 0.205) {
	fprintf(stderr,"WARNING: estimated residual soil mositure too high (%f), resetting to maximum value (%f).\n",soil_con[i].resid_moist[j],0.205);
	soil_con[i].resid_moist[j] = 0.205;
      }
    }
  }
  
  /** Write ASCII soil parameter file **/
  fprintf(fascii,"#RUN\tGRID\tLAT\tLNG\tINFILT\tDs\tDs_MAX\tWs\tC\t");
  for ( j = 0; j < Nlayers; j++ ) fprintf(fascii,"EXPT_%i\t",j);  
  for ( j = 0; j < Nlayers; j++ ) fprintf(fascii,"Ksat_%i\t",j);  
  for ( j = 0; j < Nlayers; j++ ) fprintf(fascii,"PHI_%i\t",j);  
  for ( j = 0; j < Nlayers; j++ ) fprintf(fascii,"MOIST_%i\t",j);  
  fprintf(fascii,"ELEV\t");
  for ( j = 0; j < Nlayers; j++ ) fprintf(fascii,"DEPTH_%i\t",j);  
  fprintf(fascii,"AVG_T\tDP\t");
  for ( j = 0; j < Nlayers; j++ ) fprintf(fascii,"BUBLE%i\t",j);  
  for ( j = 0; j < Nlayers; j++ ) fprintf(fascii,"QURTZ%i\t",j);  
  for ( j = 0; j < Nlayers; j++ ) fprintf(fascii,"BULKDN%i\t",j);  
  for ( j = 0; j < Nlayers; j++ ) fprintf(fascii,"PARTDN%i\t",j);  
  fprintf(fascii,"OFF_GMT\t");
  for ( j = 0; j < Nlayers; j++ ) fprintf(fascii,"WcrFRC%i\t",j);  
  for ( j = 0; j < Nlayers; j++ ) fprintf(fascii,"WpFRC%i\t",j);  
  fprintf(fascii,"Z0_SOIL\tZ0_SNOW\tAVG_PREC\t");
  for ( j = 0; j < Nlayers; j++ ) fprintf(fascii,"RESMST%i\t",j);  
  fprintf(fascii,"FSACT\n");
  for ( i = 0; i < Ncells; i++ ) {
    fprintf(fascii,"1\t%i\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t",
	    soil_con[i].gridcel, cell_lat[i], cell_lng[i], 
	    soil_con[i].b_infilt, soil_con[i].Ds, soil_con[i].Dsmax,
	    soil_con[i].Ws, soil_con[i].c);
    for ( j = 0; j < Nlayers; j++ ) 
      fprintf(fascii,"%f\t", soil_con[i].expt[j]);  
    for ( j = 0; j < Nlayers; j++ ) 
      fprintf(fascii,"%f\t", soil_con[i].Ksat[j]);  
    for ( j = 0; j < Nlayers; j++ ) 
      fprintf(fascii,"%f\t", soil_con[i].phi_s[j]);  
    for ( j = 0; j < Nlayers; j++ ) 
      fprintf(fascii,"%f\t", soil_con[i].init_moist[j]);  
    fprintf(fascii,"%f\t", soil_con[i].elevation);
    for ( j = 0; j < Nlayers; j++ ) 
      fprintf(fascii,"%f\t", soil_con[i].depth[j]);  
    fprintf(fascii,"%f\t%f\t", soil_con[i].avg_temp, soil_con[i].dp);
    for ( j = 0; j < Nlayers; j++ ) 
      fprintf(fascii,"%f\t", soil_con[i].bubble[j]);  
    for ( j = 0; j < Nlayers; j++ ) 
      fprintf(fascii,"%f\t", soil_con[i].quartz[j]);  
    for ( j = 0; j < Nlayers; j++ ) 
      fprintf(fascii,"%f\t", soil_con[i].bulk_density[j]);  
    for ( j = 0; j < Nlayers; j++ ) 
      fprintf(fascii,"%f\t", soil_con[i].soil_density[j]);  
    fprintf(fascii,"%f\t", soil_con[i].off_gmt);
    for ( j = 0; j < Nlayers; j++ ) 
      fprintf(fascii,"%f\t", soil_con[i].Wcr_FRACT[j]);  
    for ( j = 0; j < Nlayers; j++ ) 
      fprintf(fascii,"%f\t", soil_con[i].Wpwp_FRACT[j]);  
    fprintf(fascii,"%f\t%f\t%f\t", soil_con[i].rough, 
	    soil_con[i].snow_rough, soil_con[i].annual_prec);
    for ( j = 0; j < Nlayers; j++ ) 
      fprintf(fascii,"%f\t", soil_con[i].resid_moist[j]);  
    if ( soil_con[i].FS_ACTIVE ) 
      fprintf(fascii,"1\n");
    else
      fprintf(fascii,"0\n");
  }

  return(0);
  
}

void read_arcinfo_mask(FILE   *fmask, 
		       int    *Ncells, 
		       float **cell_lat, 
		       float **cell_lng, 
		       int   **cell_row, 
		       int   **cell_col) {
  
  char  tmpstr[512];
  int   ncol, nrow, NoData_value;
  int   cell, i, j;
  int   tmpvalue;
  float ll_lat, ll_lng, cellsize;
  float tmp_lat, tmp_lng;
  
  /** Skip and copy header **/
  for ( i = 0; i < 6; i++ ) {
    fgets(tmpstr,512,fmask);
    if(i==0) sscanf(tmpstr,"%*s %i",&ncol);
    if(i==1) sscanf(tmpstr,"%*s %i",&nrow);
    if(i==2) sscanf(tmpstr,"%*s %f",&ll_lat);
    if(i==3) sscanf(tmpstr,"%*s %f",&ll_lng);
    if(i==4) sscanf(tmpstr,"%*s %f",&cellsize);
    if(i==5) sscanf(tmpstr,"%*s %i",&NoData_value);
  }
  
  (*Ncells) = 0;
  for ( j = 0; j < nrow; j++) {
    for( i = 0; i < ncol; i++ ) {
      fscanf(fmask,"%i",&tmpvalue);
      if ( tmpvalue > 0 && tmpvalue != NoData_value ) (*Ncells)++;
    }
  }
  
  (*cell_lat) = (float *)calloc((*Ncells),sizeof(float));
  (*cell_lng) = (float *)calloc((*Ncells),sizeof(float));
  (*cell_row) = (int   *)calloc((*Ncells),sizeof(int  ));
  (*cell_col) = (int   *)calloc((*Ncells),sizeof(int  ));

  rewind(fmask);
  for ( i = 0; i < 6; i++ ) fgets(tmpstr,512,fmask);

  cell = 0;
  for ( j = 0; j < nrow; j++ ) {
    tmp_lat = ll_lat + (float)(nrow - j - 0.5) * cellsize;
    for ( i = 0; i < ncol; i++) {
      tmp_lng = ll_lng + (float)(i + 0.5) * cellsize;
      fscanf(fmask,"%i",&tmpvalue);
      if ( tmpvalue > 0 && tmpvalue != NoData_value ) {
	(*cell_lat)[cell] = tmp_lat;
	(*cell_lng)[cell] = tmp_lng;
	(*cell_row)[cell] = j;
	(*cell_col)[cell] = i;
	cell++;
      }
    }
  }
}

void get_arc_info_values(FILE  *farc, 
			 int    Ncells, 
			 float *cell_lat, 
			 float *cell_lng, 
			 float *values) {
  
  char  tmpstr[512];
  int   ncol, nrow, NoData_value;
  int   cell, i, j;
  float ll_lat, ll_lng, cellsize;
  float tmp_lat, tmp_lng, tmpvalue;

  /** Skip and copy header **/
  for ( i = 0; i < 6; i++ ) {
    fgets(tmpstr,512,farc);
    if(i==0) sscanf(tmpstr,"%*s %i",&ncol);
    if(i==1) sscanf(tmpstr,"%*s %i",&nrow);
    if(i==2) sscanf(tmpstr,"%*s %f",&ll_lat);
    if(i==3) sscanf(tmpstr,"%*s %f",&ll_lng);
    if(i==4) sscanf(tmpstr,"%*s %f",&cellsize);
    if(i==5) sscanf(tmpstr,"%*s %i",&NoData_value);
  }

  cell = 0;
  for ( j = 0; j < nrow; j++ ) {
    tmp_lat = ll_lat + (float)(nrow - j - 0.5) * cellsize;
    for ( i = 0; i < ncol; i++) {
      tmp_lng = ll_lng + (float)(i + 0.5) * cellsize;
      fscanf(farc,"%f",&tmpvalue);
      if ( tmp_lat == cell_lat[cell] && tmp_lng == cell_lng[cell] ) {
	values[cell] = tmpvalue;
	cell++;
      }
    }
  }
}

