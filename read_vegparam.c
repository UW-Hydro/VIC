#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
#include <string.h>

static char vcid[] = "$Id$";

veg_con_struct *read_vegparam(FILE *vegparam,
                              int   gridcel,
                              int   Nveg_type)
/**********************************************************************
  read_vegparam.c    Keith Cherkauer and Dag Lohmann       1997

  This routine reads in vegetation parameters for the current grid cell.
  It also relates each vegetation class in the cell to the appropriate
  parameters in the vegetation library.

  Modifications:
  09-24-98  Modified to read root zone distribution information so
           that soil layer root fractions can be computed for new 
	   soil layer depths - see calc_root_fractions.c           KAC
  07-15-99 Modified to read LAI values from a new line in the vegetation
           parameter file.  Added specifically to work with the new
	   global LAI files.
  11-18-02 Added code to read in blowing snow parameters.          LCB
  03-27-03 Modified code to update Wdmax based on LAI values read in
           for the current grid cell.  If LAI is not obtained from this
           function, then the values cacluated in read_veglib.c are
           left unchanged.                                   DP & KAC
  2006-Nov-07 Allocates MaxVeg+1 veg tiles.  TJB
  2007-May-11 Changed some 'fscanf' statements to 'fgets' and 'sscanf' 
           to count rootzone and BLOWING fields. Also tests for fetch < 1.GCT
**********************************************************************/
{

  void ttrim( char *string );
  extern veg_lib_struct *veg_lib;
  extern option_struct   options;
#if LINK_DEBUG
  extern debug_struct    debug;
#endif

  veg_con_struct *temp;
  int             vegcel, i, j, vegetat_type_num, skip, veg_class;
  int             MaxVeg;
  int             NF, NTOT;
  float           depth_sum;
  float           sum;
  char            str[500];
  char            ErrStr[MAXSTRING];
  char            line[MAXSTRING];
  char            tmpline[MAXSTRING];
  const char      delimiters[] = " ";
  char            *token;
  char            *vegarr[500];

  if(options.GLOBAL_LAI) skip=2;
  else skip=1;

#if !NO_REWIND
  rewind(vegparam);
#endif  
    


  while ( ( fscanf(vegparam, "%d %d", &vegcel, &vegetat_type_num) == 2 ) && vegcel != gridcel ){
    for (i = 0; i <= vegetat_type_num * skip; i++){
      if ( fgets(str, 500, vegparam) == NULL ){
        sprintf(ErrStr,"ERROR unexpected EOF for cell %i while reading root zones and LAI\n",vegcel);
        nrerror(ErrStr);
      }
    }
  }
  fgets(str, 500, vegparam); // read newline at end of veg class line to advance to next line
  if (vegcel != gridcel) {
    fprintf(stderr, "Error in vegetation file.  Grid cell %d not found\n",
            gridcel);
    exit(99);
  }
  if(vegetat_type_num >= MAX_VEG) {
    sprintf(ErrStr,"Vegetation parameter file wants more vegetation types in grid cell %i (%i) than are defined by MAX_VEG (%i) [NOTE: bare soil class is assumed].  Edit vicNl_def.h and recompile.",gridcel,vegetat_type_num+1,MAX_VEG);
    nrerror(ErrStr);
  }

  // Make sure to allocate extra memory for wetlands vegetation
  if ( options.LAKES )
    MaxVeg = vegetat_type_num+1;
  else
    MaxVeg = vegetat_type_num;

  /** Allocate memory for vegetation grid cell parameters **/
  if ( MaxVeg > 0 )
    temp = (veg_con_struct*) calloc( MaxVeg+1, sizeof(veg_con_struct));
  else
    temp = (veg_con_struct*) calloc(1, sizeof(veg_con_struct));
  temp[0].Cv_sum = 0.0;

  for (i = 0; i < vegetat_type_num; i++) {
    temp[i].zone_depth = calloc(options.ROOT_ZONES,sizeof(float));
    temp[i].zone_fract = calloc(options.ROOT_ZONES,sizeof(float));
    temp[i].vegetat_type_num = vegetat_type_num;

    // Read the root zones line
    if ( fgets( line, MAXSTRING, vegparam ) == NULL ){
      sprintf(ErrStr,"ERROR unexpected EOF for cell %i while reading vegetat_type_num %d\n",vegcel,vegetat_type_num);
      nrerror(ErrStr);
    }
    strcpy(tmpline, line);
    ttrim( tmpline );
    token = strtok (tmpline, delimiters);    /*  token => veg_class, move 'line' pointer to next field */  
    NF = 0;
    vegarr[NF] = calloc( 500, sizeof(char));
    strcpy(vegarr[NF],token);
    NF++;

    while( ( token = strtok (NULL, delimiters)) != NULL ){
      vegarr[NF] = calloc( 500, sizeof(char));      
      strcpy(vegarr[NF],token);
      NF++;
    }

    NTOT = 2 + 2 * options.ROOT_ZONES;  /* Number of expected fields this line */
    if( options.BLOWING ){
      NTOT += 3;
    }
    if ( NF != NTOT ) {
      sprintf(ErrStr,"ERROR - cell %d - expecting %d fields but found %d in veg line %s\n",gridcel,NTOT, NF, line);
      nrerror(ErrStr);
    }

    temp[i].veg_class = atoi( vegarr[0] );
    temp[i].Cv = atof( vegarr[1] );
    depth_sum = 0;
    sum = 0.;
    for(j=0;j<options.ROOT_ZONES;j++) {
      temp[i].zone_depth[j] = atof( vegarr[2 + j*2] );
      temp[i].zone_fract[j] = atof( vegarr[3 + j*2] );
      depth_sum += temp[i].zone_depth[j];
      sum += temp[i].zone_fract[j];
    }
    if(depth_sum <= 0) {
      sprintf(str,"Root zone depths must sum to a value greater than 0.");
      nrerror(str);
    }
    if(sum != 1.) {
      fprintf(stderr,"WARNING: Root zone fractions sum to more than 1 ( = %f), normalizing fractions.  If the sum is large, check that your vegetation parameter file is in the form - <zone 1 depth> <zone 1 fract> <zone 2 depth> <zone 2 fract> ...\n", sum);
      for(j=0;j<options.ROOT_ZONES;j++) {
	temp[i].zone_fract[j] /= sum;
      }
    }

    if(options.BLOWING) {
      j = 2 * options.ROOT_ZONES;
      temp[i].sigma_slope = atof( vegarr[2 + j] );
      temp[i].lag_one = atof( vegarr[3 + j] );
      temp[i].fetch = atof( vegarr[4 + j]) ;
      if( temp[i].sigma_slope <= 0. || temp[i].lag_one <= 0.) {
        sprintf(str,"Deviation of terrain slope must be greater than 0.");
        nrerror(str);
      }
      if( temp[i].fetch < 1.0  ) {
	sprintf(str,"ERROR - BLOWING parameter fetch should be >> 1 but cell %i has fetch = %.2f\n", gridcel, temp[i].fetch );
        nrerror(str);
      }
    }

    veg_class = MISSING;
    for(j=0;j<Nveg_type;j++)
      if(temp[i].veg_class == veg_lib[j].veg_class)
	veg_class = j;
    if(veg_class == MISSING) {
      sprintf(ErrStr,"Vegetation class %i from cell %i is not defined in the vegetation library file.", temp[i].veg_class, gridcel);
      nrerror(ErrStr);
    }
    else
      temp[i].veg_class = veg_class;

    temp[0].Cv_sum += temp[i].Cv;

    if ( options.GLOBAL_LAI )
      // Read the LAI line
      if ( fgets( line, MAXSTRING, vegparam ) == NULL ){
        sprintf(ErrStr,"ERROR unexpected EOF for cell %i while reading LAI for vegetat_type_num %d\n",vegcel,vegetat_type_num);
        nrerror(ErrStr);
      }      
      NF = 0;
      strcpy(tmpline, line);
      ttrim( tmpline );
      token = strtok (tmpline, delimiters); 
      strcpy(vegarr[NF],token);
      NF++;
 
      while( ( token = strtok (NULL, delimiters)) != NULL ){
        vegarr[NF] = calloc( 500, sizeof(char));      
        strcpy(vegarr[NF],token);
        NF++;
     }
     NTOT = 12; /* For LAI */
     if ( NF != NTOT ) {
       sprintf(ErrStr,"ERROR - cell %d - expecting %d LAI values but found %d in line %s\n",gridcel, NTOT, NF, line);
       nrerror(ErrStr);
     }

     for ( j = 0; j < 12; j++ ) {
       veg_lib[temp[i].veg_class].LAI[j] = atof( vegarr[j] );
       veg_lib[temp[i].veg_class].Wdmax[j] = 
	  LAI_WATER_FACTOR * veg_lib[temp[i].veg_class].LAI[j];
      }
  }
  if(temp[0].Cv_sum>1.0){
    fprintf(stderr,"WARNING: Cv exceeds 1.0 at grid cell %d, fractions being adjusted to equal 1\n", gridcel);
    for(j=0;j<vegetat_type_num;j++)
      temp[j].Cv = temp[j].Cv / temp[0].Cv_sum;
    temp[0].Cv_sum = 1.;
  }
  if(temp[0].Cv_sum>0.99 && temp[0].Cv_sum<1.0){
    fprintf(stderr,"WARNING: Cv > 0.99 and Cv < 1.0 at grid cell %d, model assuming that bare soil is not to be run - fractions being adjusted to equal 1\n", gridcel);
    for(j=0;j<vegetat_type_num;j++)
      temp[j].Cv = temp[j].Cv / temp[0].Cv_sum;
    temp[0].Cv_sum = 1.;
  }
  return temp;
} 






/* trim trailing newlines */

#define END '\0'
#define NEW '\n'

void ttrim( char *c ) 
{
  while( (*c++ != END) );
    --c;
  for( ; *--c == NEW; *c = END );

}

