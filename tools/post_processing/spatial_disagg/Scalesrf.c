/*
 * Purpose: aggregate 5 minute .srf files produced by the soil data program to
 *          other resolutions (multiple of 5 minutes)
 * Usage:   Scalesrf <input file root> <output file root> <aggregate resolution>
 * Author:  Bart Nijssen
 * Created: Thu Jun 17 14:37:31 1999
 * Last Changed: Thu Jul 15 16:04:07 1999 by Bart Nijssen <nijssen@u.washington.edu>
 *               <nijssen@u.washington.edu> 
 * Comment: Input is <input file>.srf and <input file>.doc.  Output is sent to
 *          <output file>.srf and <output file>.doc
 *          Talking about bloat.....
 */

/******************************************************************************/
/*				    includes                                  */
/******************************************************************************/
#include <ctype.h>
#include <errno.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/******************************************************************************/
/*			defines, typedefs, globals, etc.                      */
/******************************************************************************/
#define ENOERROR 0		/* no error */
#define EUSAGE 1000		/* usage error */
#define ENOELEMENT 1001		/* element not in list */
#define ENODOUBLE 1002		/* not a double */
#define ENOFLOAT 1003		/* not a float */
#define ENOINT 1004		/* not an integer */
#define ENOLONG 1005		/* not a long */
#define ENOSHORT 1006		/* not a short */
#define ENOCASE 1007		/* not a valid case in a switch statement */
#define ENOTINFILE 1008		/* not in the file */
#define ENOFLOW 1009		/* no flow data */
#define ENORES 1010		/* not a valid resolution */


#ifndef FALSE
#define FALSE 0
#endif
#ifndef TRUE
#define TRUE 1
#endif
#define MISSING -9999

#define OCEAN -2
#define VOID -1
#define EPS 1e-7

const int ARRAY_INCREMENT = 10; /* number of slots by which to increase the
				   neighbors array when it gets full */ 

typedef enum {double_f, int_f, float_f, long_f, short_f} FORMAT_SPECIFIER;

extern int errno;
char message[BUFSIZ+1] = "";
const char *docext = ".doc";
const char *srfext = ".srf";
const char *usage = 
"\nUsage: %s <5 min .srf file> <output file> <output resolution>\n\n";
int status = ENOERROR;

/******************************************************************************/
/*				   prototypes                                 */
/******************************************************************************/
int GetNumber(char *str, int format, int start, int end, void *value);
int ProcessCommandLine(int argc, char **argv, char *infilename, 
		       char *outfilename, int *resolution);
int ProcessError(void);
int ReadDocFile(char *infilename, char ***doc, int *inrows, int *incols);
int ReadWriteSrfFiles(char *infilename, char *outfilename, int resolution, 
		      int inrows, int incols, float *minval, float *maxval);
int SetToMissing(int format, void *value);
int WriteDocFile(char *outfilename, char **doc, int inrows, int incols, 
		 int resolution, float minval, float maxval);

/******************************************************************************/
/******************************************************************************/
/*				      main                                    */
/******************************************************************************/
/******************************************************************************/
int main(int argc, char **argv)
{
  char infilename[BUFSIZ+1];
  char outfilename[BUFSIZ+1];
  char **doc = NULL;
  float maxval;
  float minval;
  int i;
  int resolution;
  int incols;
  int inrows;
  
  status = ProcessCommandLine(argc, argv, infilename, outfilename, &resolution);
  if (status != ENOERROR)
    goto error;

  status = ReadDocFile(infilename, &doc, &inrows, &incols);
  if (status != ENOERROR)
    goto error;

  status = ReadWriteSrfFiles(infilename, outfilename, resolution, inrows, 
 			     incols, &minval, &maxval);  
  if (status != ENOERROR) 
    goto error; 

  status = WriteDocFile(outfilename, doc, inrows, incols, resolution, minval,
 			maxval);
  if (status != ENOERROR)
    goto error;

  if (doc != NULL) {
    i = 0;
    while (doc[i]) {
      free(doc[i]);
      i++;
    }
    free(doc);
  }

  return EXIT_SUCCESS;

 error:
  if (doc != NULL) {
    i = 0;
    while (doc[i]) {
      if (doc[i] != NULL)
	free(doc[i]);
      else
	break;
      i++;
    }
    free(doc);
  }
  ProcessError();
  exit(EXIT_FAILURE);
}

/******************************************************************************/
/*			       ReadWriteDocFiles                              */
/******************************************************************************/
int ReadDocFile(char *infilename, char ***doc, int *inrows, int *incols)
{
  char filename[BUFSIZ+1];
  char str[BUFSIZ+1];
  FILE *infile = NULL;
  int arraysize;
  int ndatalines;

  sprintf(filename, "%s%s", infilename, docext);
  infile = fopen(filename, "r");
  if (infile == NULL) {
    status = errno;
    strcpy(message, filename);
    goto error;
  }

  arraysize = ARRAY_INCREMENT;
  *doc = calloc(arraysize, sizeof(char *));
  if (*doc == NULL) {
    status = errno;
    strcpy(message, "doc array");
    goto error;
  } 
  ndatalines = 0;
  while (fgets(str, BUFSIZ, infile) != NULL) {
    if (arraysize < ndatalines + 1) {
      arraysize += ARRAY_INCREMENT;
      *doc = realloc(*doc, arraysize * sizeof(char *));
      if (*doc == NULL) {
	status = errno;
	strcpy(message, "doc array");
	goto error;
      }
    }
    (*doc)[ndatalines] = calloc(strlen(str)+1, sizeof(char));
    if ((*doc)[ndatalines] == NULL) {
      status = errno;
      strcpy(message, "doc array string");
      goto error;
    }
    strcpy((*doc)[ndatalines], str);
    ndatalines++;
  } 
  
  sscanf((*doc)[3], "%*s %*s %d", incols);
  sscanf((*doc)[4], "%*s %*s %d", inrows);

  if (ferror(infile)) {
    strcpy(message, infilename);
    status = errno;
    goto error;
  }

  fclose(infile);

  return ENOERROR;

 error:
  if (infile != NULL)
    fclose(infile);
  return status;
}

/******************************************************************************/
/*			       ReadWriteSrfFiles                              */
/******************************************************************************/
int ReadWriteSrfFiles(char *infilename, char *outfilename, int resolution, 
		      int inrows, int incols, float *minval, float *maxval)
{
  FILE *infile = NULL;
  FILE *outfile = NULL;
  char filename[BUFSIZ+1];
  float datavalue;
  float **val = NULL;
  int **count = NULL;
  int **nocean = NULL;
  int i;
  int j;
  int outcols;
  int outrows;
  
  outcols = incols/resolution;
  outrows = inrows/resolution;

  val = calloc(outrows, sizeof(float *));
  if (val == NULL) {
    status = errno;
    strcpy(message, "val array");
    goto error;
  }
  for (i = 0; i < outrows; i++) {
    val[i] = calloc(outcols, sizeof(float));
    if (val[i] == NULL) {
      status = errno;
      strcpy(message, "val array");
      goto error;
    }
  }
  count = calloc(outrows, sizeof(int *));
  if (count == NULL) {
    status = errno;
    strcpy(message, "count array");
    goto error;
  }
  for (i = 0; i < outrows; i++) {
    count[i] = calloc(outcols, sizeof(int));
    if (count[i] == NULL) {
      status = errno;
      strcpy(message, "count array");
      goto error;
    }
  }
  nocean = calloc(outrows, sizeof(int *));
  if (nocean == NULL) {
    status = errno;
    strcpy(message, "nocean array");
    goto error;
  }
  for (i = 0; i < outrows; i++) {
    nocean[i] = calloc(outcols, sizeof(int));
    if (nocean[i] == NULL) {
      status = errno;
      strcpy(message, "nocean array");
      goto error;
    }
  }
  
  sprintf(filename, "%s%s", infilename, srfext);
  infile = fopen(filename, "r");
  if (infile == NULL) {
    status = errno;
    strcpy(message, filename);
    goto error;
  }

  sprintf(filename, "%s%s", outfilename, srfext);
  outfile = fopen(filename, "w");
  if (outfile == NULL) {
    status = errno;
    strcpy(message, filename);
    goto error;
  }

  for (i = 0; i < outrows*resolution; i++) {
    for (j = 0; j < outcols*resolution; j++) {
      fscanf(infile, "%f", &datavalue);
      if (fabs(datavalue-VOID) > EPS) {
	if (fabs(datavalue-OCEAN) < EPS) 
	  nocean[i/resolution][j/resolution] += 1;
	else {
	  val[i/resolution][j/resolution] += datavalue;
	  count[i/resolution][j/resolution] += 1;
	}
      }
    }
    for ( ; j < incols; j++)
      fscanf(infile, "%f", &datavalue);
  }
  
  fclose(infile);

  *maxval = FLT_MIN;
  *minval = FLT_MAX;
  for (i = 0; i < outrows; i++) {
    for (j = 0; j < outcols; j++) {
      if (nocean[i][j] == resolution*resolution)
	val[i][j] = OCEAN;
      else if (count[i][j] == 0)
	val[i][j] = VOID;
      else {
	val[i][j] /= count[i][j];
	if (val[i][j] > *maxval)
	  *maxval = val[i][j];
	if (val[i][j] < *minval)
	  *minval = val[i][j];
      }
      fprintf(outfile, "%f\n", val[i][j]);
    }
  }
  if (*minval > *maxval) {
    *minval = VOID;
    *maxval = VOID;
  }

  fclose(outfile);

  for (i = 0; i < outrows; i++) {
    free(val[i]);
    free(count[i]);
    free(nocean[i]);
  }
  free(val);
  free(count);
  free(nocean);

  return ENOERROR;

 error:
  if (infile != NULL)
    fclose(infile);
  if (outfile != NULL)
    fclose(outfile);
  if (val != NULL) {
    for (i = 0; i < outrows; outrows++) {
      if (val[i] != NULL)
	free (val[i]);
      else
	break;
    }
    free(val);
  }
  if (count != NULL) {
    for (i = 0; i < outrows; outrows++) {
      if (count[i] != NULL)
	free(count[i]);
      else
	break;
    }
    free(count);
  }
  if (nocean != NULL) {
    for (i = 0; i < outrows; outrows++) {
      if (nocean[i] != NULL)
	free(nocean[i]);
      else
	break;
    }
    free(nocean);
  }
  return status;
}

/******************************************************************************/
/*				  WriteDocFile                                */
/******************************************************************************/
int WriteDocFile(char *outfilename, char **doc, int inrows, int incols, 
		 int resolution, float minval, float maxval)
{
  char filename[BUFSIZ+1];
  char str1[BUFSIZ+1];
  char str2[BUFSIZ+1];
  FILE *outfile = NULL;
  float minx;
  float maxy;
  int i;

  for (i = 0; doc[3][i] != '\0' && doc[3][i] != ':'; i++)
    str1[i] = doc[3][i];
  str1[i] = '\0';
  sprintf(str2, "%s: %d\n", str1, incols/resolution);
  doc[3] = realloc(doc[3], (strlen(str2)+1) * sizeof(char));
  if (doc[3] == NULL) {
    strcpy(message, "realloc doc[3]");
    status = errno;
    goto error;
  }
  strcpy(doc[3], str2);

  for (i = 0; doc[4][i] != '\0' && doc[4][i] != ':'; i++)
    str1[i] = doc[4][i];
  str1[i] = '\0';
  sprintf(str2, "%s: %d\n", str1, inrows/resolution);
  doc[4] = realloc(doc[4], (strlen(str2)+1) * sizeof(char));
  if (doc[4] == NULL) {
    strcpy(message, "realloc doc[4]");
    status = errno;
    goto error;
  }
  strcpy(doc[4], str2);
  
  for (i = 0; doc[7][i] != '\0' && doc[7][i] != ':'; i++)
    str1[i] = doc[7][i];
  str1[i] = '\0';
  sprintf(str2, "%s: %f * %f (lat * long)\n", str1, resolution/12.,
	  resolution/12.);
  doc[7] = realloc(doc[7], (strlen(str2)+1) * sizeof(char));
  if (doc[7] == NULL) {
    strcpy(message, "realloc doc[7]");
    status = errno;
    goto error;
  }
  strcpy(doc[7], str2);

  for (i = 0; doc[9][i] != '\0' && doc[9][i] != ':'; i++)
    str1[i] = doc[9][i];
  str1[i] = '\0';
  sscanf(&(doc[8][i+1]), "%f", &minx);
  sprintf(str2, "%s: %f\n", str1, minx+(incols/12));
  doc[9] = realloc(doc[9], (strlen(str2)+1) * sizeof(char));
  if (doc[9] == NULL) {
    strcpy(message, "realloc doc[9]");
    status = errno;
    goto error;
  }
  strcpy(doc[9], str2);
  
  for (i = 0; doc[10][i] != '\0' && doc[10][i] != ':'; i++)
    str1[i] = doc[10][i];
  str1[i] = '\0';
  sscanf(&(doc[11][i+1]), "%f", &maxy);
  sprintf(str2, "%s: %f\n", str1, maxy-(inrows/12));
  doc[10] = realloc(doc[10], (strlen(str2)+1) * sizeof(char));
  if (doc[10] == NULL) {
    strcpy(message, "realloc doc[11]");
    status = errno;
    goto error;
  }
  strcpy(doc[10], str2);
  
  for (i = 0; doc[14][i] != '\0' && doc[14][i] != ':'; i++)
    str1[i] = doc[14][i];
  str1[i] = '\0';
  sprintf(str2, "%s: %f\n", str1, minval);
  doc[14] = realloc(doc[14], (strlen(str2)+1) * sizeof(char));
  if (doc[14] == NULL) {
    strcpy(message, "realloc doc[14]");
    status = errno;
    goto error;
  }
  strcpy(doc[14], str2);
  
  for (i = 0; doc[15][i] != '\0' && doc[15][i] != ':'; i++)
    str1[i] = doc[15][i];
  str1[i] = '\0';
  sprintf(str2, "%s: %f\n", str1, maxval);
  doc[15] = realloc(doc[15], (strlen(str2)+1) * sizeof(char));
  if (doc[15] == NULL) {
    strcpy(message, "realloc doc[15]");
    status = errno;
    goto error;
  }
  strcpy(doc[15], str2);

  sprintf(filename, "%s%s", outfilename, docext);
  outfile = fopen(filename, "w");
  if (outfile == NULL) {
    status = errno;
    strcpy(message, filename);
    goto error;
  }

  for (i = 0; doc[i] != NULL; i++) {
    status = fputs(doc[i], outfile);
    if (status < 0) {
      strcpy(message, filename);
      status = errno;
      goto error;
    }
  }

  fclose(outfile);
  

  return ENOERROR;

 error:
  if (outfile != NULL)
    fclose(outfile);

  return status;
}

/******************************************************************************/
/*			       Processcommandline                             */
/******************************************************************************/
int ProcessCommandLine(int argc, char **argv, char *infilename, 
		       char *outfilename, int *resolution) 
{
  int errflg = 0;
  float res;

  if (argc  < 4) {
    fprintf(stderr, "\nMissing command-line arguments\n");
    errflg++;
  }
  else if (argc > 4 && errflg == 0 ) {
    fprintf(stderr, "\n\nToo many command-line arguments\n");
    errflg++;
  }
  
  if (errflg == 0) {
    strcpy(infilename, argv[1]);
    strcpy(outfilename, argv[2]);
    status = GetNumber(argv[3], float_f, 0, strlen(argv[3]), &res);
    if (status != ENOERROR)
      goto error;
  }
  
  *resolution = 12.* res;
  if (fabs(*resolution - 12.*res) > EPS) {
    strcpy(message, argv[3]);
    status = ENORES;
    goto error;
  }
    
  if (errflg) {
    strcpy(message, argv[0]); 
    status = EUSAGE;
    goto error;
  }
      
  return ENOERROR;

 error:
  return status;
}

/******************************************************************************/
/*				  ProcessError                                */
/******************************************************************************/
int ProcessError(void)
{
  if (errno) 
    perror(message);
  else {
    switch (status) {
    case EUSAGE:
      fprintf(stderr, usage, message);
      break;
    case ENOELEMENT:
      fprintf(stderr, "\nElement not found in list:\n%s\n\n", message);
      break;      
    case ENODOUBLE:
      fprintf(stderr, "\nNot a valid double: %s\n\n", message);
      break;
    case ENOFLOAT:
      fprintf(stderr, "\nNot a valid float:\n%s\n\n", message);
      break;      
    case ENOINT:
      fprintf(stderr, "\nNot a valid integer:\n%s\n\n", message);
      break;      
    case ENOLONG:
      fprintf(stderr, "\nNot a valid long: %s\n\n", message);
      break;
    case ENOSHORT:
      fprintf(stderr, "\nNot a valid short: %s\n\n", message);
      break;
    case ENOCASE:
      fprintf(stderr, "\nInvalid case:\n%s\n\n", message);
      break;      
    case ENOTINFILE:
      fprintf(stderr, "\nNot in file: %s\n\n", message);
      break;
    case ENORES:
      fprintf(stderr, "\nResolution not a multiple of 5 minutes: %s\n\n",
	      message);
      break;
    default:
      fprintf(stderr, "\nError: %s\n\n", message);
      break;
    }
  }
  status = ENOERROR;
  return status;
}

/******************************************************************************/
/*				   GetNumber                                  */
/******************************************************************************/
int GetNumber(char *str, int format, int start, int end, void *value)
{
  char valstr[BUFSIZ+1];
  int i;
  int nchar;
  char *endptr = NULL;
  
  nchar = end-start+1;
  strncpy(valstr, &(str[start]), nchar);
  valstr[nchar] = '\0';
  
  for (i = 0; i < nchar; i++) {
    if (!isspace(valstr[i]))
      break;
  }
  if (i == nchar) {
    status = SetToMissing(format, value);
    if (status != ENOERROR)
      goto error;
    return ENOERROR;
  }
  
  switch (format) {
  case double_f:
    *((double *)value) = strtod(valstr, &endptr);
    break;
  case float_f:
    *((float *)value) = (float) strtod(valstr, &endptr);
    break;
  case int_f:
    *((int *)value) = (int) strtol(valstr, &endptr, 0);
    break;
  case long_f:
    *((long *)value) = (int) strtol(valstr, &endptr, 0);
    break;
  case short_f:
    *((short *)value) = (short) strtol(valstr, &endptr, 0);
    break;
  default:
    strcpy(message, "Unknown number format");
    status = ENOCASE;
    goto error;
    break;
  }

  if (!(isspace(*endptr) || *endptr == '\0' || endptr == NULL) ) {
    /* check if the string consists entirely of whitespace.  If so, this is a
       missing value */
    strcpy(message, valstr);
    switch (format) {
    case double_f:
      status = ENODOUBLE;
      break;
    case float_f:
      status = ENOFLOAT;
      break;
    case int_f:
      status = ENOINT;
      break;
    case long_f:
      status = ENOLONG;
      break;
    case short_f:
      status = ENOSHORT;
      break;
    default:
      strcpy(message, "Unknown number format");
      status = ENOCASE;
      break;
    }
    goto error;
  }
  
  return ENOERROR;
  
 error:
  return status;
}

/******************************************************************************/
/*				  SetToMissing                                */
/******************************************************************************/
int SetToMissing(int format, void *value)
{
  switch (format) {
  case double_f:
    *((double *)value) = MISSING;
    break;
  case int_f:
    *((int *)value) = MISSING;
    break;
  case float_f:
    *((float *)value) = MISSING;
    break;
  case long_f:
    *((long *)value) = MISSING;
    break;
  case short_f:
    *((short *)value) = MISSING;
    break;
  default:
    strcpy(message, "Unknown number format");
    status = ENOCASE;
    goto error;
    break;
  }    
  return ENOERROR;
  
 error:
  return status;
}
