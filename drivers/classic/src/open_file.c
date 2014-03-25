#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <vicNl.h>

static char vcid[] = "$Id$";

FILE *open_file(char string[],char type[])

/******************************************************************/
/*  Open a file named by string and associate a stream with it.   */
/*  Type has one of the associated values with it:                */
/*                                                                */
/*  "r"    open for reading                                       */
/*                                                                */
/*  "w"    truncate or create for writing                         */
/*                                                                */
/*  "a"    append; open for writing at end of file, or create     */
/*         for writing                                            */
/*                                                                */
/*  "r+"   open for update (reading and writing)                  */
/*                                                                */
/*  "w+"   truncate or create for update                          */
/*                                                                */
/*  "a+"   append; open or create for update at end-of-file       */
/*                                                                */
/*  Open_file returns a pointer to the file structure associated  */
/*  with the stream.                                              */
/******************************************************************/
/* 30-Oct-03 Added message announcing the opening of files when type
	     is "rb".						TJB
 ******************************************************************/

{

  FILE *stream;
  char zipname[MAXSTRING],
       command[MAXSTRING],
       jnkstr[MAXSTRING];
  int  temp, headcnt, i;

  stream = fopen(string,type);

#if VERBOSE
  fprintf(stderr,"\n");
#endif

  if (stream == NULL) {

    /** Check if file is compressed **/
    strcpy(zipname,string);
    strcat(zipname,".gz");
    stream = fopen(zipname,type);
    if (stream == NULL) {
      fprintf(stderr,"\n Error opening \"%s\".",string);
      fprintf(stderr,"\n");
      nrerror("Unable to open File");
    }
    fclose(stream);

    /** uncompress and open zipped file **/
#if VERBOSE
    fprintf(stderr,"unzipping \"%s\".",string);
#endif

    sprintf(command,"gzip -d %s",zipname);
    system(command);
    stream = fopen(string,type);
    if (stream == NULL) {
      fprintf(stderr,"\n Error opening \"%s\".",string);
      fprintf(stderr,"\n");
      nrerror("Unable to Open File");
    }
  }

#if VERBOSE
  fprintf(stderr,"\n \"%s\" has been",string);
#endif

  if(strcmp(type,"r") == 0) {

#if VERBOSE
    fprintf(stderr,"\n  opened for reading.");
#endif

    temp=fgetc(stream);
    while(temp==32) temp=fgetc(stream);
    if(temp==35) {

#if VERBOSE
      fprintf(stderr,"... skipping header");
#endif

      headcnt = 0;
      while(temp==35) {
	fgets(jnkstr,MAXSTRING,stream);
	temp=fgetc(stream);
	while(temp==32) temp=fgetc(stream);
        headcnt++;
      }
      rewind(stream);
      for(i=0;i<headcnt;i++) fgets(jnkstr,MAXSTRING,stream);
    }
    else rewind(stream);
  }

#if VERBOSE
  if(strcmp(type,"rb") == 0) fprintf(stderr,"\n  opened for reading.");

  if(strcmp(type,"w") == 0) fprintf(stderr,"\n  truncated or created for writing.");

  if(strcmp(type,"wb") == 0) fprintf(stderr,"\n  truncated or created for writing.");

  if(strcmp(type,"a") == 0) {
    fprintf(stderr,"\n  opened or created for writing at");
    fprintf(stderr,"\n  the end-of-file.");
  }

  if(strcmp(type,"r+") == 0)
  fprintf(stderr,"\n  opened for update (reading and writing).");

  if(strcmp(type,"w+") == 0)
  fprintf(stderr,"\n  truncated or updated for update.");

  if(strcmp(type,"a+") == 0) {
    fprintf(stderr,"\n  opened or created for updating at");
    fprintf(stderr,"\n  the end-of-file.");
  }

  fprintf(stderr,"\n");
#endif

  fflush(stderr);

  return stream;
}
