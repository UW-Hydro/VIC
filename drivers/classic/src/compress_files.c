#include <vic_def.h>
#include <vic_run.h>
#include <vic_driver_classic.h>

void
compress_files(char string[])

/**********************************************************************
   compress_files.c	Keith Cherkauer		September 10, 1997

   This subroutine compresses the file "string" using a system call.

**********************************************************************/
{
    char command[MAXSTRING];

    /** uncompress and open zipped file **/
    sprintf(command, "nice gzip -f %s &", string);
    system(command);
}
