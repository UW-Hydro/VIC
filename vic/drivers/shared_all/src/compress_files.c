/******************************************************************************
 * @section DESCRIPTION
 *
 * This subroutine compresses the file "string" using a system call.
 *****************************************************************************/

#include <vic_driver_shared_all.h>

/******************************************************************************
 * @brief    This subroutine compresses the file "string" using a system call.
 *****************************************************************************/
void
compress_files(char      string[],
               short int level)
{
    char command[MAXSTRING];

    // Compress the file
    if (level == COMPRESSION_LVL_DEFAULT) {
        sprintf(command, "nice gzip -f %s &", string);
    }
    else if (level != COMPRESSION_LVL_UNSET) {
        sprintf(command, "nice gzip -%d -f %s &", level, string);
    }
    else if (level <= 0) {
        log_err("Invalid compression level for gzip, must be an integer 1-9");
    }

    system(command);
}
