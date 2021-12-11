/******************************************************************************
 * @section DESCRIPTION
 *
 * Open a file named by string and associate a stream with it.
 *****************************************************************************/

#include <vic_driver_shared_all.h>

/******************************************************************************
 * @brief    Open a file named by string and associate a stream with it.
 *
 * @param    string path to file
 * @type     Type has one of the associated values with it:
 *             - "r"    open for reading
 *             - "w"    truncate or create for writing
 *             - "a"    append; open for writing at end of file, or create for
 *                      writing
 *             - "r+"   open for update (reading and writing)
 *             - "w+"   truncate or create for update
 *             - "a+"   append; open or create for update at end-of-file
 * @return   a pointer to the file structure associated with the stream.
 *****************************************************************************/
FILE *
open_file(char string[],
          char type[])
{
    FILE *stream;
    char  zipname[MAXSTRING],
          command[MAXSTRING],
          jnkstr[MAXSTRING];
    int   temp, headcnt, i;

    stream = fopen(string, type);

    if (stream == NULL) {
        /** Check if file is compressed **/
        strcpy(zipname, string);
        strcat(zipname, ".gz");
        stream = fopen(zipname, type);
        if (stream == NULL) {
            log_err("Unable to open File %s", string);
        }
        fclose(stream);

        /** uncompress and open zipped file **/
        sprintf(command, "gzip -d %s", zipname);
        system(command);
        stream = fopen(string, type);
        if (stream == NULL) {
            log_err("Unable to open File %s", string);
        }
    }

    if (strcmp(type, "r") == 0) {
        temp = fgetc(stream);
        while (temp == 32) {
            temp = fgetc(stream);
        }
        if (temp == 35) {
            headcnt = 0;
            while (temp == 35) {
                fgets(jnkstr, MAXSTRING, stream);
                temp = fgetc(stream);
                while (temp == 32) {
                    temp = fgetc(stream);
                }
                headcnt++;
            }
            rewind(stream);
            for (i = 0; i < headcnt; i++) {
                fgets(jnkstr, MAXSTRING, stream);
            }
        }
        else {
            rewind(stream);
        }
    }

    fflush(stderr);

    return stream;
}
