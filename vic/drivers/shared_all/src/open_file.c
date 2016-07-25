/******************************************************************************
 * @section DESCRIPTION
 *
 * Open a file named by string and associate a stream with it.
 *
 * @section LICENSE
 *
 * The Variable Infiltration Capacity (VIC) macroscale hydrological model
 * Copyright (C) 2016 The Computational Hydrology Group, Department of Civil
 * and Environmental Engineering, University of Washington.
 *
 * The VIC model is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
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
