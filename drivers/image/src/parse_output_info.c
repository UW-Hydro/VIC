#include <vic_def.h>
#include <vic_run.h>
#include <vic_driver_image.h>

// Determine which variables will be output.
int
parse_output_info(FILE             *gp,
                  out_data_struct **out_data)
{
    char   cmdstr[MAXSTRING];
    char   optstr[MAXSTRING];
    char   varname[MAXSTRING];
    bool   found;
    int    outvarnum;
    size_t i;

    rewind(gp);
    fgets(cmdstr, MAXSTRING, gp);
    outvarnum = 0;
    while (!feof(gp)) {
        found = false;
        if (cmdstr[0] != '#' && cmdstr[0] != '\n' && cmdstr[0] != '\0') {
            sscanf(cmdstr, "%s", optstr);
            if (strcasecmp("OUTVAR", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", varname);
                for (i = 0; i < N_OUTVAR_TYPES; i++) {
                    if (strcmp(out_data[0][i].varname, varname) == 0) {
                        found = true;
                        out_data[0][i].write = true;
                    }
                }
                if (!found) {
                    fprintf(stderr, "Error: parse_output_info: \"%s\" was "
                            "not found in the list of supported output "
                            "variable names.  Please use "
                            "the exact name listed in vic_def.h.\n",
                            varname);
                }
                outvarnum++;
            }
        }
        fgets(cmdstr, MAXSTRING, gp);
    }

    return outvarnum;
}
