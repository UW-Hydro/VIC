//
// Created by sibada on 16-5-2.
//

#include <cstdio>
#include <string>
#include <vector>
#include <vicNl.h>
#include <vic.h>

using std::vector;
using std::string;

int get_soil_param(vector<soil_con_struct*> * soil_params, string soil_path){

    char read_done;
    char is_run;

    FILE * soil_in;
    soil_in = fopen(soil_path.c_str(), "r");
    if(soil_in == NULL) return ERROR;

    read_done = FALSE;
    while(!read_done){
        soil_con_struct tmp = read_soilparam(soil_in, &read_done, &is_run);

        if(!is_run){
            soil_params->push_back(NULL);
            continue;
        }
        soil_con_struct* tmpptr = new soil_con_struct(tmp);

        soil_params->push_back(tmpptr);


    }

    fclose(soil_in);
    return 0;
}
