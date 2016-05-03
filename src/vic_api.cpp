//
// Created by sibada on 16-5-2.
//

#include <vic.h>

using std::string;

int get_soil_param(soil_con_struct ** soil_params, int ncell, string soil_path){

    extern veg_lib_struct* veg_lib;
    if(veg_lib == NULL)return ERROR;

    char read_done;
    char is_run;

    FILE * soil_in;
    soil_in = open_file(soil_path.c_str(), "r");
    if(soil_in == NULL) return ERROR;

    read_done = FALSE;

    for(int i = 0; i < ncell; i ++){
        if(read_done)break;

        soil_con_struct tmp = read_soilparam(soil_in, &is_run, &read_done);
        if(!is_run){
            soil_params[i] = NULL;
            continue;
        }

        soil_con_struct* tmpptr = new soil_con_struct(tmp);
        soil_params[i] = tmpptr;
    }

    fclose(soil_in);
    return 0;
}

void destory_soil_param(soil_con_struct** soil_params, int ncell){
    for(int i = 0; i < ncell; i++){
        if(soil_params[i] != NULL) delete soil_params[i];
    }
}

int get_veg_lib(string veglib_path){

    extern veg_lib_struct* veg_lib;
    int n_veg_type = 0;

    FILE * veglib_in;
    veglib_in = open_file(veglib_path.c_str(), "r");
    if(veglib_in == NULL) return ERROR;

    veg_lib = read_veglib(veglib_in, &n_veg_type);

    fclose(veglib_in);
    return n_veg_type;
}

int get_global(string global_path){
    extern option_struct options;
    extern global_param_struct global_param;

    FILE * global_in;
    global_in = open_file(global_path.c_str(), "r");
    if(global_in == NULL) return ERROR;

    filenames_struct filenames;

    initialize_global();
    global_param = get_global_param(&filenames, global_in);

    fclose(global_in);
    return 0;
}