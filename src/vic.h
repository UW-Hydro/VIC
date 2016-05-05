//
// Created by sibada on 16-5-2.
//
/**********************************************************************

  Define VIC APIs for C++.

  Created by Sibada on 2016-05-02.


*********************************************************************/
#ifndef VIC_H
#define VIC_H

#include <string>

extern "C"{
#include <cstdio>
#include <string>
#include <vicNl.h>
#include <vicNl_def.h>
#include <global.h>
};

using std::string;

namespace vic {
    bool _n_global_inited;
    bool _soil_gotten;
    bool _veglib_gotten;
    bool _vegparam_gotten;
    bool _atmos_inited;

    static filenames_struct _filenames;
    static filep_struct _filep;
    static atmos_data_struct *atmos;
    static out_data_file_struct *_out_data_files;
    static out_data_struct *_out_data;
    static save_data_struct _save_data;

    dmy_struct* _dmy;
}

int get_veg_param(veg_con_struct **veg_params, int Nveg_type, soil_con_struct **soil_params, int ncell,
                  string veg_path);

int get_global(string global_path);

int get_veg_lib(string veglib_path);

int get_soil_params(soil_con_struct **soil_params, int ncell, string soil_path);

void destory_soil_param(soil_con_struct **soil_params, int ncell);

void destory_veg_param(veg_con_struct **veg_params, int ncell);

int run_a_cell(soil_con_struct *soil_con,
               veg_con_struct *veg_con,
               lake_con_struct *lake_con);

#endif //VIC_H
