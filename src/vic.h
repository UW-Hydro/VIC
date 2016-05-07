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

#define INFILT 0
#define DS     1
#define DS_MAX 2
#define WS     3
#define D1     4
#define D2     5
#define D3     6

using std::string;

namespace vic {
    /* **************************************************************
     *
     * VIC global state parameters and
     * some structures frequently used.
     * VIC状态参数以及一些常用结构。
     *
     * **************************************************************
     */
    static bool _n_global_inited;
    static bool _soil_gotten;
    static bool _veglib_gotten;
    static bool _vegparam_gotten;
    static bool _atmos_inited;

    static filenames_struct _filenames;
    static filep_struct _filep;
    static atmos_data_struct *atmos;
    static out_data_file_struct *_out_data_files;
    static out_data_struct *_out_data;
    static save_data_struct _save_data;

    static dmy_struct* _dmy;
}

int get_veg_param(veg_con_struct **veg_params, int Nveg_type, soil_con_struct **soil_params, int ncell,
                  string veg_path);

int get_global(string global_path);

int get_veg_lib(string veglib_path);

int get_soil_params(soil_con_struct **soil_params, int ncell, string soil_path);

void destroy_soil_param(soil_con_struct **soil_params, int ncell);

void destroy_veg_param(veg_con_struct **veg_params, int ncell);

int run_a_cell(soil_con_struct *soil_con,
               veg_con_struct *veg_con,
               lake_con_struct *lake_con,
               double* out_runoffs);

void modify_soil_params(soil_con_struct** soil_params,
                        int* cellids,
                        int cellnum,
                        int param_type,
                        double value);

#endif //VIC_H
