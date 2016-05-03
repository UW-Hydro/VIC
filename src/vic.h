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

int get_global(string global_path);

int get_veg_lib(string veglib_path);

int get_soil_param(soil_con_struct ** soil_params, int ncell, string soil_path);
void destory_soil_param(soil_con_struct** soil_params, int ncell);

#endif //VIC_H
