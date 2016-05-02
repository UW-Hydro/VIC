//
// Created by sibada on 16-5-2.
//

#ifndef VICNL_VIC_H
#define VICNL_VIC_H

#include <vicNl.h>
#include <string>
#include <vector>

using std::vector;
using std::string;

int get_soil_param(vector<soil_con_struct*> * soil_params, string soil_path);

#endif //VICNL_VIC_H
