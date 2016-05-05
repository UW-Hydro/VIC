//
// Created by sibada on 16-5-2.
//
#include <iostream>
#include <cstring>

#include <vic.h>
#include "vicNl_def.h"

using namespace std;
using namespace vic;

int get_veg_param(veg_con_struct ** veg_params, int Nveg_type, soil_con_struct ** soil_params, int ncell, string veg_path) {

    FILE * veg_in;
    veg_in = open_file(veg_path.c_str(), "r");
    if(veg_in == NULL) return ERROR;

    for(int i = 0; i < ncell; i++) {
        if(soil_params[i] == NULL){
            veg_params[i] = NULL;
            continue;
        }
        veg_params[i] = read_vegparam(veg_in, soil_params[i]->gridcel,
                                Nveg_type);
        calc_root_fractions(veg_params[i], soil_params[i]);
    }

    fclose(veg_in);
    _vegparam_gotten = TRUE;
    return 0;
}

int get_soil_params(soil_con_struct **soil_params, int ncell, string soil_path){

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
        soil_params[i] = new soil_con_struct(tmp);
    }

    fclose(soil_in);
    _soil_gotten = TRUE;
    return 0;
}

void destory_soil_param(soil_con_struct** soil_params, int ncell){
    for(int i = 0; i < ncell; i++){
        if(soil_params[i] == NULL)continue;
        free((char *)soil_params[i]->AreaFract);
        free((char *)soil_params[i]->BandElev);
        free((char *)soil_params[i]->Tfactor);
        free((char *)soil_params[i]->Pfactor);
        free((char *)soil_params[i]->AboveTreeLine);
        delete soil_params[i];
    }
    delete soil_params;
    _soil_gotten = FALSE;
}

void destory_veg_param(veg_con_struct** veg_params, int ncell){
    extern option_struct options;

    for(int c = 0; c < ncell; c++){
        if(veg_params[c] == NULL) continue;
        for(int i = 0;i<veg_params[c][0].vegetat_type_num;i++) {
            free((char *)veg_params[c][i].zone_depth);
            free((char *)veg_params[c][i].zone_fract);
            if (options.CARBON) {
                free((char *)veg_params[c][i].CanopLayerBnd);
            }
        }
        free((char *)veg_params[c]);
    }
    delete veg_params;
    _vegparam_gotten = FALSE;
}

/* **************************************************************
 *
 * Get veg library data.
 * 获取植被库数据。
 *
 * **************************************************************
 */
int get_veg_lib(string veglib_path){
    if(_n_global_inited){
        cout<<"Error: Global did not initiated.\n";
        return ERROR;
    }

    extern veg_lib_struct* veg_lib;
    int n_veg_type = 0;

    FILE * veglib_in;
    veglib_in = open_file(veglib_path.c_str(), "r");
    if(veglib_in == NULL) return ERROR;

    veg_lib = read_veglib(veglib_in, &n_veg_type);

    fclose(veglib_in);
    return n_veg_type;
}

/* **************************************************************
 *
 * Get global parameters and initiate model.
 * 获取全局参数并初步初始化。
 *
 * **************************************************************
 */
int get_global(string global_path){

    extern global_param_struct global_param;

    FILE * global_in;
    global_in = open_file(global_path.c_str(), "r");
    if(global_in == NULL) return ERROR;

    initialize_global();
    global_param = get_global_param(&_filenames, global_in);
    fclose(global_in);

    _out_data = create_output_list();
    _out_data_files = set_output_defaults(_out_data);

    global_in = open_file(global_path.c_str(),"r");
    parse_output_info(&_filenames, global_in, &_out_data_files, _out_data);

    _dmy = make_dmy(&global_param);

    // Initiate params. 初始化参数
    _n_global_inited = FALSE;
    _veglib_gotten = FALSE;
    _soil_gotten = FALSE;
    _vegparam_gotten = FALSE;
    _atmos_inited = FALSE;

    atmos = NULL;
//    _out_data_files = NULL;
//    _out_data = NULL;

    return 0;
}

/* **************************************************************
 *
 * Run at a single cell.
 * 在单个网格上运行。
 *
 * **************************************************************
 */

int run_a_cell(soil_con_struct* soil_con,
               veg_con_struct* veg_con,
               lake_con_struct* lake_con){

    if(!_soil_gotten || !_vegparam_gotten)return -1;

    extern global_param_struct global_param;
    extern veg_lib_struct *veg_lib;
    extern option_struct options;
    extern Error_struct Error;

    option_struct optclone = options;

    int ErrorFlag;
    int cellnum;

    all_vars_struct all_vars;
    veg_hist_struct **veg_hist;

    cellnum = soil_con->gridcel;

    /** 若没有为气象参数分配空间则分配 **/
    if(!_atmos_inited){
        alloc_atmos(global_param.nrecs, &atmos);
        _atmos_inited = TRUE;
    }

    /** Build Gridded Filenames, and Open **/
    make_in_and_outfiles(&_filep, &_filenames, soil_con, _out_data_files);

    /** Write output file headers
     * 写表头 **/
    if (options.PRT_HEADER) write_header(_out_data_files, _out_data, _dmy, global_param);

    if (!options.OUTPUT_FORCE) {
        /** Make Top-level Control Structure **/
        all_vars = make_all_vars(veg_con[0].vegetat_type_num);
        /** allocate memory for the veg_hist_struct **/
        alloc_veg_hist(global_param.nrecs, veg_con[0].vegetat_type_num, &veg_hist);
    } /* !OUTPUT_FORCE */

    /**************************************************
       Initialize Meteological Forcing Values That
       Have not Been Specifically Set
       初始化气象驱动数据
     **************************************************/

    initialize_atmos(atmos, _dmy, _filep.forcing, veg_lib, veg_con, veg_hist,
                     soil_con, _out_data_files, _out_data);

    if (!options.OUTPUT_FORCE) {
        /**************************************************
          Initialize Energy Balance and Snow Variables
          初始化数据
        **************************************************/
        ErrorFlag = initialize_model_state(&all_vars, _dmy[0], &global_param, _filep,
                                           soil_con->gridcel, veg_con[0].vegetat_type_num,
                                           options.Nnode,
                                           atmos[0].air_temp[NR],
                                           soil_con, veg_con, *lake_con);
        if ( ErrorFlag == ERROR ){
            return ERROR;
        }
        /** Update Error Handling Structure **/
        Error.filep = _filep;
        Error.out_data_files = _out_data_files;

        /** Initialize the storage terms in the water and energy balances **/
        /** Sending a negative record number (-global_param.nrecs) to put_data() will accomplish this **/
        ErrorFlag = put_data(&all_vars, &atmos[0], soil_con, veg_con, lake_con, _out_data_files, _out_data, &_save_data, &_dmy[0], -global_param.nrecs);

        /**************************************************
          Run at each time step.
          在每个时间步运行
        **************************************************/
        int startrec = 0;
        for (int rec = startrec ; rec < global_param.nrecs; rec++ ) {

            /* VIC核心函数 能量/水量平衡计算 */
            ErrorFlag = full_energy(cellnum, rec, &atmos[rec], &all_vars, _dmy, &global_param, lake_con, soil_con, veg_con, veg_hist);

            /**************************************************
              Write cell average values for current time step
              导出每个时间步的数据
            **************************************************/
            ErrorFlag = put_data(&all_vars, &atmos[rec], soil_con, veg_con, lake_con, _out_data_files, _out_data, &_save_data, &_dmy[rec], rec);

        } /* End Rec Loop */

    } /* !OUTPUT_FORCE */
    if (!options.OUTPUT_FORCE) {
        free_veg_hist(global_param.nrecs, veg_con[0].vegetat_type_num, &veg_hist);
        free_all_vars(&all_vars, veg_con[0].vegetat_type_num);
    }

    close_files(&_filep, _out_data_files, &_filenames);
    return 0;
}
