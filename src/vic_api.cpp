#include <iostream>
#include <cstring>

#include <vic.h>
#include "vicNl_def.h"

using namespace std;
using namespace vic;

/***************************************************************
 *
 * Internal functions.
 * 内部函数
 *
 ***************************************************************/

double _get_runoff(cell_data_struct** cell, soil_con_struct* soil_con, veg_con_struct* veg_con);
void _make_in_files(filep_struct *filep, filenames_struct *filenames, soil_con_struct *soil);
void _close_in_files(filep_struct *filep, filenames_struct *fnames);

/* **************************************************************
 *
 * 2016-05-03 Zhong Ruida
 *
 * Read in all soil parameter data and create data series.
 * 读入所有土壤参数并生成序列。
 *
 * **************************************************************
 */
int get_veg_param(veg_con_struct ** veg_params, int Nveg_type, soil_con_struct ** soil_params, int ncell, string veg_path) {

    if(!_veglib_gotten || !_soil_gotten) return ERROR;
    if(veg_path == "") veg_path = _filenames.veg;

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

/* **************************************************************
 *
 * 2016-05-03 Zhong Ruida
 *
 * Read in all vegetation parameter data and create data series.
 * 读入所有植被覆盖参数并生成序列。
 *
 * **************************************************************
 */
int get_soil_params(soil_con_struct **soil_params, int ncell, string soil_path){

    if(!_veglib_gotten) return ERROR;
    extern veg_lib_struct* veg_lib;

    if(soil_path == "") soil_path = _filenames.soil;

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

/* **************************************************************
 *
 * 2016-05-03 Zhong Ruida
 *
 * Destroy soil parameter data series.
 * 删除土壤参数数据序列。
 *
 * **************************************************************
 */
void destroy_soil_param(soil_con_struct **soil_params, int ncell){
    if(!_soil_gotten) return;
    for(int i = 0; i < ncell; i++){
        if(soil_params[i] == NULL) continue;
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

/* **************************************************************
 *
 * 2016-05-03 Zhong Ruida
 *
 * Destroy vegetation parameter data series.
 * 删除植被覆盖数据序列。
 *
 * **************************************************************
 */
void destroy_veg_param(veg_con_struct **veg_params, int ncell){
    if(!_vegparam_gotten) return;

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
 * 2016-05-03 Zhong Ruida
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
    if(veglib_path == "") veglib_path = _filenames.veglib;

    extern veg_lib_struct* veg_lib;
    int n_veg_type = 0;

    FILE * veglib_in;
    veglib_in = open_file(veglib_path.c_str(), "r");
    if(veglib_in == NULL) return ERROR;

    veg_lib = read_veglib(veglib_in, &n_veg_type);

    fclose(veglib_in);
    _veglib_gotten = TRUE;
    return n_veg_type;
}

/* **************************************************************
 *
 * 2016-05-03 Zhong Ruida
 *
 * Get global parameters and initiate.
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

    return 0;
}

/* **************************************************************
 *
 * 2016-05-05 Zhong Ruida
 *
 * Run at a single cell.
 * 在单个网格上运行。
 *
 * **************************************************************
 */
int run_a_cell(soil_con_struct* soil_con,
               veg_con_struct* veg_con,
               lake_con_struct* lake_con, double* out_runoffs) {

    if(!_soil_gotten || !_vegparam_gotten)return -1;

    extern global_param_struct global_param;
    extern veg_lib_struct *veg_lib;
    extern option_struct options;
    extern Error_struct Error;

    int ErrorFlag;
    int cellnum;

    all_vars_struct all_vars;
    veg_hist_struct **veg_hist;

    cellnum = soil_con->gridcel;

    /** Allocate space for atmos structures if not allocated **/
    /** 若没有为气象参数分配空间则分配 **/
    if(!_atmos_inited){
        alloc_atmos(global_param.nrecs, &atmos);
        _atmos_inited = TRUE;
    }

    /** Build Gridded Filenames, and Open **/
    /** 生成各网格的输入和输出文件名(若无提供产流量待写数组)并打开文件 **/
    if(out_runoffs != NULL) _make_in_files(&_filep, &_filenames, soil_con);
    else {
        make_in_and_outfiles(&_filep, &_filenames, soil_con, _out_data_files);
        if (options.PRT_HEADER) write_header(_out_data_files, _out_data, _dmy, global_param);
    }

    /** Make Top-level Control Structure **/
    all_vars = make_all_vars(veg_con[0].vegetat_type_num);

    /** allocate memory for the veg_hist_struct **/
    /** 为植被时间序列结构分配内存空间 **/
    alloc_veg_hist(global_param.nrecs, veg_con[0].vegetat_type_num, &veg_hist);

    /**************************************************
        Initialize Variables
        初始化变量
    **************************************************/
    initialize_atmos(atmos, _dmy, _filep.forcing, veg_lib, veg_con, veg_hist,
                     soil_con, _out_data_files, _out_data);

    ErrorFlag = initialize_model_state(&all_vars, _dmy[0], &global_param, _filep,
                                       soil_con->gridcel, veg_con[0].vegetat_type_num,
                                       options.Nnode,
                                       atmos[0].air_temp[NR],
                                       soil_con, veg_con, *lake_con);
    if ( ErrorFlag == ERROR ){
        return ERROR;
    }
    /** Update Error Handling Structure **/
    /** 更新错误控制结构 **/
    Error.filep = _filep;
    Error.out_data_files = _out_data_files;

    /** Initialize the storage terms in the water and energy balances **/
    /** Sending a negative record number (-global_param.nrecs) to put_data() will accomplish this **/
    if(out_runoffs == NULL)
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
          Write cell average values for current time step.
          导出每个时间步的数据。若提供的待写入的产流量数组为NULL则为输出到文件。
        **************************************************/
        if(out_runoffs == NULL)ErrorFlag = put_data(&all_vars, &atmos[rec], soil_con, veg_con, lake_con, _out_data_files, _out_data, &_save_data, &_dmy[rec], rec);
        else out_runoffs[rec] = _get_runoff(all_vars.cell, soil_con, veg_con);
    }

    free_veg_hist(global_param.nrecs, veg_con[0].vegetat_type_num, &veg_hist);
    free_all_vars(&all_vars, veg_con[0].vegetat_type_num);

    if(out_runoffs != NULL) _close_in_files(&_filep, &_filenames);
    else close_files(&_filep, _out_data_files, &_filenames);
    return 0;
}

/* **************************************************************
 *
 * 2016-05-06 Zhong Ruida
 *
 * Get runoff value from each time step without other output data
 * replacing function put_data.
 * 在每个时间步只获取产流量数据。用于取代put_data。
 *
 * (Situation with Lake/Wetland Parameter File input is not supported
 * at preasent.)
 * (目前尚不支持带有湖泊/湿地参数输入的情况。)
 *
 * **************************************************************
 */
double _get_runoff(cell_data_struct** cell,
                   soil_con_struct* soil_con,
                   veg_con_struct* veg_con) {

    double Cv;
    double treeAdj;
    double areaFract;
    double runoff = 0.0;

    int Nbands = options.SNOW_BAND;
    int Nvegs = veg_con[0].vegetat_type_num;

    for(int band = 0; band < Nbands; band++){
        areaFract = soil_con->AreaFract[band];
        if(!soil_con->AboveTreeLine[band]) treeAdj = 1.0;
        else{
            Cv = 0;
            for(int veg = 0; veg < Nvegs; veg++)
                if(veg_lib[veg_con[veg].veg_class].overstory ) Cv += veg_con[veg].Cv;
            treeAdj = Cv / (1.0 - Cv);
        }
        for(int veg = 0; veg <= Nvegs; veg++){
            Cv = veg_con[veg].Cv;
            double i_runoff = cell[veg][band].baseflow + cell[veg][band].runoff;
            runoff += i_runoff * Cv * areaFract * treeAdj;
        }
    }
    return runoff;
}

/* **************************************************************
 *
 * 2016-05-07 Zhong Ruida
 *
 * Modify soil parameters to be calibrated, including INFILT, DS,
 * DS_MAX, WS and D1~D3 by a series of cell ids.
 * 根据一系列网格点编号修改率定参数,包括INFILT, DS, DS_MAX, WS和D1~D3。
 *
 *
 * **************************************************************
 */

void modify_soil_params(soil_con_struct** soil_params,
                        int* cellids,
                        int cellnum,
                        int param_type,
                        double value) {

    extern option_struct options;

    if(param_type == INFILT) {
        for (int i = 0; i < cellnum; i++)
            soil_params[cellids[i]]->b_infilt = value;
        return;
    }
    if(param_type == DS) {
        for (int i = 0; i < cellnum; i++)
            soil_params[cellids[i]]->Ds = value;
        return;
    }
    if(param_type == DS_MAX) {
        for (int i = 0; i < cellnum; i++)
            soil_params[cellids[i]]->Dsmax = value;
        return;
    }
    if(param_type == WS) {
        for (int i = 0; i < cellnum; i++)
            soil_params[cellids[i]]->Ws = value;
        return;
    }
    if(param_type == D1) {
        for (int i = 0; i < cellnum; i++)
            soil_params[cellids[i]]->depth[0] = value;
        return;
    }
    if(param_type == D2 && options.Nlayer >= 2) {
        for (int i = 0; i < cellnum; i++)
            soil_params[cellids[i]]->depth[1] = value;
        return;
    }
    if(param_type == D3 && options.Nlayer >= 3) {
        for (int i = 0; i < cellnum; i++)
            soil_params[cellids[i]]->depth[2] = value;
        return;
    }
}

/***************************************************************
 *
 * Internal functions.
 * 以下为内部函数
 *
 ***************************************************************/

/* **************************************************************
 *
 * 2016-05-07 Zhong Ruida
 *
 * Builds the files names for input and output of grided data files and open.
 * 生成气象驱动数据的输入文件的文件名并打开文件。
 *
 * (Code copied from make_in_and_outfile)
 *
 * **************************************************************
 */
void _make_in_files(filep_struct *filep,
                    filenames_struct *filenames,
                    soil_con_struct *soil)
{

    extern option_struct    options;
    extern param_set_struct param_set;

    char   latchar[20], lngchar[20], junk[6];

    sprintf(junk, "%%.%if", options.GRID_DECIMAL);
    sprintf(latchar, junk, soil->lat);
    sprintf(lngchar, junk, soil->lng);

    /********************************
    Input Forcing Files
    ********************************/

    strcpy(filenames->forcing[0], filenames->f_path_pfx[0]);
    strcat(filenames->forcing[0], latchar);
    strcat(filenames->forcing[0], "_");
    strcat(filenames->forcing[0], lngchar);
    if(param_set.FORCE_FORMAT[0] == BINARY)
        filep->forcing[0] = open_file(filenames->forcing[0], "rb");
    else
        filep->forcing[0] = open_file(filenames->forcing[0], "r");

    filep->forcing[1] = NULL;
    if(strcasecmp(filenames->f_path_pfx[1],"MISSING")!=0) {
        strcpy(filenames->forcing[1], filenames->f_path_pfx[1]);
        strcat(filenames->forcing[1], latchar);
        strcat(filenames->forcing[1], "_");
        strcat(filenames->forcing[1], lngchar);
        if(param_set.FORCE_FORMAT[0] == BINARY)
            filep->forcing[1] = open_file(filenames->forcing[1], "rb");
        else
            filep->forcing[1] = open_file(filenames->forcing[1], "r");
    }
}
/* **************************************************************
 *
 * 2016-05-07 Zhong Ruida
 *
 * Close input files.
 * 关闭气象驱动数据输入文件。
 *
 * (Code copied from close_files)
 *
 * **************************************************************
 */
void _close_in_files(filep_struct *filep,
                     filenames_struct *fnames) {
    extern option_struct options;

    fclose(filep->forcing[0]);
    if (options.COMPRESS) compress_files(fnames->forcing[0]);
    if (filep->forcing[1] != NULL) {
        fclose(filep->forcing[1]);
        if (options.COMPRESS) compress_files(fnames->forcing[1]);
    }
}