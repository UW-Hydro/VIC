//
// Created by victi on 16-5-3.
//
/*
 * Testing APIs.
 */
#include <iostream>
#include <vic.h>

using namespace std;

int main(int argc, char** argv){
    extern veg_lib_struct* veg_lib;

    string globalpath = "/home/victi/VIClab/testdata/global.v";


    string soilpath = "/home/victi/VIClab/testdata/soil.v";
    string vegpath = "/home/victi/VIClab/testdata/veg.v";
    int nvegtype = 0;
    int ncell = 17723;

//    _run_all(globalpath);

    soil_con_struct** soilps = new soil_con_struct* [ncell];
    veg_con_struct** vegps = new veg_con_struct* [ncell];

    lake_con_struct lake_con;

    get_global(globalpath);

    nvegtype = get_veg_lib("/home/victi/VIClab/testdata/veglib.LDAS");
    cout<<nvegtype<<endl;

    veg_lib_struct *veglibclone = veg_lib;

    get_soil_params(soilps, ncell, soilpath);
    get_veg_param(vegps, nvegtype, soilps, ncell, vegpath);

    int cell = 8;
    soil_con_struct* si = soilps[cell];
    veg_con_struct* vi = vegps[cell];

    double rs[global_param.nrecs];

    for(int i = 0; i < 8; i++){
        run_a_cell(soilps[i], vegps[i], &lake_con, rs);
        for(int d = 0; d < 10; d++)cout<<rs[d]<<endl;
        cout<<endl;
    }


    destroy_veg_param(vegps, ncell);

    destroy_soil_param(soilps, ncell);

    cout<<"done\n";
    return 0;
}