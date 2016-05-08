/*
 * Testing APIs.
 */
#include <iostream>
#include <vic.h>

using namespace std;

int main(int argc, char** argv){
    string globalpath = "/home/victi/VIClab/testdata/global.v";

    string soilpath = "/home/victi/VIClab/testdata/soil.v";
    string vegpath = "/home/victi/VIClab/testdata/veg.v";
    int nvegtype = 0;
    int ncell = 17723;

    soil_con_struct** soilps = new soil_con_struct* [ncell];
    veg_con_struct** vegps = new veg_con_struct* [ncell];

    lake_con_struct lake_con;

    get_global(globalpath);

    nvegtype = get_veg_lib("/home/victi/VIClab/testdata/veglib.LDAS");
    cout<<nvegtype<<endl;

    get_soil_params(soilps, ncell, soilpath);
    get_veg_param(vegps, nvegtype, soilps, ncell, vegpath);

    int cell = 8;
    soil_con_struct* si = soilps[cell];
    veg_con_struct* vi = vegps[cell];

    double rs[get_timesteps()];

    run_a_cell(soilps[6], vegps[6], &lake_con, rs);
    for(int d = 0; d < 10; d++)cout<<rs[d]<<endl;
    cout<<endl;

    int basin0[] = {0, 5, 6, 12, 13};
    modify_soil_params(soilps, basin0, 5, D2, 0.15);

    run_a_cell(soilps[6], vegps[6], &lake_con, rs);
    for(int d = 0; d < 10; d++)cout<<rs[d]<<endl;
    cout<<endl;

    destroy_veg_param(vegps, ncell);
    destroy_soil_param(soilps, ncell);

    cout<<"done\n";
    return 0;
}