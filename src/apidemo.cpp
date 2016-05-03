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

    string soilpath = "/home/victi/VIClab/testdata/soil.v";
    string globalpath = "/home/victi/VIClab/testdata/global.v";
    int nvegtype = 0;
    int ncell = 17723;
    soil_con_struct** soilps = new soil_con_struct* [ncell];;

    get_global(globalpath);

    nvegtype = get_veg_lib("/home/victi/VIClab/testdata/veglib.LDAS");
    cout<<nvegtype<<endl;

    veg_lib_struct *veglibclone = veg_lib;

    get_soil_param(soilps, ncell, soilpath);
    cout<<ncell<<endl;

    destory_soil_param(soilps,ncell);

    cout<<"done\n";
    return 0;
}