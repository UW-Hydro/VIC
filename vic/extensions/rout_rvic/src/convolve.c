#include <vic_def.h>
#include <rout.h>
#include <vic_driver_image.h>

// function call:
// convolve_new(rout.rout_param.iSources,               /*scalar - number of sources*/
// rout.rout_param.iOutlets,            /*scalar - length of subset*/
// rout.rout_param.iSubsetLength,       /*scalar - length of subset*/
// global_domain.n_nx,
// rout.rout_param.source2outlet_ind,   /*1d array - source to outlet mapping*/
// rout.rout_param.source_y_ind,        /*1d array - source y location*/
// rout.rout_param.source_x_ind,        /*1d array - source x location*/
// rout.rout_param.source_time_offset,  /*1d array - source time offset*/
// rout.rout_param.unit_hydrograph,     /*2d array[times][sources] - unit hydrographs*/
// //rout.rout_param.aggrunin,            /*2d array[ysize][xsize] - vic runoff flux*/
// rout.ring);                          /*2d array[times][outlets] - convolution ring*/


void
convolve_new(const size_t nsources,                    /*scalar - number of sources*/
             const size_t noutlets,                /*scalar - length of subset*/
             const size_t subset_length,           /*scalar - length of subset*/
// const size_t x_size,
             const int   *source2outlet_ind,    /*1d array - source to outlet mapping*/
// const int* source_y_ind,        /*1d array - source y location*/
// const int* source_x_ind,        /*1d array - source x location*/
             const int   *source_time_offset,   /*1d array - source time offset*/
             const double*unit_hydrograph,      /*2d array[times][sources] - unit hydrographs*/
// const double* aggrunin,         /*2d array[ysize][xsize] - vic runoff flux*/
             double      *ring)                 /*2d array[times][outlets] - convolution ring*/
{
    log_info("Doing convolution test!");
    size_t                   s, i, j;              /*counters*/
    // int y, x;
    int                      offset, outlet;      /*2d indicies*/
// int xyind;
    int                      rind, uhind; /*1d indicies*/
    extern out_data_struct **out_data;
    extern rout_struct       rout;

    /*Loop through all sources*/
    for (s = 0; s < nsources; s++) {
        outlet = source2outlet_ind[s];
// y = source_y_ind[s];
// x = source_x_ind[s];
        offset = source_time_offset[s];

        // 1d index location
        // 2d-->1d indexing goes like this:  ind = y*x_size + x
        // xyind = y*x_size + x;

        // printf("source: %4i: , outlet: %4i: ,x: %4i: ,y: %4i ,offset: %4i ,xyind: %4i\n",s,outlet,x,y,offset,xyind);

        /* Do the convolution */
        // i is the position in the unit hydrograph
        // j is the position in the ring
        for (i = 0; i < subset_length; i++) {
            j = i + offset;

            // 1d index locations
            rind = j * noutlets + outlet;
            uhind = i * nsources + s;


            // printf("%4i: ,%4i: ,%4i\n",rind,uhind,xyind);
// ring[rind] += unit_hydrograph[uhind] * aggrunin[xyind];
            // printf("     %4i: ,%4i: ,%4i: ,%4i\n",i,rind,uhind,xyind);
            // ring[i] += unit_hydrograph[uhind] * aggrunin[xyind];      //      rind = noutlets + outlet + j;

            // uhind = s * subset_length + i;
            // printf("new: %4i: ,  %4i: , %4i: ,%4i: ,%4i: ,%4i\n", s, outlet, i, rind, uhind, xyind);

// ring[rind] += unit_hydrograph[uhind] * aggrunin[xyind];
// ring[rind] += unit_hydrograph[uhind];
            // ring[rind] += out_data[rout.rout_param.source_VIC_index[s]][OUT_RUNOFF].data[0] + out_data[rout.rout_param.source_VIC_index[s]][OUT_BASEFLOW].data[0];
            ring[rind] +=
                (unit_hydrograph[uhind] *
                 (out_data[rout.rout_param.source_VIC_index[s]][OUT_RUNOFF].
                  data[0]
                  +
                  out_data[rout.rout_param.source_VIC_index[s]][OUT_BASEFLOW].
                  data[0
                  ]));
        }
    }
}

void
convolve(const size_t nsources,                    /*scalar - number of sources*/
         const size_t noutlets,                    /*scalar - length of subset*/
         const size_t subset_length,               /*scalar - length of subset*/
         const size_t x_size,
         const int   *source2outlet_ind,        /*1d array - source to outlet mapping*/
         const int   *source_y_ind,             /*1d array - source y location*/
         const int   *source_x_ind,             /*1d array - source x location*/
         const int   *source_time_offset,       /*1d array - source time offset*/
         const double*unit_hydrograph,          /*2d array[times][sources] - unit hydrographs*/
         const double*aggrunin,                 /*2d array[ysize][xsize] - vic runoff flux*/
         double      *ring)                     /*2d array[times][outlets] - convolution ring*/
{
    log_info("Doing convolution test!");
    size_t s, i, j;                                /*counters*/
    int    y, x, offset, outlet;                /*2d indicies*/
    int    xyind, rind, uhind;                  /*1d indicies*/

    /*Loop through all sources*/
    for (s = 0; s < nsources; s++) {
        outlet = source2outlet_ind[s];
        y = source_y_ind[s];
        x = source_x_ind[s];
        offset = source_time_offset[s];

        // 1d index location
        // 2d-->1d indexing goes like this:  ind = y*x_size + x
        xyind = y * x_size + x;

        /* Do the convolution */
        // i is the position in the unit hydrograph
        // j is the position in the ring
        for (i = 0; i < subset_length; i++) {
            j = i + offset;

            // 1d index locations
            rind = j * noutlets + outlet;
            uhind = i * nsources + s;
            printf("%4i: ,%4i: ,%4i\n", rind, uhind, xyind);

            ring[rind] += unit_hydrograph[uhind] * aggrunin[xyind];
        }
    }
}
