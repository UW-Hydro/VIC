#include <vic_def.h>
#include <rout.h>
#include <vic_driver_image.h>

void convolve(const size_t nsources,               /*scalar - number of sources*/
                const size_t noutlets,             /*scalar - length of subset*/
                const size_t subset_length,        /*scalar - length of subset*/
                const size_t x_size,
                const int* source2outlet_ind,   /*1d array - source to outlet mapping*/
                const int* source_y_ind,        /*1d array - source y location*/
                const int* source_x_ind,        /*1d array - source x location*/
                const int* source_time_offset,  /*1d array - source time offset*/
                const double* unit_hydrograph,  /*2d array[times][sources] - unit hydrographs*/
                const double* aggrunin,         /*2d array[ysize][xsize] - vic runoff flux*/
                double* ring)                   /*2d array[times][outlets] - convolution ring*/
{
    log_info("Doing convolution test!");
    size_t s, i, j;                                /*counters*/
    int y, x, offset, outlet;                   /*2d indicies*/
    int xyind, rind, uhind;                     /*1d indicies*/

    /*Loop through all sources*/
    for (s = 0; s < nsources; s++) {

        outlet = source2outlet_ind[s];
        y = source_y_ind[s];
        x = source_x_ind[s];
        offset = source_time_offset[s];

        //1d index location
        //2d-->1d indexing goes like this:  ind = y*x_size + x
        xyind = y*x_size + x;

        /* Do the convolution */
        // i is the position in the unit hydrograph
        // j is the position in the ring
        for (i = 0; i < subset_length; i++) {
            j = i + offset;

            //1d index locations
            rind = j * noutlets + outlet;
            uhind = i * nsources + s;

            ring[rind] += unit_hydrograph[uhind] * aggrunin[xyind];
        }
    }
}
