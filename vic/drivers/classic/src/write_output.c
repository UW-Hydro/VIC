/******************************************************************************
 * @section DESCRIPTION
 *
 * Write output data.
 *****************************************************************************/

#include <vic_driver_classic.h>

/******************************************************************************
 * @brief    Write output data and convert units if necessary.
 *****************************************************************************/
void
write_output(stream_struct **streams,
             dmy_struct     *dmy)
{
    extern option_struct options;

    size_t               stream_idx;

    // Write data
    for (stream_idx = 0; stream_idx < options.Noutstreams; stream_idx++) {
        if (raise_alarm(&(*streams)[stream_idx].agg_alarm, dmy)) {
            write_data(&((*streams)[stream_idx]));
            reset_stream((&(*streams)[stream_idx]), dmy);
        }
    }
}
