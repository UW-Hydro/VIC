/******************************************************************************
 * @section DESCRIPTION
 *
 * Initialize output structures.
 *****************************************************************************/

#include <vic_driver_shared_all.h>

/******************************************************************************
 * @brief    Perform temporal aggregation on stream data
 *****************************************************************************/
void
agg_stream_data(stream_struct *stream,
                dmy_struct    *dmy_current,
                double      ***out_data)
{
    extern metadata_struct out_metadata[N_OUTVAR_TYPES];

    alarm_struct          *alarm;
    size_t                 i;
    size_t                 j;
    size_t                 k;
    size_t                 nelem;
    unsigned int           varid;
    bool                   alarm_now;

    alarm = &(stream->agg_alarm);
    alarm->count++;
    alarm_now = raise_alarm(alarm, dmy_current);

    if (alarm->count == 1) {
        stream->time_bounds[0] = *dmy_current;
    }

    if (alarm_now) {
        stream->time_bounds[1] = *dmy_current;
    }

    for (i = 0; i < stream->ngridcells; i++) {
        for (j = 0; j < stream->nvars; j++) {
            varid = stream->varid[j];
            nelem = out_metadata[varid].nelem;

            // Instantaneous at the beginning of the period
            if ((stream->aggtype[j] == AGG_TYPE_END) && (alarm_now)) {
                for (k = 0; k < nelem; k++) {
                    stream->aggdata[i][j][k][0] = out_data[i][varid][k];
                }
            }
            // Instantaneous at the end of the period
            else if ((stream->aggtype[j] == AGG_TYPE_BEG) &&
                     (alarm->count == 1)) {
                for (k = 0; k < nelem; k++) {
                    stream->aggdata[i][j][k][0] = out_data[i][varid][k];
                }
            }
            // Sum over the period
            else if ((stream->aggtype[j] == AGG_TYPE_SUM) ||
                     (stream->aggtype[j] == AGG_TYPE_AVG)) {
                for (k = 0; k < nelem; k++) {
                    stream->aggdata[i][j][k][0] += out_data[i][varid][k];
                }
            }
            // Maximum over the period
            else if (stream->aggtype[j] == AGG_TYPE_MAX) {
                for (k = 0; k < nelem; k++) {
                    stream->aggdata[i][j][k][0] =
                        max(stream->aggdata[i][j][k][0], out_data[i][varid][k]);
                }
            }
            // Minimum over the period
            else if (stream->aggtype[j] == AGG_TYPE_MIN) {
                for (k = 0; k < nelem; k++) {
                    stream->aggdata[i][j][k][0] =
                        min(stream->aggdata[i][j][k][0], out_data[i][varid][k]);
                }
            }
            // Average over the period if counter is full
            if ((stream->aggtype[j] == AGG_TYPE_AVG) && (alarm_now)) {
                for (k = 0; k < nelem; k++) {
                    stream->aggdata[i][j][k][0] /= (double) alarm->count;
                }
            }
        }
    }
}
