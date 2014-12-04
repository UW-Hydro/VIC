#include <vic_def.h>
#include <vic_run.h>
#include <vic_driver_classic.h>

void
write_output(out_data_struct      *out_data,
             out_data_file_struct *out_data_files,
             dmy_struct           *dmy,
             int                   rec,
             int                   nveg)
{
    extern global_param_struct global_param;
    extern option_struct       options;

    int                        i;
    int                        index;
    int                        out_dt_sec;
    int                        v;

    out_dt_sec = global_param.out_dt * SECPHOUR;

    /***********************************************
       Change of units for ALMA-compliant output
    ***********************************************/
    if (options.ALMA_OUTPUT) {
        out_data[OUT_BASEFLOW].aggdata[0] /= out_dt_sec;
        out_data[OUT_EVAP].aggdata[0] /= out_dt_sec;
        out_data[OUT_EVAP_BARE].aggdata[0] /= out_dt_sec;
        out_data[OUT_EVAP_CANOP].aggdata[0] /= out_dt_sec;
        out_data[OUT_INFLOW].aggdata[0] /= out_dt_sec;
        out_data[OUT_LAKE_BF_IN].aggdata[0] /= out_dt_sec;
        out_data[OUT_LAKE_BF_IN_V].aggdata[0] /= out_dt_sec;
        out_data[OUT_LAKE_BF_OUT].aggdata[0] /= out_dt_sec;
        out_data[OUT_LAKE_BF_OUT_V].aggdata[0] /= out_dt_sec;
        out_data[OUT_LAKE_CHAN_IN].aggdata[0] /= out_dt_sec;
        out_data[OUT_LAKE_CHAN_IN_V].aggdata[0] /= out_dt_sec;
        out_data[OUT_LAKE_CHAN_OUT].aggdata[0] /= out_dt_sec;
        out_data[OUT_LAKE_CHAN_OUT_V].aggdata[0] /= out_dt_sec;
        out_data[OUT_LAKE_DSTOR].aggdata[0] /= out_dt_sec;
        out_data[OUT_LAKE_DSTOR_V].aggdata[0] /= out_dt_sec;
        out_data[OUT_LAKE_DSWE].aggdata[0] /= out_dt_sec;
        out_data[OUT_LAKE_DSWE_V].aggdata[0] /= out_dt_sec;
        out_data[OUT_LAKE_EVAP].aggdata[0] /= out_dt_sec;
        out_data[OUT_LAKE_EVAP_V].aggdata[0] /= out_dt_sec;
        out_data[OUT_LAKE_ICE_TEMP].aggdata[0] += KELVIN;
        out_data[OUT_LAKE_PREC_V].aggdata[0] /= out_dt_sec;
        out_data[OUT_LAKE_RCHRG].aggdata[0] /= out_dt_sec;
        out_data[OUT_LAKE_RCHRG_V].aggdata[0] /= out_dt_sec;
        out_data[OUT_LAKE_RO_IN].aggdata[0] /= out_dt_sec;
        out_data[OUT_LAKE_RO_IN_V].aggdata[0] /= out_dt_sec;
        out_data[OUT_LAKE_VAPFLX].aggdata[0] /= out_dt_sec;
        out_data[OUT_LAKE_VAPFLX_V].aggdata[0] /= out_dt_sec;
        out_data[OUT_LAKE_SURF_TEMP].aggdata[0] += KELVIN;
        out_data[OUT_PREC].aggdata[0] /= out_dt_sec;
        out_data[OUT_RAINF].aggdata[0] /= out_dt_sec;
        out_data[OUT_REFREEZE].aggdata[0] /= out_dt_sec;
        out_data[OUT_RUNOFF].aggdata[0] /= out_dt_sec;
        out_data[OUT_SNOW_MELT].aggdata[0] /= out_dt_sec;
        out_data[OUT_SNOWF].aggdata[0] /= out_dt_sec;
        out_data[OUT_SUB_BLOWING].aggdata[0] /= out_dt_sec;
        out_data[OUT_SUB_CANOP].aggdata[0] /= out_dt_sec;
        out_data[OUT_SUB_SNOW].aggdata[0] /= out_dt_sec;
        out_data[OUT_SUB_SNOW].aggdata[0] += out_data[OUT_SUB_CANOP].aggdata[0];
        out_data[OUT_SUB_SURFACE].aggdata[0] /= out_dt_sec;
        out_data[OUT_TRANSP_VEG].aggdata[0] /= out_dt_sec;
        out_data[OUT_BARESOILT].aggdata[0] += KELVIN;
        out_data[OUT_SNOW_PACK_TEMP].aggdata[0] += KELVIN;
        out_data[OUT_SNOW_SURF_TEMP].aggdata[0] += KELVIN;
        for (index = 0; index < options.Nlayer; index++) {
            out_data[OUT_SOIL_TEMP].aggdata[index] += KELVIN;
        }
        for (index = 0; index < options.Nnode; index++) {
            out_data[OUT_SOIL_TNODE].aggdata[index] += KELVIN;
            out_data[OUT_SOIL_TNODE_WL].aggdata[index] += KELVIN;
        }
        out_data[OUT_SURF_TEMP].aggdata[0] += KELVIN;
        out_data[OUT_VEGT].aggdata[0] += KELVIN;
        out_data[OUT_FDEPTH].aggdata[0] /= 100;
        out_data[OUT_TDEPTH].aggdata[0] /= 100;
        out_data[OUT_DELTACC].aggdata[0] *= out_dt_sec;
        out_data[OUT_DELTAH].aggdata[0] *= out_dt_sec;
        out_data[OUT_AIR_TEMP].aggdata[0] += KELVIN;
        out_data[OUT_PRESSURE].aggdata[0] *= 1000;
        out_data[OUT_VP].aggdata[0] *= 1000;
        out_data[OUT_VPD].aggdata[0] *= 1000;
    }

    /*************
       Write Data
    *************/
    if (rec >= global_param.skipyear) {
        if (options.BINARY_OUTPUT) {
            for (v = 0; v < N_OUTVAR_TYPES; v++) {
                for (i = 0; i < out_data[v].nelem; i++) {
                    out_data[v].aggdata[i] *= out_data[v].mult;
                }
            }
        }
        write_data(out_data_files, out_data, dmy, global_param.out_dt, nveg);
    }

    // Reset the agg data
    for (v = 0; v < N_OUTVAR_TYPES; v++) {
        for (i = 0; i < out_data[v].nelem; i++) {
            out_data[v].aggdata[i] = 0;
        }
    }
}
