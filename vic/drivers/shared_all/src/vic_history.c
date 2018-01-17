/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine contains routines for computing and saving VIC history files.
 *
 * @section LICENSE
 *
 * The Variable Infiltration Capacity (VIC) macroscale hydrological model
 * Copyright (C) 2016 The Computational Hydrology Group, Department of Civil
 * and Environmental Engineering, University of Washington.
 *
 * The VIC model is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 *****************************************************************************/

#include <vic_driver_shared_all.h>

/******************************************************************************
 * @brief    This routine creates the list of output data.
 *****************************************************************************/
void
alloc_out_data(size_t    ngridcells,
               double ***out_data)
{
    extern metadata_struct out_metadata[N_OUTVAR_TYPES];

    size_t                 i;
    size_t                 j;

    for (i = 0; i < ngridcells; i++) {
        out_data[i] = calloc(N_OUTVAR_TYPES, sizeof(*(out_data[i])));
        check_alloc_status(out_data[i], "Memory allocation error.");
        // Allocate space for data
        for (j = 0; j < N_OUTVAR_TYPES; j++) {
            out_data[i][j] =
                calloc(out_metadata[j].nelem, sizeof(*(out_data[i][j])));
            check_alloc_status(out_data[i][j], "Memory allocation error.");
        }
    }
}

/******************************************************************************
 * @brief    This routine creates the list of output streams.
 *****************************************************************************/
void
setup_stream(stream_struct *stream,
             size_t         nvars,
             size_t         ngridcells)
{
    size_t     i;
    int        default_n = 1;
    dmy_struct dmy_junk;

    // Set stream scalars
    stream->nvars = nvars;
    stream->ngridcells = ngridcells;
    stream->file_format = UNSET_FILE_FORMAT;
    stream->compress = false;

    // Initialize dmy_junk - this step is to avoid time-related error caused
    // by junk dmy; the date set here does not matter and will be overwritten
    // later
    dmy_junk.day = 1;
    dmy_junk.day_in_year = 1;
    dmy_junk.month = 12;
    dmy_junk.year = 1900;
    dmy_junk.dayseconds = 0;

    // Set default agg alarm
    set_alarm(&dmy_junk, FREQ_NDAYS, &default_n, &(stream->agg_alarm));

    // Set default write alarm
    set_alarm(&dmy_junk, FREQ_END, &default_n, &(stream->write_alarm));

    // Allocate stream members of shape [nvars]
    stream->varid = calloc(nvars, sizeof(*(stream->varid)));
    check_alloc_status(stream->varid, "Memory allocation error.");

    stream->aggtype = calloc(nvars, sizeof(*(stream->aggtype)));
    check_alloc_status(stream->aggtype, "Memory allocation error.");

    stream->type = calloc(nvars, sizeof(*(stream->type)));
    check_alloc_status(stream->type, "Memory allocation error.");

    stream->mult = calloc(nvars, sizeof(*(stream->mult)));
    check_alloc_status(stream->mult, "Memory allocation error.");

    // Question: do we have to dynamically allocate the length of each string
    stream->format = calloc(nvars, sizeof(*(stream->format)));
    check_alloc_status(stream->format, "Memory allocation error.");

    for (i = 0; i < nvars; i++) {
        stream->format[i] = calloc(MAXSTRING, sizeof(*(stream->format[i])));
        check_alloc_status(stream->format[i], "Memory allocation error.");
    }
    // Initialize some of the stream members
    // these will be overwritten in set_output_var
    for (i = 0; i < nvars; i++) {
        stream->type[i] = OUT_TYPE_DEFAULT;
        stream->mult[i] = OUT_MULT_DEFAULT;
        stream->aggtype[i] = AGG_TYPE_DEFAULT;
    }
}

/******************************************************************************
 * @brief    This routine validates the streams.
 *****************************************************************************/
void
validate_streams(stream_struct **streams)
{
    extern option_struct options;

    size_t               streamnum;

    // validate stream settings
    for (streamnum = 0; streamnum < options.Noutstreams; streamnum++) {
        if ((*streams)[streamnum].ngridcells < 1) {
            log_err("Number of gridcells in stream is less than 1");
        }
        if ((*streams)[streamnum].nvars < 1) {
            log_err("Number of variables in stream is less than 1");
        }
        if (strcasecmp("", ((*streams)[streamnum].prefix)) == 0) {
            log_err("Stream prefix not set");
        }
        if ((*streams)[streamnum].file_format == UNSET_FILE_FORMAT) {
            log_err("Stream file_format not set");
        }
        if ((*streams)[streamnum].type == NULL) {
            log_err("Stream type array not allocated");
        }
        if ((*streams)[streamnum].mult == NULL) {
            log_err("Stream mult array not allocated");
        }
        if ((*streams)[streamnum].varid == NULL) {
            log_err("Stream varid array not allocated");
        }
        if ((*streams)[streamnum].aggtype == NULL) {
            log_err("Stream aggtype array not allocated");
        }
        if ((*streams)[streamnum].aggdata == NULL) {
            log_err("Stream agg_data array not allocated");
        }
    }
}

/******************************************************************************
 * @brief   This routine allocates memory for the stream aggdata array.  The
            shape of this array is [ngridcells, nvars, nelems, nbins].
 *****************************************************************************/
void
alloc_aggdata(stream_struct *stream)
{
    extern metadata_struct out_metadata[N_OUTVAR_TYPES];

    size_t                 i;
    size_t                 j;
    size_t                 k;
    size_t                 nelem;

    stream->aggdata = calloc(stream->ngridcells, sizeof(*(stream->aggdata)));
    check_alloc_status(stream->aggdata, "Memory allocation error.");

    for (i = 0; i < stream->ngridcells; i++) {
        stream->aggdata[i] =
            calloc(stream->nvars, sizeof(*(stream->aggdata[i])));
        check_alloc_status(stream->aggdata[i], "Memory allocation error.");
        for (j = 0; j < stream->nvars; j++) {
            nelem = out_metadata[stream->varid[j]].nelem;
            stream->aggdata[i][j] =
                calloc(nelem, sizeof(*(stream->aggdata[i][j])));
            check_alloc_status(stream->aggdata[i][j],
                               "Memory allocation error.");

            for (k = 0; k < nelem; k++) {
                // TODO: Also allocate for nbins, for now just setting to size 1
                stream->aggdata[i][j][k] =
                    calloc(1, sizeof(*(stream->aggdata[i][j][k])));
                check_alloc_status(stream->aggdata[i][j][k],
                                   "Memory allocation error.");
            }
        }
    }
}

/******************************************************************************
 * @brief   This routine resets an output stream
 *****************************************************************************/
void
reset_stream(stream_struct *stream,
             dmy_struct    *dmy_current)
{
    extern metadata_struct out_metadata[N_OUTVAR_TYPES];

    size_t                 i;
    size_t                 j;
    size_t                 k;
    size_t                 varid;

    // Reset alarm to next agg period
    reset_alarm(&(stream->agg_alarm), dmy_current);

    // Set aggdata to zero
    for (i = 0; i < stream->ngridcells; i++) {
        for (j = 0; j < stream->nvars; j++) {
            varid = stream->varid[j];
            for (k = 0; k < out_metadata[varid].nelem; k++) {
                stream->aggdata[i][j][k][0] = 0.;
            }
        }
    }
}

/******************************************************************************
 * @brief   This routine sets the default aggregation type for variables in an
            output stream
 *****************************************************************************/
unsigned int
get_default_outvar_aggtype(unsigned int varid)
{
    unsigned int agg_type;

    switch (varid) {
    // AGG_TYPE_END
    case OUT_ASAT:
    case OUT_LAKE_AREA_FRAC:
    case OUT_LAKE_DEPTH:
    case OUT_LAKE_ICE:
    case OUT_LAKE_ICE_FRACT:
    case OUT_LAKE_ICE_HEIGHT:
    case OUT_LAKE_MOIST:
    case OUT_LAKE_SURF_AREA:
    case OUT_LAKE_SWE:
    case OUT_LAKE_SWE_V:
    case OUT_LAKE_VOLUME:
    case OUT_ROOTMOIST:
    case OUT_SMFROZFRAC:
    case OUT_SMLIQFRAC:
    case OUT_SNOW_CANOPY:
    case OUT_SNOW_COVER:
    case OUT_SNOW_DEPTH:
    case OUT_SOIL_ICE:
    case OUT_SOIL_LIQ:
    case OUT_SOIL_MOIST:
    case OUT_SOIL_WET:
    case OUT_SURFSTOR:
    case OUT_SURF_FROST_FRAC:
    case OUT_SWE:
    case OUT_WDEW:
    case OUT_ZWT:
    case OUT_ZWT_LUMPED:
    case OUT_SNOW_CANOPY_BAND:
    case OUT_SNOW_COVER_BAND:
    case OUT_SNOW_DEPTH_BAND:
    case OUT_SWE_BAND:
        agg_type = AGG_TYPE_END;
        break;
    // AGG_TYPE_SUM
    case OUT_BASEFLOW:
    case OUT_DELINTERCEPT:
    case OUT_DELSOILMOIST:
    case OUT_DELSWE:
    case OUT_DELSURFSTOR:
    case OUT_EVAP:
    case OUT_EVAP_BARE:
    case OUT_EVAP_CANOP:
    case OUT_INFLOW:
    case OUT_LAKE_BF_IN:
    case OUT_LAKE_BF_IN_V:
    case OUT_LAKE_BF_OUT:
    case OUT_LAKE_BF_OUT_V:
    case OUT_LAKE_CHAN_IN:
    case OUT_LAKE_CHAN_IN_V:
    case OUT_LAKE_CHAN_OUT:
    case OUT_LAKE_CHAN_OUT_V:
    case OUT_LAKE_DSTOR:
    case OUT_LAKE_DSTOR_V:
    case OUT_LAKE_DSWE:
    case OUT_LAKE_DSWE_V:
    case OUT_LAKE_EVAP:
    case OUT_LAKE_EVAP_V:
    case OUT_LAKE_PREC_V:
    case OUT_LAKE_RCHRG:
    case OUT_LAKE_RCHRG_V:
    case OUT_LAKE_RO_IN:
    case OUT_LAKE_RO_IN_V:
    case OUT_LAKE_VAPFLX:
    case OUT_LAKE_VAPFLX_V:
    case OUT_PET:
    case OUT_PREC:
    case OUT_RAINF:
    case OUT_REFREEZE:
    case OUT_RUNOFF:
    case OUT_SNOWF:
    case OUT_SUB_BLOWING:
    case OUT_SUB_CANOP:
    case OUT_SUB_SNOW:
    case OUT_SUB_SURFACE:
    case OUT_TRANSP_VEG:
    case OUT_DELTACC_BAND:
    case OUT_SNOW_MELT:
    case OUT_SNOWT_FBFLAG:
    case OUT_SOILT_FBFLAG:
    case OUT_SURFT_FBFLAG:
    case OUT_TCAN_FBFLAG:
    case OUT_TFOL_FBFLAG:
        agg_type = AGG_TYPE_SUM;
        break;
    default:
        agg_type = AGG_TYPE_AVG;
    }
    return agg_type;
}

/******************************************************************************
 * @brief    This routine updates the output information for a given output
 *           variable.
 *****************************************************************************/
void
set_output_var(stream_struct     *stream,
               char              *varname,
               size_t             varnum,
               char              *format,
               unsigned short int type,
               double             mult,
               unsigned short int aggtype)
{
    extern metadata_struct out_metadata[N_OUTVAR_TYPES];

    int                    varid;
    int                    found = false;

    if (varnum >= stream->nvars) {
        log_err("Invalid varnum %zu, must be less than the number of variables "
                "in the stream %zu", varnum, stream->nvars);
    }
    // Find the output varid by looping through out_metadata, comparing to varname
    for (varid = 0; varid < N_OUTVAR_TYPES; varid++) {
        if (strcmp(out_metadata[varid].varname, varname) == 0) {
            found = true;
            break;
        }
    }
    if (!found) {
        log_err("set_output_var: \"%s\" was not found in the list of "
                "supported output variable names.  Please use the exact name "
                "listed in vic_driver_shared.h.", varname);
    }
    // Set stream members
    stream->varid[varnum] = varid;
    // Format (ASCII only)
    if ((strcmp(format, "*") != 0) || (strcmp(format, "") != 0)) {
        strcpy(stream->format[varnum], format);
    }
    else {
        strcpy(stream->format[varnum], "%.4f");
    }
    // Output type (BINARY and netCDF)
    if (type != OUT_TYPE_DEFAULT) {
        stream->type[varnum] = type;
    }
    else {
        stream->type[varnum] = OUT_TYPE_FLOAT;
    }
    // Mult (BINARY and netCDF)
    if (mult != OUT_MULT_DEFAULT) {
        stream->mult[varnum] = mult;
    }
    else {
        stream->mult[varnum] = 1.;
    }
    // Aggregation type
    if (aggtype != AGG_TYPE_DEFAULT) {
        stream->aggtype[varnum] = aggtype;
    }
    else {
        stream->aggtype[varnum] = get_default_outvar_aggtype(varid);
    }
}

/******************************************************************************
 * @brief    This routine frees the memory in the streams array.
 *****************************************************************************/
void
free_streams(stream_struct **streams)
{
    extern option_struct   options;
    extern metadata_struct out_metadata[N_OUTVAR_TYPES];

    size_t                 streamnum;
    size_t                 i;
    size_t                 j;
    size_t                 k;
    size_t                 varid;

    // free output streams
    for (streamnum = 0; streamnum < options.Noutstreams; streamnum++) {
        // Free aggdata first
        for (i = 0; i < (*streams)[streamnum].ngridcells; i++) {
            for (j = 0; j < (*streams)[streamnum].nvars; j++) {
                varid = (*streams)[streamnum].varid[j];
                for (k = 0; k < out_metadata[varid].nelem; k++) {
                    free((*streams)[streamnum].aggdata[i][j][k]);
                }
                free((*streams)[streamnum].aggdata[i][j]);
            }
            free((*streams)[streamnum].aggdata[i]);
        }
        for (j = 0; j < (*streams)[streamnum].nvars; j++) {
            free((*streams)[streamnum].format[j]);
        }
        free((*streams)[streamnum].aggdata);
        // free remaining arrays
        free((*streams)[streamnum].type);
        free((*streams)[streamnum].mult);
        free((*streams)[streamnum].format);
        free((*streams)[streamnum].varid);
        free((*streams)[streamnum].aggtype);
    }
    free(*streams);
}

/******************************************************************************
 * @brief    This routine frees the memory in the out_data array.
 *****************************************************************************/
void
free_out_data(size_t    ngridcells,
              double ***out_data)
{
    size_t i;
    size_t j;


    if (out_data == NULL) {
        return;
    }

    for (i = 0; i < ngridcells; i++) {
        for (j = 0; j < N_OUTVAR_TYPES; j++) {
            free(out_data[i][j]);
        }
        free(out_data[i]);
    }

    free(out_data);
}
