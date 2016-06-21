#!/usr/bin/env python
''' VIC Image Driver testing '''
import os
import re

import xarray as xr
import numpy as np
import numpy.testing as npt


def test_image_driver_no_output_file_nans(fnames, domain_file):
    '''
    Test that all VIC image driver output files have the same nan structure as
    the domain file
    '''
    ds_domain = xr.open_dataset(domain_file)
    for fname in fnames:
        ds_output = xr.open_dataset(fname)
        assert_nan_equal(ds_domain, ds_output)


def assert_nan_equal(ds_domain, ds_output):
    """
    test to see if nans in two data arrays are the same in order to check
    domain first dataarray is domain file, second dataarray is a variable in
    the output dataset for image driver
    """

    # get lat/lng coordinate names from domain file mask
    match_vars = ds_domain['mask'].coords

    # check to be sure that mask, lats and lons match between domain file and
    # output file
    for var in match_vars:
        # raise AssertionError if they don't match
        npt.assert_allclose(ds_output[var], ds_domain[var], equal_nan=True)

    # check that nans are occurring in the same place in the arrays
    # check all variables in the dataset
    for da in ds_output.data_vars:
        # get dimensions to reduce DataArray on
        dim_diff = set(ds_domain['mask'].dims).symmetric_difference(
                                                    set(ds_output[da].dims))
        # reduce DataArray
        da_null_reduced = ds_output[da].isnull().all(dim=dim_diff)
        # raise AssertionError if NaNs do not match
        npt.assert_array_equal(da_null_reduced.values,
                               np.isnan(ds_domain['mask']))


def check_multistream_image(fnames):
    '''
    Test the multistream aggregation in the image driver
    '''

    how_dict = {'OUT_ALBEDO': 'max',
                'OUT_SOIL_TEMP': 'min',
                'OUT_PRESSURE': 'sum',
                'OUT_AIR_TEMP': 'first',
                'OUT_SWDOWN': 'mean',
                'OUT_LWDOWN': 'last'}

    streams = {}  # Dictionary to store parsed stream names
    stream_fnames = {}

    for path in fnames:
        # split up the path name to get info about the stream
        resultdir, fname = os.path.split(path)
        pieces = re.split('[_.]', fname)
        stream = '_'.join(pieces[:-2])  # stream name
        freq_n = pieces[2]  # stream frequency n

        stream_fnames[stream] = path

        # set the stream frequency for resample
        if 'NSTEPS' in stream:
            inst_stream = stream
        else:
            if 'NDAYS' in stream:
                streams[stream] = '{}D'.format(freq_n)
            elif 'NHOURS' in stream:
                streams[stream] = '{}H'.format(freq_n)
            elif 'NMINUTES' in stream:
                streams[stream] = '{}min'.format(freq_n)
            elif 'NSECONDS' in stream:
                streams[stream] = '{}S'.format(freq_n)
            else:
                ValueError('stream %s not supported in this test' % stream)

    # Loop over all streams
    instant_ds = xr.open_dataset(stream_fnames[inst_stream])

    # Loop over all streams
    for stream, freq in streams.items():
        agg_ds = xr.open_dataset(stream_fnames[stream])

        # Loop over the variables in the stream
        for key, how in how_dict.items():

            # Resample of the instantaneous data
            expected = instant_ds[key].resample(freq, dim='time', how=how)

            # Get the aggregated values (from VIC)
            actual = agg_ds[key].values

            # Compare the actual and expected (with tolerance)
            npt.assert_array_equal(actual, expected,
                                   err_msg='Variable=%s, freq=%s, how=%s: failed'
                                           ' comparison' % (key, freq, how))
