#!/usr/bin/env python
''' VIC Image Driver testing '''
import os
import re

import xarray as xr
import pandas as pd
import numpy as np
import numpy.testing as npt
import glob
import warnings


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
    mask_dims = set(ds_domain['mask'].dims)
    for var, da in ds_output.data_vars.items():
        # skip variables that don't have mask_dims (e.g. time_bnds)
        if not all(d in da.dims for d in mask_dims):
            continue
        # get dimensions to reduce DataArray on
        dim_diff = mask_dims.symmetric_difference(set(da.dims))
        # reduce DataArray
        da_null_reduced = da.isnull().all(dim=dim_diff)
        # raise AssertionError if NaNs do not match
        npt.assert_array_equal(da_null_reduced.values,
                               np.isnan(ds_domain['mask']))


def check_multistream_image(fnames):
    '''
    Test the multistream aggregation in the image driver
    '''

    def get_pandas_how_from_cell_method(cell_method):
        if cell_method == 'time: end':
            how = 'last'
        elif cell_method == 'time: beg':
            how = 'first'
        elif cell_method == 'time: mean':
            how = np.mean
        elif cell_method == 'time: minimum':
            how = 'min'
        elif cell_method == 'time: maximum':
            how = 'max'
        elif cell_method == 'time: sum':
            how = 'sum'
        else:
            raise ValueError('Unknown cell method argument: %s', cell_method)
        return how

    def reindex_xr_obj_timedim(obj, freq):
        # Here we're basically rounding the timestamp in the time index to even
        # values.
        new = pd.date_range(obj.time[0].values, freq=freq,
                            periods=len(obj.time))
        return instant_ds.reindex({'time': new}, method='nearest')

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

    # Open the instantaneous stream
    instant_ds = xr.open_dataset(stream_fnames[inst_stream])
    instant_ds = reindex_xr_obj_timedim(instant_ds, '1H')  # TODO infer freq

    # Loop over all streams
    for stream, freq in streams.items():
        print(stream, freq)
        agg_ds = xr.open_dataset(stream_fnames[stream])
        agg_ds = reindex_xr_obj_timedim(agg_ds, freq)

        # Loop over the variables in the stream
        for key, agg_da in agg_ds.data_vars.items():
            how = get_pandas_how_from_cell_method(agg_da.attrs['cell_methods'])
            # Resample of the instantaneous data
            expected = instant_ds[key].resample(freq, dim='time', how=how,
                                                label='left', closed='left')

            # Get the aggregated values (from VIC)
            actual = agg_da

            # Compare the actual and expected (with tolerance)
            try:
                npt.assert_array_equal(actual.values, expected.values)
            except AssertionError as e:
                print('Variable=%s, freq=%s, how=%s: failed comparison' %
                      (key, freq, how))
                print('actual=%s\nexpected=%s' % (actual, expected))
                print(np.abs(actual - expected).max())
                raise e


def setup_subdirs_and_fill_in_global_param_mpi_test(
        s, list_n_proc, result_basedir, state_basedir, test_data_dir):
    ''' Fill in global parameter output directories for multiple runs for mpi
        testing, image driver

    Parameters
    ----------
    s: <string.Template>
        Template of the global param file to be filled in
    list_n_proc: <list>
        A list of number of processors to run and compare
    result_basedir: <str>
        Base directory of output fluxes results; runs with different number of
        processors are output to subdirectories under the base directory
    state_basedir: <str>
        Base directory of output state results; runs with different number of
        processors are output to subdirectories under the base directory
    test_data_dir: <str>
        Base directory of test data

    Returns
    ----------
    list_global_param: <list>
        A list of global parameter strings to be run with parameters filled in

    Require
    ----------
    os
    '''

    list_global_param = []
    for j, n_proc in enumerate(list_n_proc):
        # Set up subdirectories for results and states
        result_dir = os.path.join(result_basedir,
                                  'processors_{}'.format(n_proc))
        state_dir = os.path.join(state_basedir,
                                 'processors_{}'.format(n_proc))
        os.makedirs(result_dir, exist_ok=True)
        os.makedirs(state_dir, exist_ok=True)

        # Fill in global parameter options
        list_global_param.append(s.safe_substitute(test_data_dir=test_data_dir,
                                                   result_dir=result_dir,
                                                   state_dir=state_dir))

    return(list_global_param)


def check_mpi_fluxes(result_basedir, list_n_proc):
    ''' Check whether all the fluxes are the same with different number of
        processors, image driver

    Parameters
    ----------
    result_basedir: <str>
        Base directory of output fluxes results; runs with different number of
        processors are output to subdirectories under the base directory
    list_n_proc: <list>
        A list of number of processors to run and compare

    Require
    ----------
    os
    glob
    numpy
    warnings
    '''

    # Read the first run - as base
    n_proc = list_n_proc[0]
    result_dir = os.path.join(
        result_basedir,
        'processors_{}'.format(n_proc))
    if len(glob.glob(os.path.join(result_dir, '*.nc'))) > 1:
        warnings.warn(
            'More than one netCDF file found under directory {}'.
            format(result_dir))
    fname = glob.glob(os.path.join(result_dir, '*.nc'))[0]
    ds_first_run = xr.open_dataset(fname)

    # Loop over all rest runs and compare fluxes with the base run
    for i, n_proc in enumerate(list_n_proc):
        # Skip the first run
        if i == 0:
            continue
        # Read flux results for this run
        result_dir = os.path.join(
            result_basedir,
            'processors_{}'.format(n_proc))
        if len(glob.glob(os.path.join(result_dir, '*.nc'))) > 1:
            warnings.warn('More than one netCDF file found under '
                          'directory {}'.format(result_dir))
        fname = glob.glob(os.path.join(result_dir, '*.nc'))[0]
        ds_current_run = xr.open_dataset(fname)
        # Compare current run with base run
        for var in ds_first_run.data_vars:
            npt.assert_array_equal(ds_current_run[var].values,
                                   ds_first_run[var].values,
                                   err_msg='Fluxes are not an exact match')


def check_mpi_states(state_basedir, list_n_proc):
    ''' Check whether all the output states are the same with different number
        of processors, image driver

    Parameters
    ----------
    state_basedir: <str>
        Base directory of output states; runs with different number of
        processors are output to subdirectories under the base directory
    list_n_proc: <list>
        A list of number of processors to run and compare

    Require
    ----------
    os
    glob
    numpy
    warnings
    '''

    # Read the first run - as base
    n_proc = list_n_proc[0]
    state_dir = os.path.join(
        state_basedir,
        'processors_{}'.format(n_proc))
    if len(glob.glob(os.path.join(state_dir, '*.nc'))) > 1:
        warnings.warn('More than one netCDF file found under '
                      'directory {}'.format(state_dir))
    fname = glob.glob(os.path.join(state_dir, '*.nc'))[0]
    ds_first_run = xr.open_dataset(fname)

    # Loop over all rest runs and compare fluxes with the base run
    for i, n_proc in enumerate(list_n_proc):
        # Skip the first run
        if i == 0:
            continue
        # Read output states for this run
        state_dir = os.path.join(
            state_basedir,
            'processors_{}'.format(n_proc))
        if len(glob.glob(os.path.join(state_dir, '*.nc'))) > 1:
            warnings.warn('More than one netCDF file found under '
                          'directory {}'.format(state_dir))
        fname = glob.glob(os.path.join(state_dir, '*.nc'))[0]
        ds_current_run = xr.open_dataset(fname)
        # Compare current run with base run
        for var in ds_first_run.data_vars:
            npt.assert_array_equal(ds_current_run[var].values,
                                   ds_first_run[var].values,
                                   err_msg='States are not an exact match '
                                   'for variable: {}'.format(var))
