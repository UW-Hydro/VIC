#!/usr/bin/env python
''' VIC exact restart testing '''

import os
from collections import OrderedDict

import traceback

import numpy as np
import pandas as pd

from tonic.models.vic.vic import (VICRuntimeError,
                                  default_vic_valgrind_error_code)
from tonic.testing import check_completed, check_for_nans, VICTestError

OUTPUT_WIDTH = 100
ERROR_TAIL = 20  # lines


class VICReturnCodeError(Exception):
    pass


class VICValgrindError(Exception):
    pass


def setup_test_dirs(testname, out_dir, mkdirs=['results', 'state',
                                               'logs', 'plots']):
    '''create test directories for testname'''

    dirs = OrderedDict()
    dirs['test'] = os.path.join(out_dir, testname)
    for d in mkdirs:
        dirs[d] = os.path.join(dirs['test'], d)

    for dirname in dirs.values():
        os.makedirs(dirname, exist_ok=True)

    return dirs


def print_test_dict(d):
    '''print a nicely formatted set of test results'''
    print('{0: <48} | {1: <6} | {2}'.format('Test Name', 'Passed', 'Comment'))
    print('-'.ljust(OUTPUT_WIDTH, '-'))
    for k, v in d.items():
        print('{0: <48} | {1: <6} | {2}'.format(clip_string(v.name, 48),
                                                str(bool(v.passed)),
                                                v.comment))
        print('-'.ljust(OUTPUT_WIDTH, '-'))


def clip_string(string, length=50):
    if len(string) > length:
        string = string[:length - 3] + '...'
    return string


def print_tail(string, n=20, indent='\t--->'):
    '''print tail of multiline string'''
    lines = string.decode().splitlines()
    for l in lines[-n:]:
        print('{0}{1}'.format(indent, l))


def replace_global_values(gp, replace):
    '''given a multiline string that represents a VIC global parameter file,
       loop through the string, replacing values with those found in the
       replace dictionary'''
    gpl = []
    for line in iter(gp.splitlines()):
        line_list = line.split()
        if line_list:
            key = line_list[0]
            if key in replace:
                value = replace.pop(key)
                val = list([str(value)])
            else:
                val = line_list[1:]
            gpl.append('{0: <20} {1}\n'.format(key, ' '.join(val)))

    if replace:
        for key, val in replace.items():
            try:
                value = ' '.join(val)
            except:
                value = val
            gpl.append('{0: <20} {1}\n'.format(key, value))

    return gpl


def drop_tests(config, driver):
    '''helper function to remove tests that should not be run for driver'''
    new = {}
    for key, test_cfg in config.items():
        try:
            if test_cfg['driver'].lower() == driver.lower():
                new[key] = test_cfg
        except KeyError:
            raise KeyError('test configuration must specify driver')
    return new


def pop_run_kwargs(config):
    '''pop run kwargs for VIC executable'''
    run_kwargs = {}
    run_kwargs['valgrind'] = config.pop('valgrind', False)
    run_kwargs['mpi_proc'] = config.pop('mpi_proc', None)
    return run_kwargs


def check_returncode(exe, expected=0):
    '''check return code given by VIC, raise error if appropriate'''
    if exe.returncode == expected:
        return None
    elif exe.returncode == default_vic_valgrind_error_code:
        raise VICValgrindError(
            'Valgrind raised an error when running: "{}"'.format(exe.argstring))
    else:
        raise VICReturnCodeError(
            'VIC return code ({0}) did not match expected ({1}) when running '
            '"{2}"'.format(exe.returncode, expected, exe.argstring))


def process_error(error, vic_exe):
    '''Helper function to process possible error raised during testing'''
    tail = None
    if isinstance(error, VICRuntimeError):
        test_comment = 'Test failed during simulation'
        tail = vic_exe.stderr
    elif isinstance(error, VICTestError):
        test_comment = 'Test failed during testing of output files'
    elif isinstance(error, VICValgrindError):
        test_comment = 'Test failed due to memory error detected by valgrind'
        tail = vic_exe.stderr
    elif isinstance(error, VICReturnCodeError):
        test_comment = 'Test failed due to incorrect return code'
        tail = vic_exe.stderr
    elif isinstance(error, AssertionError):
        test_comment = 'AssertionError raised during testing'
    else:
        test_comment = 'Unknown test failure'
        traceback.print_stack()

    print('\t{0}'.format(test_comment))
    print('\t{0}'.format(error))
    if tail is not None:
        print('\tLast {0} lines of standard out:'.format(ERROR_TAIL))
        print_tail(tail, n=ERROR_TAIL)

    return test_comment, error


def test_classic_driver_all_complete(fnames):
    '''
    Test that all VIC files in fnames have the same first and last index
    position
    '''
    start = None
    end = None
    for fname in fnames:
        df = read_vic_ascii(fname)

        # check that each dataframe includes all timestamps
        if (start is not None) and (end is not None):
            check_completed(df, start, end)
        else:
            start = df.index[0]
            end = df.index[-1]


def test_classic_driver_no_output_file_nans(fnames):
    '''Test that all VIC classic driver output files in fnames have no nans'''
    for fname in fnames:
        df = read_vic_ascii(fname)
        check_for_nans(df)


# TODO: Update tonic version of this function, need to check that subdaily works
def read_vic_ascii(filepath, parse_dates=True, datetime_index=None, sep='\t',
                   comment='#', **kwargs):
    '''Generic reader function for VIC ASCII output with a standard header
    filepath: path to VIC output file
    header (True or False):  Standard VIC header is present
    parse_dates (True or False): Parse dates from file
    datetime_index (Pandas.tseries.index.DatetimeIndex):  Index to use as
    datetime index names (list like): variable names
    **kwargs: passed to Pandas.read_table
    returns Pandas.DataFrame
    '''

    df = pd.read_table(filepath, sep=sep, comment=comment, **kwargs)
    # Strip extra whitespace around variable names
    df.rename(columns=lambda x: x.strip(), inplace=True)

    if parse_dates and datetime_index:
        raise ValueError('cannot specify both parse_dates and datetime_index')

    if parse_dates:
        time_cols = ['YEAR', 'MONTH', 'DAY']
        df.index = pd.to_datetime(df[time_cols])
        if 'SEC' in df:
            df.index += pd.Series([pd.Timedelta(s, unit='s') for s in df['SEC']],
                                  index=df.index)
            time_cols.append('SEC')
        df.drop(time_cols, axis=1)

    if datetime_index is not None:
        df.index = datetime_index

    return df


def find_global_param_value(gp, param_name):
    ''' Return the value of a global parameter

    Parameters
    ----------
    gp: <str>
        Global parameter file, read in by read()
    param_name: <str>
        The name of the global parameter to find

    Returns
    ----------
    line_list[1]: <str>
        The value of the global parameter
    '''
    for line in iter(gp.splitlines()):
        line_list = line.split()
        if line_list == []:
            continue
        key = line_list[0]
        if key == param_name:
            return line_list[1]


def check_multistream_classic(fnames):
    '''
    Test the multistream aggregation in the classic driver
    '''

    how_dict = {'OUT_ALBEDO': 'max',
                'OUT_SOIL_TEMP_1': 'min',
                'OUT_PRESSURE': 'sum',
                'OUT_AIR_TEMP': 'first',
                'OUT_SWDOWN': 'mean',
                'OUT_LWDOWN': 'last'}

    streams = {}  # Dictionary to store parsed stream names
    gridcells = []  # list of gridcells (lat_lon)

    for path in fnames:
        # split up the path name to get info about the stream
        resultdir, fname = os.path.split(path)
        pieces = os.path.splitext(fname)[0].split('_')
        gridcells.append('_'.join(pieces[-2:]))
        stream = '_'.join(pieces[:-2])  # stream name
        freq_n = pieces[-3]  # stream frequency n

        # set the stream frequency for pandas resample
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

    # unique gridcells
    gridcells = list(set(gridcells))

    # Loop over all grid cells in result dir
    for gridcell in gridcells:
        fname = os.path.join(resultdir, '{}_{}.txt'.format(inst_stream, gridcell))
        instant_df = read_vic_ascii(fname)

        # Loop over all streams
        for stream, freq in streams.items():
            fname = os.path.join(resultdir, '{}_{}.txt'.format(stream, gridcell))
            agg_df = read_vic_ascii(fname)

            # Setup the resample of the instantaneous data
            rs = instant_df.resample(freq)

            # Loop over the variables in the stream
            for key, how in how_dict.items():
                # Get the aggregated values (from VIC)
                actual = agg_df[key].values
                # Calculated the expected values based on the resampling from pandas
                expected = rs[key].aggregate(how).values

                # Compare the actual and expected (with tolerance)
                np.testing.assert_almost_equal(
                    actual, expected, decimal=4,
                    err_msg='Variable=%s, freq=%s, how=%s: '
                            'failed comparison' % (key, freq, how))
