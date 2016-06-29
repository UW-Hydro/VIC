#!/usr/bin/env python

import os
from collections import OrderedDict

import traceback

import numpy as np
import pandas as pd
import glob
import re

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
        raise error
    error_message = error
    traceback.print_stack()
    print('\t{0}'.format(test_comment))
    print('\t{0}'.format(error_message))
    if tail is not None:
        print('\tLast {0} lines of standard out:'.format(ERROR_TAIL))
        print_tail(tail, n=ERROR_TAIL)

    return test_comment, error_message


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


def check_multistream(fnames, driver):
    '''

    '''

    if driver.lower() != 'classic':
        raise ValueError('only classic driver is supported in this test')

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

def tsplit(string, delimiters):
    '''Behaves like str.split but supports multiple delimiters.

    '''

    delimiters = tuple(delimiters)
    stack = [string,]

    for delimiter in delimiters:
        for i, substring in enumerate(stack):
            substack = substring.split(delimiter)
            stack.pop(i)
            for j, _substring in enumerate(substack):
                stack.insert(i+j, _substring)

    return stack

def plot_science_tests(driver, testname, result_dir, plot_dir, vic_42_dir, vic_50_dir,
                        obs_dir, plots_to_make):

    ''' makes science test figures

    Parameters
    ----------
    driver: <str>
        Name of Driver
    testname: <str>
        Name of test
    result_dir: <str>
        Result directory
    plot_dir: <str>
        Directory for output plots
    vic_42_dir: <str>
        Directory for VIC 4.2 archive
    vic_50_dir <str>
        Directory for VIC 5.0 archive
    obs_dir <str>
        Directory for observations archive
    plots_to_make <Dict>
        Keys that indicate which plots should be made

    Returns
    ----------
    '''
    if testname == "science_test_snotel":
        plot_snotel_comparison(driver, testname, result_dir, plot_dir, vic_42_dir, vic_50_dir,
                                obs_dir, plots_to_make)
    elif testname == "science_test_fluxnet"
        plot_fluxnet_comparison(driver, testname, result_dir, plot_dir, vic_42_dir, vic_50_dir,
                                obs_dir, plots_to_make)
    else:
        print("this has not yet been implemented in the VIC 5.0 science test suite")

def plot_snotel_comparison(driver, testname, result_dir, plot_dir, vic_42_dir, vic_50_dir,
                        obs_dir, plots_to_make):
    ''' makes snotel figures

    '''
    for filename in os.listdir(obs_dir):
        filename_fullpath = os.path.join(obs_dir, filename)
        file_split = re.split('_', filename)
        lng = file_split[3].split('.txt')[0]
        lat = file_split[2]

        # load snotel obs
        snotel_swe = pd.read_csv(filename_fullpath,
                                skiprows=0,
                                delim_whitespace=True,
                                names=['YEAR', 'MONTH', 'DAY', 'OUT_SWE'])
        # add datetime object column to snotel DataFrame
        snotel_swe['DATES'] = pd.to_datetime(snotel_swe.YEAR * 10000 +
                                            snotel_swe.MONTH * 100 +
                                            snotel_swe.DAY,
                                            format='%Y%m%d')

        # load VIC 4.2 data
        vic_42_file = '%s_%s' % (lat, lng)

        vic_42 = pd.read_csv(os.path.join(vic_42_dir, vic_42_file),
                                sep='\t',
                                skiprows=5)

        # remove comment sign from column names in DataFrames
        vic_42 = vic_42.rename(columns=lambda x: x.replace('#', ''))

        # load VIC 5.0 data

        vic_50_file = '%s_%s.txt' %(lat, lng)

        vic_50 = pd.read_csv(os.path.join(vic_50_dir, vic_50_file),
                                skiprows=3,
                                sep='\t')

        vic_5x_file = '%s_%s.txt' %(lat, lng)

        vic_50x = pd.read_csv(os.path.join(result_dir, vic_5x_file),
                                skiprows=3,
                                sep='\t')

        # variables to plot
        plot_variables = ['OUT_SWE', 'OUT_ALBEDO', 'OUT_SALBEDO', 'OUT_SNOW_DEPTH',
                        'OUT_SNOW_CANOPY', 'OUT_SNOW_PACK_TEMP', 'OUT_SNOW_MELT',
                        'OUT_R_NET', 'OUT_LATENT', 'OUT_SENSIBLE']

        plot_units = ['mm', 'fraction', 'fraction', 'mm', '%', 'degrees C',
                        'mm', 'W/$m^2$', 'W/$m^2$', 'W/$m^2$']

        len_dates = len(snotel_swe['DATES'])

        # plotting preferences
        lw = 4.0
        # loop over variables to plot

        for i, plot_variable in enumerate(plot_variables):

            if 'water_year' in plots_to_make:

                plt.figure(figsize=(10,10))

                if plot_variable == "OUT_SWE":

                    # plot SnoTel SWE observations
                    plt.plot(snotel_swe['DATES'], snotel_swe[plot_variable],
                            'k', label='Snotel', linewidth=lw)

                # plot VIC 4.2 simulations
                plt.plot(snotel_swe['DATES'], vic_42_snow[plot_variable][:len_dates],
                        'b', label='VIC 4.2', linewidth=lw)

                # plot VIC 5.0 simulations
                plt.plot(snotel_swe['DATES'], vic_50_snow[plot_variable][:len_dates],
                        'r', label='VIC 5.0', linewidth=lw)

                # plot VIC 5.0.x simulations
                plt.plot(snotel_swe['DATES'], vic_50x_snow[snow_variable][:len_dates]),
                        'y', label='VIC 5.0.x', linewidth=lw)

                plt.title(plot_variable)
                plt.legend(loc='upper left')
                plt.ylabel(plot_units[i])

                # save figure
                os.makedirs(os.path.join(plot_dir, 'plot_variable'), exist_ok=True)
                plotname = '%s_%s.png'
                savepath = os.path.join(plot_dir, 'plot_variable', plotname)
                plt.savefig(savepath)


def plot_fluxnet_comparison(driver, testname, result_dir, plot_dir, vic_42_dir, vic_50_dir,
                        obs_dir, plots_to_make):

    ''' makes Ameriflux figures

    '''

    # loop over Ameriflux sites
    for subdir in os.listdir(obs_dir):

        lat_issue = False

        if os.listdir(os.path.join(obs_dir, subdir) != []:

            subdir_files = True

            # get CSV file from site directory to get lat/lng for site
            try:
                site_csv_file = glob.glob(os.path.join(obs_dir, subdir, 'AMF*.csv'))[0]
            except IndexError:
                site_csv_file = glob.glob(os.path.join(obs_dir, subdir, 'us*.csv'))[0]
            with open(site_csv_file) as file:
                second_line = list(file)[1]

            # parse line from header to get lat/lng
            str_split = tsplit(second_line, ('Latitude: ', 'Longitude: ', 'Elevation (masl): '))
            lat = str_split[1].strip()
            lng = str_split[2].strip()

            # load Ameriflux data
            filename = '%s.stdfmt.hourly.local.txt' %subdir

            # column names for DataFrame
            # column names
            names = ['Year', 'Month', 'Day', 'Hour', 'P', 'Tair', 'SWdown', 'LWdown', 'RH', 'Patm',
                    'Wind', 'ET', 'Tsdep1', 'Tsdep2', 'Tsdep3', 'Tsdep4', 'Tsdep5', 'SM1', 'SM2',
                    'SM3', 'SM4', 'SM5', 'Tsoil1', 'Tsoil2', 'Tsoil3', 'Tsoil4', 'Tsoil5', 'SWnet',
                    'LWnet', 'H', 'LE', 'FG']

            # read in data with -9999.0000 as NaNs
            ecflux_df = pd.read_csv(os.path.join(obs_dir, subdir, filename),
                    skiprows=0,
                    delim_whitespace=True,
                    header=None,
                    names=names,
                    na_values=-9999.0000)

            ecflux_df['DATES'] = pd.to_datetime(ecflux_df.Year * 10000 +
                                         ecflux_df.Month * 100 +
                                         ecflux_df.Day,format='%Y%m%d')

        else:
            subdir_files = False

        # load VIC 4.2 simulations
        ebal_42_file = 'en_bal_%s_%s' % (lat, lng)

        try:

            ebal_42 = pd.read_csv(os.path.join(vic_42_dir,
                                            ebal_42_file),
                                            sep='\t',
                                            skiprows=5)

            # remove comment sign from column names in DataFrame
            ebal_42 = ebal_42.rename(columns=lambda x: x.replace('#', ''))
            ebal_42 = ebal_42.rename(columns=lambda x: x.replace(' ', ''))

            ebal_42['DATES'] = pd.to_datetime(ebal_42.YEAR * 10000 +
                                     ebal_42.MONTH * 100 +
                                     ebal_42.DAY, format='%Y%m%d')
        except OSError:
            lat_issue=True

        # load VIC 5.0 simulations
        ebal_50_file = 'en_bal_%s_%s.txt' %(lat, lng)

        try:

            # load VIC 5 simulations
            ebal_50 = pd.read_csv(os.path.join(vic_50_dir,
                                 ebal_50_file),
                                 skiprows=3,
                                 sep='\t')

            # remove space from column names in DataFrame
            ebal_50 = ebal_50.rename(columns=lambda x: x.replace(' ', ''))

            ebal_50['DATES'] = pd.to_datetime(ebal_50.YEAR * 10000 +
                                         ebal_50.MONTH * 100 +
                                         ebal_50.DAY, format='%Y%m%d')

        except OSError:
            # To-Do: deal with lat/lng precision issue

        try:

            # load VIC 5.x simulations
            ebal_50x = pd.read_csv(os.path.join(result_dir,
                                 ebal_50_file),
                                 skiprows=3,
                                 sep='\t')
            # remove space from column names in DataFrame
            ebal_50x = ebal_50x.rename(columns=lambda x: x.replace(' ', ''))

            ebal_50x['DATES'] = pd.to_datetime(ebal_50x.YEAR * 10000 +
                                             ebal_50x.MONTH * 100 +
                                             ebal_50x.DAY, format='%Y%m%d')

        except:
            OSError
            # To-Do: deal with lat/lng precision issue for some sites

        # make figures

        ecflux_vars = ['LE', 'H']
        vic_vars = ['OUT_LATENT', 'OUT_SENSIBLE']
        variable_names = ['Latent Heat', 'Sensible Heat']

        # plot preferences
        lw = 4.0
        fs = 15

        if not lat_issue:
            if 'annual_mean_diurnal_cycle' in plots_to_make:

                # make annual mean diurnal cycle plots
                f, axarr = plt.subplots(2, 1, figsize=(8,8))

                for i, ecflux_var in enumerate(ecflux_vars):

                    ebal_42.index = ebal_42.DATES

                    # plt VIC 4.2
                    axarr[i].plot(ebal_42[vic_vars[i]].groupby(ebal_42['HOUR']).mean(),
                                'b', label='VIC 4.2', linewidth=lw)

                    # plot VIC 5.0
                    ebal_50.index = ebal_50.DATES

                    # convert seconds column to hours for groupby
                    ebal_50['HOUR'] = (ebal_50['SEC'] * (1/3600)).astype(int)

                    axarr[i].plot(ebal_50[vic_vars[i]].groupby(ebal_50['HOUR']).mean(),
                                'r', label='VIC 5.0', linewidth=lw)

                    # plot VIC 5.0.x
                    ebal_50x.index = ebal50x.DATES

                    # convert seconds column to hours for groupby
                    ebal_50x['HOUR'] = (ebal_50x['SEC'] * (1/3600)).astype(int)

                    axarr[i].plot(ebal_50x[vic_vars[i]].groupby(ebal_50x['HOUR']).mean(),
                                'y', label='VIC 5.0.x', linewidth=lw)

                    axarr[i].legend(loc='upper left')
                    axarr[i].set_title(variable_names[i])
                    axarr[i].set_ylabel('W / $m^2$', size=fs))
                    axarr[i].set_xlim([0,24])

                # save plot
                plotname = '%s_%s.png'
                os.makedirs(os.path.join(plot_dir, 'annual_mean'), exist_ok=True)
                savepath = os.path.join(plot_dir, 'annual_mean', plotname)
                plt.savefig(savepath, bbox_inches='tight')

            elif 'monthly_mean_diurnal_cycle' in plots_to_make:

                for i, ecflux_var in enumerate(ecflux_vars):
                    for j in range(12):

                        # plot VIC 4.2
                        # make datetime column the DataFrame index
                        ebal_42.index = ebal_42.DATES

                        # calculate monthly mean diurnal cycle
                        vic_monthly_mean_42 = ebal_42[vic_vars[i]].groupby(
                                                    [ebal_42['MONTH'], ebal_42['HOUR']]).mean()

                        axarr[i,j].plot(vic_monthly_mean_42[j+1], 'b',
                                        label='VIC 4.2', linewidth=lw)

                        # plot VIC 5.0

                        # make datetime column the DataFrame index
                        ebal_50.index = ebal_50.DATES

                        # convert seconds column to hours for groupby
                        ebal_50['HOUR'] = (ebal_50['SEC'] * (1/3600)).astype(int)

                        # calculate monthly mean diurnal cycle
                        vic_monthly_mean_50 = ebal_50[vic_vars[i]].groupby(
                                                    [ebal_50['MONTH'], ebal_50['HOUR']]).mean()

                        axarr[i,j].plot(vic_monthly_mean_50[j+1], 'r',
                                        label='VIC 5.0', linewidth=lw)

                        # plot VIC 5.0.x

                        # make datetime column the DataFrame index
                        ebal_50x.index = ebal_50x.DATES

                        # convert seconds column to hours for groupby
                        ebal_50x['HOUR'] = (ebal_50x['SEC'] * (1/3600)).astype(int)

                        # calculate monthly mean diurnal cycle
                        vic_monthly_mean_50x = ebal_50x[vic_vars[i]].groupby(
                                                    [ebal_50x['MONTH'], ebal_50x['HOUR']]).mean()

                        axarr[i,j].plot(vic_monthly_mean_50x[j+1], 'y',
                                        label='VIC 5.0.x', linewidth=lw)

                # save plot
                plotname = '%s_%s.png'
                os.makedirs(os.path.join(plot_dir, 'monthly_mean'), exist_ok=True)
                savepath = os.path.join(plot_dir, 'montly_mean', plotname)
                plt.savefig(savepath, bbox_inches='tight')

            else:
                print("this has not yet been implemented")
