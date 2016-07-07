#!/usr/bin/env python

import os
from collections import OrderedDict

import traceback

import numpy as np
import pandas as pd
import glob
import re
import matplotlib.pyplot as plt
import warnings

from tonic.models.vic.vic import VICRuntimeError, default_vic_valgrind_error_code
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
    try:
        lines = string.decode().splitlines()
    except UnicodeDecodeError:
        lines = string.splitlines()
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
        print(fname)
        df = read_vic_ascii(fname)

        # check that each dataframe includes all timestamps
        if (start is not None) and (end is not None):
            check_completed(df, start, end)
        else:
            start = df.index[0]
            end = df.index[-1]
            print("start = %s" %str(start))
            print("end = %s" %str(end))


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
        time_cols = ['YEAR', 'MONTH', 'DAY', 'SEC']
        df.index = pd.to_datetime(df.YEAR * 10000 + df.MONTH * 100 +
                                    df.DAY, format='%Y%m%d')
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

def read_snotel_swe_obs(filename, science_test_data_dir, items):

    ''' reads in Snotel SWE obs and returns DataFrame
    '''

    filename_fullpath = os.path.join(science_test_data_dir,
                                    'datasets',
                                    items['archive'],
                                    'observations',
                                    filename)

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

    # remove year, month, day columns of DataFrame
    snotel_swe.drop(['YEAR', 'MONTH', 'DAY'], inplace=True, axis=1)

    snotel_swe.index = snotel_swe['DATES']

    return(snotel_swe)

def read_vic_42(lat, lng, science_test_data_dir, items):
    ''' reads output from VIC 4.2
    '''

    if 'ecflux' in items['compare_to']:
        vic_42_file = 'en_bal_%s_%s' % (lat, lng)
        vic_42_dir = os.path.join(science_test_data_dir, 'test_runs',
                                items['archive'], 'ecflux', 'results')

    elif 'snotel' in items['compare_to']:
        vic_42_file = 'outfile_%s_%s' % (lat, lng)
        vic_42_dir = os.path.join(science_test_data_dir, 'test_runs',
                                items['archive'], 'snotel', 'results')

    else:
        raise ValueError("this option has not yet been implemented")

    vic_42 = pd.read_csv(os.path.join(vic_42_dir, vic_42_file),
                            sep='\t',
                            skiprows=5)

    # remove comment sign from column names in DataFrame
    vic_42 = vic_42.rename(columns=lambda x: x.replace('#', ''))

    # remove spaces from column names in DataFrame
    vic_42 = vic_42.rename(columns=lambda x: x.replace(' ', ''))

    # add datetime object column to snotel DataFrame
    vic_42['DATES'] = pd.to_datetime(vic_42.YEAR * 10000 +
                                            vic_42.MONTH * 100 +
                                            vic_42.DAY,
                                            format='%Y%m%d')

    # remove year, day columns of DataFrame
    vic_42.drop(['YEAR', 'DAY'], inplace=True, axis=1)

    if ('HOUR' not in vic_42) and ('SEC' in vic_42):
        # add hour column for groupby
        vic_42['HOUR'] = (vic_42['SEC'] * (1/3600)
                                        ).astype(int)

    vic_42.index = vic_42.DATES

    return(vic_42)

def read_vic_5(lat, lng, science_test_data_dir, items):

    ''' read VIC 5.0.x outout
    '''

    if 'ecflux' in items['compare_to']:
        vic_5_file = 'en_bal_%s_%s.txt' % (lat, lng)
        vic_5_dir = os.path.join(science_test_data_dir, 'test_runs',
                                items['archive'], 'ecflux', 'results')

    elif 'snotel' in items['compare_to']:
        vic_5_file = 'outfile_%s_%s.txt' % (lat, lng)
        vic_5_dir = os.path.join(science_test_data_dir, 'test_runs',
                                items['archive'], 'snotel', 'results')

    else:
        raise ValueError("this option has not yet been implemented")

    vic_5 = pd.read_csv(os.path.join(vic_5_dir, vic_5_file),
                            skiprows=3,
                            sep='\t')

    # remove spaces from column names
    vic_5 = vic_5.rename(columns=lambda x: x.replace(' ', ''))

    # add datetime object column to snotel DataFrame
    vic_5['DATES'] = pd.to_datetime(vic_5.YEAR * 10000 +
                                            vic_5.MONTH * 100 +
                                            vic_5.DAY,
                                            format='%Y%m%d')

    # remove year, day columns of DataFrame
    vic_5.drop(['YEAR', 'DAY'], inplace=True, axis=1)

    if ('HOUR' not in vic_5) and ('SEC' in vic_5):
        # add hour column for groupby
        vic_5['HOUR'] = (vic_5['SEC'] * (1/3600)
                                        ).astype(int)

    vic_5.index = vic_5.DATES

    return(vic_5)

def plot_science_tests(driver, test_type, science_test_data_dir, result_dir,
                        plot_dir, plots_to_make, compare_data):

    ''' makes science test figures

    Parameters
    ----------
    driver: <str>
        Name of Driver
    test_type: <str>
        Name of test
    science_test_data_dir: <str>
        Science test data directory
    result_dir: <str>
        Result directory
    plot_dir: <str>
        Directory for output plots
    plots_to_make <Dict>
        Keys that indicate which plots should be made
    compare_data <Dict>
        Keys that indicate which datasets and model output to use for comparison

    Returns
    ----------
    '''
    if test_type == "science_test_snotel":
        plot_snotel_comparison(driver,
                                science_test_data_dir,
                                compare_data,
                                result_dir,
                                plot_dir,
                                plots_to_make)

    elif test_type == "science_test_fluxnet":
        plot_fluxnet_comparison(driver,
                                science_test_data_dir,
                                compare_data,
                                result_dir,
                                plot_dir,
                                plots_to_make)
    else:
        print("this has not been implemented in the VIC 5.0 science test suite")

def plot_snotel_comparison(driver, science_test_data_dir,
                            compare_data_dict,
                            result_dir, plot_dir,
                            plots_to_make):
    ''' makes snotel figures

    '''

    # plot settings
    plot_variables = ['OUT_SWE', 'OUT_ALBEDO', 'OUT_SALBEDO', 'OUT_SNOW_DEPTH',
                    'OUT_SNOW_CANOPY', 'OUT_SNOW_PACK_TEMP', 'OUT_SNOW_MELT',
                    'OUT_R_NET', 'OUT_LATENT', 'OUT_SENSIBLE']

    plot_units = ['mm', 'fraction', 'fraction', 'mm', '%', 'degrees C',
                    'mm', 'W/$m^2$', 'W/$m^2$', 'W/$m^2$']
    lw = 4.0

    for filename in os.listdir(os.path.join(science_test_data_dir,
                                            'datasets',
                                            'snotel',
                                            'observations')):

        # get lat/lng from filename
        file_split = re.split('_', filename)
        lng = file_split[3].split('.txt')[0]
        lat = file_split[2]

        # loop over data to compare
        data = {}
        for key, items in compare_data_dict.items():

            # read in data

            if key == "snotel":

                data[key] = read_snotel_swe_obs(filename,
                                                science_test_data_dir,
                                                items)

            elif key == "VIC.4.2.d":

                data[key] = read_vic_42(lat, lng,
                                        science_test_data_dir,
                                        items)

            else:

                data[key] = read_vic_5(lat, lng,
                                        science_test_data_dir,
                                        items)

        # loop over variables to plot
        for i, plot_variable in enumerate(plot_variables):

            if 'water_year' in plots_to_make:

                fig, ax = plt.subplots(figsize=(10,10))

                df = pd.DataFrame({key: d[plot_variable] for key, d in
                                            data.items() if plot_variable in d})

                # merge DateTime indexes (time series lengths are different)
                if 'snotel' in data.keys():
                    df_join = data.pop('snotel').merge(next(iter(data.values())),
                                                    how='outer', on='DATES')
                    df.index = df_join['DATES']
                else:
                    df.index = next(iter(data.values()))['DATES']

                for key, series in df.iteritems():
                    series.plot(linewidth=lw, ax=ax, color=compare_data_dict[
                                                            key]['color'])

                ax.legend(loc='upper left')
                ax.set_ylabel("%s [%s]" % (plot_variable, plot_units[i]))

                # save figure
                os.makedirs(os.path.join(plot_dir, plot_variable), exist_ok=True)
                plotname = '%s_%s.png' % (lat, lng)
                savepath = os.path.join(plot_dir, plot_variable, plotname)
                plt.savefig(savepath, bbox_inches='tight')

                plt.clf()
                plt.close()

def check_site_files(obs_dir, subdir):
    if os.listdir(os.path.join(obs_dir, subdir)) != []:
        return(True)
    else:
        return(False)

def get_fluxnet_lat_lon(obs_dir, subdir):

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

    return(lat, lng)

def read_fluxnet_obs(subdir, science_test_data_dir, items):

    filename = '%s.stdfmt.hourly.local.txt' %subdir

    # column names for DataFrame (same as VIC variable names)
    names = ['YEAR', 'MONTH', 'DAY', 'HOUR', 'PREC', 'AIR_TEMP', 'SWDOWN',
            'LWDOWN','OUT_REL_HUMID', 'PRESSURE', 'WIND', 'OUT_EVAP',
            'SOIL_TEMP_DEPTH1', 'SOIL_TEMP_DEPTH2', 'SOIL_TEMP_DEPTH3',
            'SOIL_TEMP_DEPTH4', 'SOIL_TEMP_DEPTH5', 'OUT_SOIL_MOIST1',
            'OUT_SOIL_MOIST2', 'OUT_SOIL_MOIST3', 'OUT_SOIL_MOIST4',
            'OUT_SOIL_MOIST5', 'OUT_SOIL_TEMP1', 'OUT_SOIL_TEMP2',
            'OUT_SOIL_TEMP3', 'OUT_SOIL_TEMP4', 'OUT_SOIL_TEMP5', 'SWNET',
            'LWNET', 'OUT_SENSIBLE', 'OUT_LATENT', 'OUT_GRND_FLUX']

    # read in data with -9999.0000 as NaNs
    obs_dir = os.path.join(science_test_data_dir, 'datasets',
                            'ec_flux_towers', 'obs')
    ecflux_df = pd.read_csv(os.path.join(obs_dir, subdir, filename),
            skiprows=0,
            delim_whitespace=True,
            header=None,
            names=names,
            na_values=-9999.0000)

    ecflux_df['DATES'] = pd.to_datetime(ecflux_df.YEAR * 10000 +
                                 ecflux_df.MONTH * 100 +
                                 ecflux_df.DAY,format='%Y%m%d')

    ecflux_df.index = ecflux_df.DATES

    ecflux_df.drop(['YEAR', 'DAY'], inplace=True, axis=1)

    if ('HOUR' not in ecflux_df) and ('SEC' in ecflux_df):
        # add hour column for groupby
        ecflux_df['HOUR'] = (ecflux_df['SEC'] * (1/3600)
                                        ).astype(int)

    return(ecflux_df)

def plot_fluxnet_comparison(driver, science_test_data_dir,
                            compare_data_dict,
                            result_dir, plot_dir,
                            plots_to_make):

    ''' makes Ameriflux figures
    '''

    # loop over Ameriflux sites
    obs_dir = os.path.join(science_test_data_dir,
                            'datasets',
                            'ec_flux_towers',
                            'obs')

    for subdir in os.listdir(obs_dir):

        if check_site_files(obs_dir, subdir):
            # get CSV file from site directory to get lat/lng for site
            lat, lng = get_fluxnet_lat_lon(obs_dir, subdir)

            # loop over data to compare
            data = {}
            for key, items in compare_data_dict.items():

                if key == "ecflux":
                    try:
                        # load Ameriflux data
                        data[key] = read_fluxnet_obs(subdir,
                                                science_test_data_dir, items)
                    except OSError:
                        warnings.warn("this site does not have data")

                elif key == "VIC.4.2.d":
                    try:
                        # load VIC 4.2 simulations
                        data[key] = read_vic_42(lat, lng, science_test_data_dir,
                                                items)

                    except OSError:
                        warnings.warn("this site has a lat/lng precision issue")

                    #except AssertionError:
                        #print("blah")
                else:
                    try:
                        # load VIC 5 simulations
                        data[key] = read_vic_5(lat, lng, science_test_data_dir,
                                                items)
                    except OSError:
                        warnings.warn("this site has a lat/lng precision issue")

            # make figures
            vic_vars = ['OUT_LATENT', 'OUT_SENSIBLE']
            variable_names = ['Latent Heat', 'Sensible Heat']

            # plot preferences
            lw = 4.0
            fs = 15
            alpha = 1.0
            dpi = 150

            if 'annual_mean_diurnal_cycle' in plots_to_make:

                # make annual mean diurnal cycle plots
                f, axarr = plt.subplots(2, 1, figsize=(8,8), sharex=True)

                for i, vic_var in enumerate(vic_vars):

                    # calculate annual mean diurnal cycle for each DataFrame
                    annual_mean = {}
                    for (key, df) in data.items():
                        annual_mean[key] = pd.DataFrame(df[vic_var].groupby(
                                                        df['HOUR']).mean())

                    df = pd.DataFrame({key: d[vic_var] for key, d in
                                    annual_mean.items() if vic_var in d})

                    for key, series in df.iteritems():
                        series.plot(linewidth=lw,
                                    ax=axarr[i],
                                    color=compare_data_dict[key]['color'])

                    axarr[i].legend(loc='upper left')
                    axarr[i].set_ylabel('%s W / $m^2$' %variable_names[i],
                                        size=fs)
                    axarr[i].set_xlim([0,24])

                # save plot
                plotname = '%s_%s.png' % (lat, lng)
                os.makedirs(os.path.join(plot_dir, 'annual_mean'),
                                        exist_ok=True)
                savepath = os.path.join(plot_dir, 'annual_mean', plotname)
                plt.savefig(savepath, bbox_inches='tight', dpi=dpi)

                plt.clf()
                plt.close()

            if 'monthly_mean_diurnal_cycle' in plots_to_make:

                # make monthly mean diurnal cycle plots
                f, axarr = plt.subplots(2, 12, figsize=(35,7), sharex=True,
                                        sharey=True)

                months = ['January', 'February', 'March', 'April', 'May',
                            'June', 'July', 'August', 'September',
                            'October', 'November', 'December']

                for i, vic_var in enumerate(vic_vars):
                    for j in range(12):

                        # calculate monthly mean diurnal cycle
                        monthly_mean = {}
                        for (key, df) in data.items():
                            monthly_mean[key] = pd.DataFrame(df[vic_var].
                                                    groupby([df['MONTH'],
                                                    df['HOUR']]).mean())

                        df = pd.DataFrame({key: d[vic_var] for key, d in
                                        monthly_mean.items() if vic_var in d})

                        for key, series in df.iteritems():
                            series[j+1].plot(linewidth=lw,
                                        ax=axarr[i,j],
                                        color=compare_data_dict[key]
                                        ['color'],
                                        alpha=alpha)

                        if i == 0 and j == 11:
                            axarr[i,j].legend(loc='center left',
                                                bbox_to_anchor=(1, 0.5))
                        axarr[i,j].set_ylabel('%s W/$m^2$' %variable_names[i],
                                            size=fs)
                        axarr[i,j].set_xlim([0,24])
                        if i == 0:
                            axarr[i,j].set_title(months[j])

                # save plot
                plotname = '%s_%s.png' % (lat, lng)
                os.makedirs(os.path.join(plot_dir, 'monthly_mean'),
                            exist_ok=True)
                savepath = os.path.join(plot_dir,
                                        'monthly_mean', plotname)
                plt.savefig(savepath, bbox_inches='tight', dpi=dpi)

                plt.clf()
                plt.close()

        else:
            subdir_files = False
            warnings.warn("this site does not have data")
