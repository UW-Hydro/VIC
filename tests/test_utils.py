# Builtin libs
import os
import re
import glob
import traceback
import warnings
from collections import OrderedDict, namedtuple
import multiprocessing as mp

# Computation libs
import numpy as np
import pandas as pd
import xarray as xr

# Plotting libs
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

# Tools from tonic
from tonic.models.vic.vic import (VICRuntimeError,
                                  default_vic_valgrind_error_code)
from tonic.testing import check_completed, check_for_nans, VICTestError

OUTPUT_WIDTH = 100
ERROR_TAIL = 20  # lines

VICOutFile = namedtuple('vic_out_file',
                        ('dirpath', 'prefix', 'lat', 'lon', 'suffix'))


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

    if not isinstance(driver, list):  # if single driver
        for key, test_cfg in config.items():
            try:
                if not isinstance(test_cfg['driver'], list):
                    if test_cfg['driver'].lower() == driver.lower():
                        new[key] = test_cfg
            except KeyError:
                raise KeyError('test configuration must specify driver')
    else:  # if multiple drivers
        for key, test_cfg in config.items():
            try:
                if isinstance(test_cfg['driver'], list):
                    # check whether the test has the same number of drivers
                    if len(test_cfg['driver']) == len(driver):
                        # check whether the test wants to test the same drivers
                        flag = 1
                        for d in driver:
                            if d not in test_cfg['driver']:
                                flag = 0
                        if flag == 1:
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
            'Valgrind raised an error when running: \
            "{}"'.format(exe.argstring))
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
        print(fname)
        df = read_vic_ascii(fname)
        check_for_nans(df)


# TODO: Update tonic version of this function,
# need to check that subdaily works
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
        # add datetime index
        time_cols = ['YEAR', 'MONTH', 'DAY']
        df.index = pd.to_datetime(df[time_cols])
        if 'SEC' in df:
            df.index += pd.Series(
                [pd.Timedelta(s, unit='s') for s in df['SEC']], index=df.index)
            time_cols.append('SEC')
        df.drop(time_cols, inplace=True, axis=1)

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
    Test the multistream aggregation in the classic driver '''

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
        fname = os.path.join(resultdir,
                             '{}_{}.txt'.format(inst_stream, gridcell))
        instant_df = read_vic_ascii(fname)

        # Loop over all streams
        for stream, freq in streams.items():
            fname = os.path.join(resultdir,
                                 '{}_{}.txt'.format(stream, gridcell))
            agg_df = read_vic_ascii(fname)

            # Setup the resample of the instantaneous data
            rs = instant_df.resample(freq)

            # Loop over the variables in the stream
            for key, how in how_dict.items():
                # Get the aggregated values (from VIC)
                actual = agg_df[key].values
                # Calculated the expected values based on the resampling from
                # pandas
                expected = rs[key].aggregate(how).values

                # Compare the actual and expected (with tolerance)
                np.testing.assert_almost_equal(
                    actual, expected, decimal=4,
                    err_msg='Variable=%s, freq=%s, how=%s: '
                            'failed comparison' % (key, freq, how))


def setup_subdirs_and_fill_in_global_param_driver_match_test(
        dict_s, result_basedir, state_basedir, test_data_dir):
    ''' Fill in global parameter output directories for multiple driver runs
        for driver-match testing

    Parameters
    ----------
    dict_s: <dict of string.Template>
        A dict of template of the global param file to be filled in
        Keys: driver name
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
    dict_global_param: <dict>
        A dict of global parameter strings to be run with parameters filled in

    Require
    ----------
    os
    '''

    dict_global_param = {}
    for driver in dict_s.keys():
        # Set up subdirectories for results and states
        result_dir = os.path.join(result_basedir, driver)
        state_dir = os.path.join(state_basedir, driver)
        os.makedirs(result_dir, exist_ok=True)
        os.makedirs(state_dir, exist_ok=True)

        # Fill in global parameter options
        s = dict_s[driver]
        dict_global_param[driver] = s.safe_substitute(
            test_data_dir=test_data_dir,
            result_dir=result_dir,
            state_dir=state_dir)

    return(dict_global_param)


def parse_classic_driver_outfile_name(fname):
    '''helper function to parse VIC classic driver output file name'''
    resultdir, filename = os.path.split(fname)
    prefix, suffix = os.path.splitext(filename)
    pieces = prefix.split('_')
    lat, lon = map(float, pieces[-2:])
    return VICOutFile(resultdir, prefix, lat, lon, suffix)


def check_drivers_match_fluxes(list_drivers, result_basedir):
    ''' Check whether the flux results are similar cross multiple drivers

    Parameters
    ----------
    list_drivers: <list>
        A list of driver names to be compared
        e.g., ['classic'; 'image']
        NOTE: must have classic driver; classic driver will be the base for
        comparison
    result_basedir: <str>
        Base directory of output fluxes results; results for drivers are
        subdirectories under the base directory

    Require
    ----------
    glob
    xarray
    numpy
    warnings
    collections.namedtuple
    parse_classic_driver_outfile_name
    VICOutFile
    read_vic_ascii
    '''

    # Identify all classic driver output flux files
    try:
        list_fnames_classic = glob.glob(
            os.path.join(result_basedir, 'classic', '*'))
    except:
        raise ValueError('incorrect classic driver output for driver-match '
                         'test')

    # Loop over all other drivers and compare with classic driver
    for driver in list_drivers:
        # skip classic driver
        if driver == 'classic':
            continue
        # if image driver
        if driver == 'image':
            # load flux file
            if len(glob.glob(os.path.join(
                    result_basedir, driver, '*.nc'))) > 1:
                warnings.warn('More than one netCDF file found under'
                              'directory {}'.format(result_basedir))
            fname = glob.glob(os.path.join(result_basedir, driver, '*.nc'))[0]
            ds_image = xr.open_dataset(fname)

            # loop over each grid cell from classic driver
            for fname in list_fnames_classic:
                gcell = parse_classic_driver_outfile_name(fname)
                df_classic = read_vic_ascii(fname)
                ds_image_cell = ds_image.sel(lat=gcell.lat, lon=gcell.lon,
                                             method='nearest')
                # compare each variable
                for var in ds_image_cell.data_vars:
                    # if one [time] dimension
                    if len(ds_image_cell[var].coords) == 3:
                        # determine precision for comparison
                        # --- if all zeros for this variable, set
                        # --- decimal = 2 --- #
                        if np.sum(np.absolute(ds_image_cell[var].values)) == 0:
                            decimal = 2
                        # --- if not all zeros, set decimal depending on the
                        # maximum aboslute value of this variable so that the
                        # comparison has a reasonable precision. Specifically,
                        # decimal ~= - log10(max_abs_value) + 1 --- #
                        else:
                            decimal = int(round(- np.log10(np.max(np.absolute(
                                ds_image_cell[var].values))) + 1))
                        # --- keep decimal to be no greater than 4 --- #
                        if decimal > 4:
                            decimal = 4
                        # assert almost equal
                        np.testing.assert_almost_equal(
                            ds_image_cell[var].values, df_classic[var].values,
                            decimal=decimal,
                            err_msg='Variable {} is different in the classic '
                                    'and image drivers'.format(var))
                    # if [time, nlayer]
                    elif len(ds_image_cell[var].coords) == 4:
                        for l in ds_image['nlayer']:
                            s_classic = df_classic['{}_{}'.format(var,
                                                                  l.values)]
                            s_image = ds_image_cell[var].sel(
                                nlayer=l).to_series()
                            # determine precision for comparison
                            if np.mean(s_image.values) == 0:
                                decimal = 2
                            else:
                                decimal = int(round(- np.log10(np.max(
                                    np.absolute(s_image.values))) + 1))
                            if decimal > 4:
                                decimal = 4
                            # assert almost eqaul
                            np.testing.assert_almost_equal(
                                s_image.values,
                                s_classic.values,
                                decimal=decimal,
                                err_msg='Variable {} is different in '
                                        'the classic and image '
                                        'drivers'.format(var))


def tsplit(string, delimiters):
    '''Behaves like str.split but supports multiple delimiters. '''

    delimiters = tuple(delimiters)
    stack = [string]

    for delimiter in delimiters:
        for i, substring in enumerate(stack):
            substack = substring.split(delimiter)
            stack.pop(i)
            for j, _substring in enumerate(substack):
                stack.insert(i + j, _substring)

    return stack


def read_snotel_swe_obs(filename, science_test_data_dir, items):
    '''Reads in Snotel SWE obs and returns DataFrame. '''

    filename_fullpath = os.path.join(science_test_data_dir,
                                     'inputdata',
                                     items['archive'],
                                     'observations',
                                     filename)

    # load snotel obs
    snotel_swe = pd.read_csv(filename_fullpath,
                             skiprows=0,
                             delim_whitespace=True,
                             names=['YEAR', 'MONTH', 'DAY', 'OUT_SWE'])

    # add datetime index
    time_cols = ['YEAR', 'MONTH', 'DAY']
    snotel_swe.index = pd.to_datetime(snotel_swe[time_cols])

    # remove year, day columns of DataFrame
    snotel_swe.drop(time_cols, inplace=True, axis=1)

    return snotel_swe


def read_vic_42_output(lat, lng, science_test_data_dir, items):
    ''' Reads output from VIC 4.2. '''

    if items['compare_to'] == 'ecflux':
        vic_42_file = 'en_bal_%s_%s' % (lat, lng)
        vic_42_dir = os.path.join(science_test_data_dir, 'archive',
                                  items['archive'], 'ecflux', 'results')

    elif items['compare_to'] == 'snotel':
        vic_42_file = 'outfile_%s_%s' % (lat, lng)

        vic_42_dir = os.path.join(science_test_data_dir, 'archive',
                                  items['archive'], 'snotel', 'results')

    else:
        raise ValueError("this option (%s) has not yet been implemented"
                         % items['compare_to'])

    vic_42 = pd.read_csv(os.path.join(vic_42_dir, vic_42_file),
                         sep='\t',
                         skiprows=5)

    # remove comment sign from column names in DataFrame
    vic_42 = vic_42.rename(columns=lambda x: x.replace('#', ''))

    # remove spaces from column names in DataFrame
    vic_42 = vic_42.rename(columns=lambda x: x.replace(' ', ''))

    # rename radiation variables to be consistent with VIC 5
    if items['compare_to'] == 'ecflux':
        vic_42 = vic_42.rename(columns=lambda x: x.replace('OUT_NET_SHORT',
                                                           'OUT_SWNET'))
        vic_42 = vic_42.rename(columns=lambda x: x.replace('OUT_NET_LONG',
                                                           'OUT_LWNET'))

    # add datetime index
    time_cols = ['YEAR', 'MONTH', 'DAY']
    vic_42.index = pd.to_datetime(vic_42[time_cols])

    if 'HOUR' in vic_42:
        vic_42.index += pd.Series(
            [pd.Timedelta(s, unit='h') for s in vic_42['HOUR']],
            index=vic_42.index)
        time_cols.append('HOUR')

    # remove year, day columns of DataFrame
    vic_42.drop(time_cols, inplace=True, axis=1)

    return vic_42


def read_vic_5_output(lat, lng, result_dir, items):
    ''' Read VIC 5.0.x output. '''

    if items['compare_to'] == 'ecflux':
        vic_5_file = 'en_bal_%s_%s.txt' % (lat, lng)
        vic_5_dir = result_dir

    elif items['compare_to'] == 'snotel':
        vic_5_file = 'outfile_%s_%s.txt' % (lat, lng)
        vic_5_dir = result_dir

    else:
        raise ValueError("this option (%s) has not yet been implemented"
                         % items['compare_to'])

    vic_5 = pd.read_csv(os.path.join(vic_5_dir, vic_5_file),
                        skiprows=2,
                        sep='\t')

    # remove spaces from column names
    vic_5.rename(columns=lambda x: x.replace(' ', ''), inplace=True)

    # add datetime index
    time_cols = ['YEAR', 'MONTH', 'DAY']
    vic_5.index = pd.to_datetime(vic_5[time_cols])

    if 'SEC' in vic_5:
        vic_5.index += pd.Series([pd.Timedelta(s, unit='s')
                                  for s in vic_5['SEC']], index=vic_5.index)
        time_cols.append('SEC')

    # remove year, day columns of DataFrame
    vic_5.drop(time_cols, inplace=True, axis=1)

    return vic_5


def plot_science_tests(driver, test_type, science_test_data_dir, result_dir,
                       plot_dir, plots_to_make, compare_data, nproc):
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
        Keys that indicate which datasets and model output to use for
        comparison.
    nproc <int>
        Number of processors to use

    Returns
    ----------
    '''
    if test_type == "science_test_snotel":
        plot_snotel_comparison(driver,
                               science_test_data_dir,
                               compare_data,
                               result_dir,
                               plot_dir,
                               plots_to_make,
                               nproc)

    elif test_type == "science_test_fluxnet":
        plot_fluxnet_comparison(driver,
                                science_test_data_dir,
                                compare_data,
                                result_dir,
                                plot_dir,
                                plots_to_make,
                                nproc)
    else:
        raise ValueError("this option %s has not been implemented in the \
                         VIC 5.0 science test suite" % test_type)


def plot_snotel_comparison(driver, science_test_data_dir,
                           compare_data_dict,
                           result_dir, plot_dir,
                           plots_to_make, nproc):
    ''' makes snotel figures '''

    # plot settings
    plot_variables = {'OUT_SWE': 'mm', 'OUT_ALBEDO': 'fraction',
                      'OUT_SALBEDO': 'fraction', 'OUT_SNOW_DEPTH': 'mm',
                      'OUT_SNOW_CANOPY': '%', 'OUT_SNOW_PACK_TEMP':
                      'degrees C', 'OUT_SNOW_MELT': 'mm', 'OUT_R_NET':
                      '$W/{m^2}$', 'OUT_LATENT': '$W/{m^2}$',
                      'OUT_SENSIBLE': '$W/{m^2}$'}
    context = "paper"
    style = "whitegrid"

    # --- Set up multiprocessing --- #
    pool = mp.Pool(processes=nproc)

    for filename in os.listdir(os.path.join(science_test_data_dir,
                                            'inputdata',
                                            'snotel',
                                            'observations')):
        pool.apply_async(plot_snotel_comparison_one_site,
                         (driver, science_test_data_dir,
                          compare_data_dict,
                          result_dir, plot_dir,
                          plots_to_make,
                          plot_variables, context, style, filename,))

    # --- Finish multiprocessing --- #
    pool.close()
    pool.join()


def plot_snotel_comparison_one_site(
        driver, science_test_data_dir,
        compare_data_dict,
        result_dir, plot_dir,
        plots_to_make,
        plot_variables, context, style, filename):
    
    print(plots_to_make)
    
    # get lat/lng from filename
    file_split = re.split('_', filename)
    lng = file_split[3].split('.txt')[0]
    lat = file_split[2]
    print('Plotting {} {}'.format(lat, lng))

    # loop over data to compare
    data = {}
    for key, items in compare_data_dict.items():

        # read in data
        if key == "snotel":
            data[key] = read_snotel_swe_obs(filename,
                                            science_test_data_dir,
                                            items)

        elif key == "VIC.4.2.d":
            data[key] = read_vic_42_output(lat, lng,
                                           science_test_data_dir,
                                           items)

        else:
            data[key] = read_vic_5_output(lat, lng,
                                          result_dir,
                                          items)

    # loop over variables to plot
    for plot_variable, units in plot_variables.items():

        if 'water_year' in plots_to_make:

            with plt.rc_context(dict(sns.axes_style(style),
                                     **sns.plotting_context(context))):
                fig, ax = plt.subplots(figsize=(10, 10))

                df = pd.DataFrame({key: d[plot_variable] for key, d in
                                   data.items() if plot_variable in d})

                for key, series in df.iteritems():
                    series.plot(
                        use_index=True,
                        linewidth=compare_data_dict[key]['linewidth'],
                        ax=ax,
                        color=compare_data_dict[key]['color'],
                        linestyle=compare_data_dict[key]
                        ['linestyle'],
                        zorder=compare_data_dict[key]['zorder'])

                ax.legend(loc='upper left')
                ax.set_ylabel("%s [%s]" % (plot_variable, units))

                # save figure
                os.makedirs(os.path.join(plot_dir, plot_variable),
                            exist_ok=True)
                plotname = '%s_%s.png' % (lat, lng)
                savepath = os.path.join(plot_dir, plot_variable, plotname)
                plt.savefig(savepath, bbox_inches='tight')
                print(savepath)
                plt.clf()
                plt.close()


def check_site_files(obs_dir, subdir):
    return len(os.listdir(os.path.join(obs_dir, subdir))) > 0


def get_fluxnet_lat_lon(obs_dir, subdir):

    # get CSV file from site directory to get lat/lng for site
    try:
        site_csv_file = glob.glob(os.path.join(obs_dir, subdir, 'AMF*.csv'))[0]
    except IndexError:
        site_csv_file = glob.glob(os.path.join(obs_dir, subdir, 'us*.csv'))[0]
    with open(site_csv_file) as f:
        second_line = list(f)[1]

    # parse line from header to get lat/lng
    str_split = tsplit(second_line,
                       ('Latitude: ', 'Longitude: ', 'Elevation (masl): '))
    lat = str_split[1].strip()
    lng = str_split[2].strip()

    return lat, lng


def read_fluxnet_obs(subdir, science_test_data_dir, items):

    # column names for DataFrame (same as VIC variable names)
    fluxnet_names = ['YEAR', 'MONTH', 'DAY', 'HOUR', 'PREC', 'AIR_TEMP',
                     'SWDOWN', 'LWDOWN', 'OUT_REL_HUMID', 'PRESSURE', 'WIND',
                     'OUT_EVAP', 'SOIL_TEMP_DEPTH1', 'SOIL_TEMP_DEPTH2',
                     'SOIL_TEMP_DEPTH3', 'SOIL_TEMP_DEPTH4',
                     'SOIL_TEMP_DEPTH5', 'OUT_SOIL_MOIST1', 'OUT_SOIL_MOIST2',
                     'OUT_SOIL_MOIST3', 'OUT_SOIL_MOIST4', 'OUT_SOIL_MOIST5',
                     'OUT_SOIL_TEMP1', 'OUT_SOIL_TEMP2', 'OUT_SOIL_TEMP3',
                     'OUT_SOIL_TEMP4', 'OUT_SOIL_TEMP5', 'OUT_SWNET',
                     'OUT_LWNET', 'OUT_SENSIBLE', 'OUT_LATENT',
                     'OUT_GRND_FLUX']

    filename = '%s.stdfmt.hourly.local.txt' % subdir
    # read in data with -9999.0000 as NaNs
    obs_dir = os.path.join(science_test_data_dir, 'inputdata',
                           'ec_flux_towers', 'obs')
    ecflux_df = pd.read_csv(os.path.join(obs_dir, subdir, filename),
                            skiprows=0,
                            delim_whitespace=True,
                            header=None,
                            names=fluxnet_names,
                            na_values=-9999.0000)

    # add datetime index
    time_cols = ['YEAR', 'MONTH', 'DAY']
    ecflux_df.index = pd.to_datetime(ecflux_df[time_cols])

    if 'HOUR' in ecflux_df:
        ecflux_df.index += pd.Series(
            [pd.Timedelta(s, unit='h') for s in ecflux_df['HOUR']],
            index=ecflux_df.index)
        time_cols.append('HOUR')

    # remove year, day columns of DataFrame
    ecflux_df.drop(time_cols, inplace=True, axis=1)

    return ecflux_df


def plot_fluxnet_comparison(driver, science_test_data_dir,
                            compare_data_dict,
                            result_dir, plot_dir,
                            plots_to_make, nproc):
    ''' makes Ameriflux figures
    '''

    context = "paper"
    style = "whitegrid"
    var_names = {'OUT_LATENT': 'LH', 'OUT_SENSIBLE': 'H', 'OUT_SWNET':
                 'SW_NET', 'OUT_LWNET': 'LW NET'}

    months = ['January', 'February', 'March', 'April', 'May',
              'June', 'July', 'August', 'September',
              'October', 'November', 'December']

    # loop over Ameriflux sites
    obs_dir = os.path.join(science_test_data_dir,
                           'inputdata',
                           'ec_flux_towers',
                           'obs')

    # --- Set up multiprocessing --- #
    pool = mp.Pool(processes=nproc)

    for subdir in os.listdir(obs_dir):
        pool.apply_async(plot_fluxnet_comparison_one_site,
                         (driver, science_test_data_dir,
                          compare_data_dict,
                          result_dir, plot_dir,
                          plots_to_make,
                          context, style, var_names, months, obs_dir, subdir,))

    # --- Finish multiprocessing --- #
    pool.close()
    pool.join()


def plot_fluxnet_comparison_one_site(driver, science_test_data_dir,
                                     compare_data_dict, result_dir, plot_dir,
                                     plots_to_make, context, style, var_names,
                                     months, obs_dir, subdir):

    if check_site_files(obs_dir, subdir):
        # get CSV file from site directory to get lat/lng for site
        lat, lng = get_fluxnet_lat_lon(obs_dir, subdir)
        print(lat, lng)

        # loop over data to compare
        data = {}
        for key, items in compare_data_dict.items():

            if key == "ecflux":
                try:
                    # load Ameriflux data
                    data[key] = read_fluxnet_obs(subdir,
                                                 science_test_data_dir,
                                                 items)
                except OSError:
                    warnings.warn(
                        "this %s site does not have data" % subdir)

            elif key == "VIC.4.2.d":
                try:
                    # load VIC 4.2 simulations
                    data[key] = read_vic_42_output(lat, lng,
                                                   science_test_data_dir,
                                                   items)

                except OSError:
                    warnings.warn(
                        "this site has a lat/lng precision issue")

            else:
                try:
                    # load VIC 5 simulations
                    data[key] = read_vic_5_output(lat, lng,
                                                  result_dir,
                                                  items)
                except OSError:
                    warnings.warn(
                        "this site has a lat/lng precision issue")

        # make figures

        # plot preferences
        fs = 15
        dpi = 150

        if 'annual_mean_diurnal_cycle' in plots_to_make:

            # make annual mean diurnal cycle plots
            with plt.rc_context(dict(sns.axes_style(style),
                                     **sns.plotting_context(context))):
                f, axarr = plt.subplots(4, 1, figsize=(8, 8), sharex=True)

                for i, (vic_var, variable_name) in enumerate(
                        var_names.items()):

                    # calculate annual mean diurnal cycle for each
                    # DataFrame
                    annual_mean = {}
                    for key, df in data.items():
                        annual_mean[key] = pd.DataFrame(
                            df[vic_var].groupby(df.index.hour).mean())

                    df = pd.DataFrame(
                        {key: d[vic_var] for key, d in annual_mean.items()
                         if vic_var in d})

                    for key, series in df.iteritems():
                        series.plot(
                            linewidth=compare_data_dict[key]['linewidth'],
                            ax=axarr[i],
                            color=compare_data_dict[key]['color'],
                            linestyle=compare_data_dict[key]['linestyle'],
                            zorder=compare_data_dict[key]['zorder'])

                    axarr[i].legend(loc='upper left')
                    axarr[i].set_ylabel(
                        '%s ($W/{m^2}$)' % variable_name,
                        size=fs)
                    axarr[i].set_xlabel('Time of Day (Hour)', size=fs)
                    axarr[i].set_xlim([0, 24])
                    axarr[i].xaxis.set_ticks(np.arange(0, 24, 3))

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
            with plt.rc_context(dict(sns.axes_style(style),
                                     **sns.plotting_context(context))):
                f, axarr = plt.subplots(4, 12, figsize=(35, 7),
                                        sharex=True,
                                        sharey=True)

                for i, (vic_var, variable_name) in enumerate(
                        var_names.items()):

                    # calculate monthly mean diurnal cycle
                    monthly_mean = {}
                    for (key, df) in data.items():
                        monthly_mean[key] = pd.DataFrame(
                            df[vic_var].groupby([df.index.month,
                                                 df.index.hour]).mean())

                    df = pd.DataFrame(
                        {key: d[vic_var] for key, d in monthly_mean.items()
                         if vic_var in d})

                    for j, month in enumerate(months):

                        for key, series in df.iteritems():
                            series[j + 1].plot(
                                linewidth=compare_data_dict[key]['linewidth'],
                                ax=axarr[i, j],
                                color=compare_data_dict[key]['color'],
                                linestyle=compare_data_dict[key]['linestyle'],
                                zorder=compare_data_dict[key]['zorder'])

                        axarr[i, j].set_ylabel(
                            '%s \n ($W/{m^2}$)' % variable_name,
                            size=fs)
                        axarr[i, j].set_xlabel('', size=fs)
                        axarr[i, j].set_xlim([0, 24])
                        axarr[i, j].xaxis.set_ticks(np.arange(0, 24, 3))
                        if i == 0:
                            axarr[i, j].set_title(month, size=fs)

                # add legend
                axarr[0, -1].legend(loc='center left',
                                    bbox_to_anchor=(1, 0.5))

                # add common x label
                f.text(0.5, 0.04, 'Time of Day (Hour)', ha='center',
                       size=fs)

                # save plot
                plotname = '%s_%s.png' % (lat, lng)
                os.makedirs(os.path.join(plot_dir, 'monthly_mean'),
                            exist_ok=True)
                savepath = os.path.join(plot_dir,
                                        'monthly_mean', plotname)
                plt.savefig(savepath, bbox_inches='tight', dpi=dpi)

                plt.clf()
                plt.close()
