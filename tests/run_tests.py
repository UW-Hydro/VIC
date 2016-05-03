#!/usr/bin/env python
'''VIC testing command line interface'''

from __future__ import print_function
import os
import sys
import glob
import argparse
import datetime
from collections import OrderedDict
import string

import pytest
import pandas as pd
import numpy as np

from tonic.models.vic.vic import VIC, VICRuntimeError  # , read_vic_ascii
from tonic.io import read_config, read_configobj
from tonic.testing import check_completed, check_for_nans, VICTestError

test_dir = os.path.dirname(os.path.abspath(__file__))


OUTPUT_WIDTH = 100

description = '''
                            VIC Test Suite
-------------------------------------------------------------------------------
This is the VIC Test Suite. There are six main test types:

    1. unit: function level tests.
    2. system: tests that aim to address model runtime issues. These tests
            are generally very quick.
        * configuration errors - tests that address model startup and error
                checking.
        * restart: tests that address model state and restart capacity.
        * I/O: tests that address model input and output functionality.
            * forcings come out the way the come in
            * parameter files are appropriately read in and allocated.
    3. science: tests that aim to assess the model's scientific skill.
            Many of these tests are compared to observations of some kind.
    4. examples: a set of examples that users may download and run.
    5. release: longer, full domain simulations performed prior to release
            demonstrating model output for a final release.
-------------------------------------------------------------------------------
'''

epilog = '''
-------------------------------------------------------------------------------
For questions about the development or use of VIC or use of this test module,
please email the VIC users list serve at vic_users@u.washington.edu.
-------------------------------------------------------------------------------
'''


class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter,
                      argparse.RawDescriptionHelpFormatter):
    pass


class TestResults(object):
    def __init__(self, name, test_complete=False, passed=False,
                 comment='', error_message='', returncode=None):
        self.name = name
        self.test_complete = test_complete
        self.passed = passed
        self.comment = comment
        self.error_message = error_message
        self.returncode = returncode

    def __repr__(self):
        r = '''
{0} test results:
    Passed: {1}
    Comment:{2}
    Return Code: {3}
    '''.format(self.name, self.passed, self.comment, self.returncode)
        return r

    def __str__(self):
        return '{0} test results:Passed: {1}, Comment:{2}'.format(self.name,
                                                                  self.passed,
                                                                  self.comment)


# -------------------------------------------------------------------- #
class VICReturnCodeError(BaseException):
    pass
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
def main():
    '''
    Run VIC tests
    '''

    # ---------------------------------------------------------------- #
    # dates and times
    starttime = datetime.datetime.now()
    ymd = starttime.strftime('%Y%m%d')
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Parse arguments
    test_results = OrderedDict()

    parser = argparse.ArgumentParser(description=description, epilog=epilog,
                                     formatter_class=CustomFormatter)

    parser.add_argument('tests', type=str,
                        help='Test sets to run',
                        choices=['all', 'unit', 'system', 'science',
                                 'examples', 'release'],
                        default=['unit', 'system'], nargs='+')
    parser.add_argument('--system', type=str,
                        help='system tests configuration file',
                        default='./system/system_tests.cfg')
    parser.add_argument('--science', type=str,
                        help='science tests configuration file',
                        default='./science/science.cfg')
    parser.add_argument('--examples', type=str,
                        help='examples tests configuration file',
                        default='./examples/examples.cfg')
    parser.add_argument('--release', type=str,
                        help='release tests configuration file',
                        default='./release/release.cfg')
    parser.add_argument('--vic_exe', type=str,
                        help='VIC executable to test',
                        default='../vic/drivers/classic/vic_classic.exe')
    parser.add_argument('--driver', type=str,
                        help='VIC driver to test',
                        choices=['classic', 'image'],
                        default='classic')
    parser.add_argument('--output_dir', type=str,
                        help='directory to write test output to',
                        default='$WORKDIR/VIC_tests_{0}'.format(ymd))
    parser.add_argument('--data_dir', type=str,
                        help='directory to find test data',
                        default='./samples/VIC_sample_data')
    args = parser.parse_args()
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Define test directories
    data_dir = args.data_dir
    out_dir = os.path.expandvars(args.output_dir)
    os.makedirs(out_dir, exist_ok=True)
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Validate input directories
    if not (len(args.tests) == 1 and args.tests[0] == 'unit'):
        for d in [data_dir, test_dir]:
            if not os.path.exists(d):
                raise VICTestError('Directory: {0} does not exist'.format(d))
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Print welcome information
    print(description)
    print('\nStarting tests now...Start Time: {0}\n'.format(starttime))
    print('Running Test Set: {0}'.format(', '.join(args.tests)))
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Setup VIC executable
    if not (len(args.tests) == 1 and args.tests[0] == 'unit'):
        vic_exe = VIC(args.vic_exe)
        print('VIC version string:\n{0}'.format(vic_exe.version.decode()))
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # run test sets
    # unit
    if any(i in ['all', 'unit'] for i in args.tests):
        test_results['unit'] = run_unit_tests()

    # system
    if any(i in ['all', 'system'] for i in args.tests):
        test_results['system'] = run_system(args.system, vic_exe, data_dir,
                                            os.path.join(out_dir, 'system'),
                                            args.driver)
    # science
    if any(i in ['all', 'science'] for i in args.tests):
        test_results['science'] = run_science(args.science, vic_exe, data_dir,
                                              os.path.join(out_dir, 'science'))
    # examples
    if any(i in ['all', 'examples'] for i in args.tests):
        test_results['examples'] = run_examples(args.examples, vic_exe, data_dir,
                                                os.path.join(out_dir, 'examples'),
                                                args.driver)
    # release
    if any(i in ['all', 'release'] for i in args.tests):
        test_results['release'] = run_release(args.release)
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Print test results
    summary = OrderedDict()
    failed = 0
    print('\nTest Results:')
    for test_set, results in test_results.items():
        print('-'.ljust(OUTPUT_WIDTH, '-'))
        print(test_set.center(OUTPUT_WIDTH))
        print('-'.ljust(OUTPUT_WIDTH, '-'))
        print_test_dict(results)

        summary[test_set] = 0
        for r in results.values():
            if not r.passed:
                summary[test_set] += 1

        failed += summary[test_set]

    print('\nTest Summary:')
    print('-'.ljust(OUTPUT_WIDTH, '-'))
    for test_set, r in summary.items():
        print('Failed tests in {0}: {1}'.format(test_set, r))
    print('-'.ljust(OUTPUT_WIDTH, '-'))
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # end date and times
    endtime = datetime.datetime.now()
    elapsed = endtime - starttime
    print('\nFinished testing VIC. Endtime: {0}'.format(endtime))
    print('Time elapsed during testing:  {0}\n'.format(elapsed))
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # return exit code
    sys.exit(failed)
    # ---------------------------------------------------------------- #

    return
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
def run_unit_tests():
    '''Run unittests from config file'''

    print('\n-'.ljust(OUTPUT_WIDTH + 1, '-'))
    print('Running Unit Tests')
    print('-'.ljust(OUTPUT_WIDTH, '-'))

    retcode = pytest.main(['-x',  os.path.join(test_dir, 'unit')])
    return {'unittests': TestResults('unittests',
                                     test_complete=True,
                                     passed=retcode == 0,
                                     comment='see stdout from pytest',
                                     returncode=retcode)}

    print('\n-'.ljust(OUTPUT_WIDTH + 1, '-'))
    print('Finished unit tests.')
    print('-'.ljust(OUTPUT_WIDTH, '-'))
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
def run_system(config_file, vic_exe, test_data_dir, out_dir, driver):
    '''Run system tests from config file'''

    # ---------------------------------------------------------------- #
    # Print test set welcome
    print('\n-'.ljust(OUTPUT_WIDTH + 1, '-'))
    print('Running System Tests')
    print('-'.ljust(OUTPUT_WIDTH, '-'))
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Get setup
    config = read_configobj(config_file)
    test_results = OrderedDict()
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Run individual system tests
    for i, (testname, test_dict) in enumerate(config.items()):

        # ------------------------------------------------------------ #
        # print out status info
        print('Running test {0}/{1}: {2}'.format(i + 1, len(config.items()),
                                                 testname))
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # Setup directories for test
        dirs = setup_test_dirs(testname, out_dir,
                               mkdirs=['results', 'state', 'logs', 'plots'])
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # read template global parameter file
        infile = os.path.join(test_dir, 'system',
                              test_dict['global_parameter_file'])

        with open(infile, 'r') as global_file:
            global_param = global_file.read()
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # If restart test, prepare running periods
        # (1) Find STATESEC option (and STATE_FORMAT option for later use)
        statesec = find_global_param_value(global_param, 'STATESEC')
        state_format = find_global_param_value(global_param, 'STATE_FORMAT')
        # (2) Prepare running periods and initial state file info for restart test
        if 'restart' in test_dict:
            run_periods = prepare_restart_run_periods(
                                test_dict['restart'],
                                dirs['state'], statesec)
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # create template string
        s = string.Template(global_param)
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # fill in global parameter options
        #--- if restart test, multiple runs ---#
        if 'restart' in test_dict:
            # Set up subdirectories and fill in global parameter options 
            # for restart testing
            list_global_param = setup_subdirs_and_fill_in_global_param_restart_test(
                s, run_periods, driver, dirs['results'], dirs['state'], test_data_dir)

        #--- Else, single run ---#
        else:
            global_param = s.safe_substitute(test_data_dir=test_data_dir,
                                             result_dir=dirs['results'],
                                             state_dir=dirs['state'])
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # replace global options from config file
        # extract global options to be substitute
        if 'options' in test_dict:
            replacements = test_dict['options']
        else:
            replacements = OrderedDict()

        # if STATE_FORMAT is specified, te the specified value (instead of 
        # the one in the global template file)
        if 'STATE_FORMAT' in replacements:
            state_format = replacements.pop('STATE_FORMAT')

        # replace global options
        if 'restart' in test_dict:
            for j, gp in enumerate(list_global_param):
                list_global_param[j] = replace_global_values(gp, replacements)
        else:
            global_param = replace_global_values(global_param, replacements)
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # write global parameter file
        if 'restart' in test_dict:
            list_test_global_file = []
            for j, gp in enumerate(list_global_param):
                test_global_file = os.path.join(
                            dirs['test'],
                            '{}_globalparam_{}_{}.txt'.format(
                                    testname,
                                    run_periods[j]['start_date'].strftime("%Y%m%d"),
                                    run_periods[j]['end_date'].strftime("%Y%m%d")))
                list_test_global_file.append(test_global_file)
                with open(test_global_file, mode='w') as f:
                    for line in gp:
                        f.write(line)
        else:
            test_global_file = os.path.join(dirs['test'],
                                            '{0}_globalparam.txt'.format(testname))
            with open(test_global_file, mode='w') as f:
                for line in global_param:
                    f.write(line)
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # run VIC
        test_complete = False
        test_passed = False
        test_comment = ''
        error_message = ''

        try:
            expectation = int(test_dict['expected_retval'])
            returncode_flag = True

            if 'restart' in test_dict:
                for j, test_global_file in enumerate(list_test_global_file):
                    returncode = vic_exe.run(test_global_file, logdir=dirs['logs'])
                    if returncode != expectation:
                        returncode_flag = False
            else:
                returncode = vic_exe.run(test_global_file, logdir=dirs['logs'])
                if returncode != expectation:
                    returncode_flag = False
            test_complete = True

            if returncode_flag == False:
                raise VICReturnCodeError('VIC return code ({0}) '
                                         'does not match expected '
                                         '({1})'.format(returncode,
                                                        expectation))

            # -------------------------------------------------------- #
            # check output files
            if test_dict['check']:
                if 'complete' in test_dict['check']:
                    if 'restart' in test_dict:
                        pass
                    else:  # single run
                        start = None
                        end = None
                        for fname in glob.glob(os.path.join(dirs['results'], '*')):
                            df = read_vic_ascii(fname, header=True)
                            # check that each dataframe includes all timestamps
                            if 'complete' in test_dict['check']:
                                if (start is not None) and (end is not None):
                                    check_completed(df, start, end)
                                else:
                                    start = df.index[0]
                                    end = df.index[-1]

                # check for nans in the df
                if 'nonans' in test_dict['check']:
                    if 'restart' in test_dict:
                        pass
                    else:  # single run
                        for fname in glob.glob(os.path.join(dirs['results'], '*')):
                            df = read_vic_ascii(fname, header=True)
                            check_for_nans(df)

                # check for exact restarts
                if 'exact_restart' in test_dict['check']:
                    #check_exact_restart_fluxes(dirs['results'], driver, run_periods)
                    check_exact_restart_states(dirs['state'], driver,
                                               run_periods, statesec, state_format)
            # -------------------------------------------------------- #

            # -------------------------------------------------------- #
            # if we got this far, the test passed.
            test_passed = True
            # -------------------------------------------------------- #

        # Handle errors
        except VICRuntimeError as e:
            test_comment = 'Test failed during simulation'
            error_message = e
        except VICTestError as e:
            test_comment = 'Test failed during testing of output files'
            error_message = e
        except VICReturnCodeError as e:
            test_comment = 'Test failed due to incorrect return code'
            error_message = e

        if test_comment or error_message:
            print('\t{0}'.format(test_comment))
            print('\t{0}'.format(error_message))
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # record the test results
        test_results[testname] = TestResults(testname,
                                             test_complete=test_complete,
                                             passed=test_passed,
                                             comment=test_comment,
                                             error_message=error_message,
                                             returncode=returncode)
        # ------------------------------------------------------------ #

    print('-'.ljust(OUTPUT_WIDTH, '-'))
    print('Finished testing system tests.')
    print('-'.ljust(OUTPUT_WIDTH, '-'))
    # ---------------------------------------------------------------- #

    return test_results
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
def run_science(config_file, vic_exe, test_data_dir, out_dir):
    '''Run science tests from config file'''

    # ---------------------------------------------------------------- #
    # Print test set welcome
    print('\n-'.ljust(OUTPUT_WIDTH + 1, '-'))
    print('Running Science Tests')
    print('-'.ljust(OUTPUT_WIDTH, '-'))
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Get setup
    config = read_config(config_file)

    test_results = OrderedDict()
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Run individual tests
    for i, (testname, test_dict) in enumerate(config.items()):

        # ------------------------------------------------------------ #
        # print out status info
        print('Running test {0}/{1}: {2}'.format(i + 1, len(config.items()),
                                                 testname))
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # Setup directories for test
        dirs = setup_test_dirs(testname, out_dir,
                               mkdirs=['results', 'state', 'logs', 'plots'])
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # read template global parameter file
        infile = os.path.join(test_dir, test_dict['global_parameter_file'])

        with open(infile, 'r') as global_file:
            global_param = global_file.read()
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # create template string
        s = string.Template(global_param)
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # fill in global parameter options
        global_param = s.safe_substitute(test_data_dir=test_data_dir,
                                         result_dir=dirs['results'],
                                         state_dir=dirs['state'],
                                         testname=testname,
                                         test_root=test_dir)

        test_global_file = os.path.join(dirs['test'],
                                        '{0}_globalparam.txt'.format(testname))
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # write global parameter file
        with open(test_global_file, 'w') as f:
            f.write(global_param)
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # run VIC
        test_complete = False
        test_passed = False
        test_comment = ''
        error_message = ''

        try:
            returncode = vic_exe.run(test_global_file, logdir=dirs['logs'])
            test_complete = True

            if returncode != 0:
                raise VICReturnCodeError('VIC return code not equal to zero, '
                                         'check VIC logs for test')

            # -------------------------------------------------------- #
            # check output files
            if test_dict['check']:
                if 'complete' in test_dict['check']:
                    start = None
                    end = None

                for fname in glob.glob(os.path.join(dirs['results'], '*')):
                    df = read_vic_ascii(fname, header=True)
                    # do numeric tests

                    # check that each dataframe includes all timestamps
                    if 'complete' in test_dict['check']:
                        if (start is not None) and (end is not None):
                            check_completed(df, start, end)
                        else:
                            start = df.index[0]
                            end = df.index[-1]

                    # check for nans in the df
                    if 'nonans' in test_dict['check']:
                        check_for_nans(df)
            # -------------------------------------------------------- #

            # -------------------------------------------------------- #
            # if we got this far, the test passed.
            test_passed = True
            # -------------------------------------------------------- #

        # Handle errors
        except VICRuntimeError as e:
            test_comment = 'Test failed during simulation'
            error_message = e
        except VICTestError as e:
            test_comment = 'Test failed during testing of output files'
            error_message = e
        except VICReturnCodeError as e:
            test_comment = 'Test failed due to incorrect return code'
            error_message = e

        if test_comment or error_message:
            print('\t{0}'.format(test_comment))
            print('\t{0}'.format(error_message))
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # record the test results
        test_results[testname] = TestResults(testname,
                                             test_complete=test_complete,
                                             passed=test_passed,
                                             comment=test_comment,
                                             error_message=error_message,
                                             returncode=returncode)
        # ------------------------------------------------------------ #

    print('-'.ljust(OUTPUT_WIDTH, '-'))
    print('Finished testing science tests.')
    print('-'.ljust(OUTPUT_WIDTH, '-'))
    # ---------------------------------------------------------------- #

    return test_results
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
def run_examples(config_file, vic_exe, test_data_dir, out_dir, driver):
    '''Run examples tests from config file '''

    # ---------------------------------------------------------------- #
    # Print test set welcome
    print('\n-'.ljust(OUTPUT_WIDTH + 1, '-'))
    print('Running Examples')
    print('-'.ljust(OUTPUT_WIDTH, '-'))
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Get setup
    config = read_config(config_file)

    # drop invalid driver tests
    drop_tests = []
    for key, test_cfg in config.items():
        if test_cfg['driver'].lower() != driver.lower():
            drop_tests.append(key)
    for test in drop_tests:
        del config[test]

    test_results = OrderedDict()
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Run individual examples
    for i, (testname, test_dict) in enumerate(config.items()):

        # ------------------------------------------------------------ #
        # print out status info
        print('Running test {0}/{1}: {2}'.format(i + 1, len(config.items()),
                                                 testname))
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # Setup directories for test
        dirs = setup_test_dirs(testname, out_dir,
                               mkdirs=['results', 'state', 'logs', 'plots'])
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # read template global parameter file
        infile = os.path.join(test_dir, 'examples', test_dict['global_parameter_file'])

        with open(infile, 'r') as global_file:
            global_param = global_file.read()
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # create template string
        s = string.Template(global_param)
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # fill in global parameter options
        global_param = s.safe_substitute(test_data_dir=test_data_dir,
                                         result_dir=dirs['results'],
                                         state_dir=dirs['state'],
                                         testname=testname,
                                         test_root=test_dir)

        test_global_file = os.path.join(dirs['test'],
                                        '{0}_globalparam.txt'.format(testname))
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # write global parameter file
        with open(test_global_file, 'w') as f:
            f.write(global_param)
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # run VIC
        test_complete = False
        test_passed = False
        test_comment = ''
        error_message = ''

        try:
            returncode = vic_exe.run(test_global_file, logdir=dirs['logs'])
            test_complete = True

            if returncode != 0:
                raise VICReturnCodeError('VIC return code not equal to zero, '
                                         'check VIC logs for test')

            # -------------------------------------------------------- #
            # check output files
            if test_dict['check']:
                if 'complete' in test_dict['check']:
                    start = None
                    end = None

                for fname in glob.glob(os.path.join(dirs['results'], '*')):
                    df = read_vic_ascii(fname, header=True)
                    # do numeric tests

                    # check that each dataframe includes all timestamps
                    if 'complete' in test_dict['check']:
                        if (start is not None) and (end is not None):
                            check_completed(df, start, end)
                        else:
                            start = df.index[0]
                            end = df.index[-1]

                    # check for nans in the df
                    if 'nonans' in test_dict['check']:
                        check_for_nans(df)
            # -------------------------------------------------------- #

            # -------------------------------------------------------- #
            # if we got this far, the test passed.
            test_passed = True
            # -------------------------------------------------------- #

        # Handle errors
        except VICRuntimeError as e:
            test_comment = 'Test failed during simulation'
            error_message = e
            tail = vic_exe.stderr
        except VICTestError as e:
            test_comment = 'Test failed during testing of output files'
            error_message = e
            tail = None
        except VICReturnCodeError as e:
            test_comment = 'Test failed due to incorrect return code'
            error_message = e
            tail = vic_exe.stderr

        if test_comment or error_message:
            print('\t{0}'.format(test_comment))
            print('\t{0}'.format(error_message))
            if tail is not None:
                print('\tLast 10 lines of standard out:')
                print_tail(tail)
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # record the test results
        test_results[testname] = TestResults(testname,
                                             test_complete=test_complete,
                                             passed=test_passed,
                                             comment=test_comment,
                                             error_message=error_message,
                                             returncode=returncode)
        # ------------------------------------------------------------ #

    print('-'.ljust(OUTPUT_WIDTH, '-'))
    print('Finished testing examples.')
    print('-'.ljust(OUTPUT_WIDTH, '-'))
    # ---------------------------------------------------------------- #

    return test_results
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
def run_release(config_file):
    '''Run release from config file'''
    return OrderedDict()
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
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
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
def print_test_dict(d):
    '''print a nicely formatated set of test results'''
    print('{0: <48} | {1: <6} | {2}'.format('Test Name', 'Passed', 'Comment'))
    print('-'.ljust(OUTPUT_WIDTH, '-'))
    for k, v in d.items():
        print('{0: <48} | {1: <6} | {2}'.format(v.name, str(bool(v.passed)),
                                                v.comment))
        print('-'.ljust(OUTPUT_WIDTH, '-'))
    return
# -------------------------------------------------------------------- #


def print_tail(string, n=10, indent='\t--->'):
    lines = string.decode().splitlines()
    for l in lines[-10:]:
        print('{0}{1}'.format(indent, l))


# -------------------------------------------------------------------- #
def replace_global_values(gp, replace):
    '''given a multiline string that represents a VIC global parameter file,
       loop through the string, replacing values with those found in the
       replace dictionary'''
    gpl = []
    for line in iter(gp.splitlines()):
        line_list = line.split()
        if line_list==[]:
            continue
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
# -------------------------------------------------------------------- #


# TODO: Update tonic version of this function, need to check that subdaily works
# -------------------------------------------------------------------- #
def read_vic_ascii(filepath, header=True, parse_dates=True,
                   datetime_index=None, names=None, **kwargs):
    """Generic reader function for VIC ASCII output with a standard header
    filepath: path to VIC output file
    header (True or False):  Standard VIC header is present
    parse_dates (True or False): Parse dates from file
    datetime_index (Pandas.tseries.index.DatetimeIndex):  Index to use as
    datetime index names (list like): variable names
    **kwargs: passed to Pandas.read_table
    returns Pandas.DataFrame
    """
    kwargs['header'] = None

    if header:
        kwargs['skiprows'] = 6

        # get names
        if names is None:
            with open(filepath) as f:
                # skip lines 0 through 3
                for _ in range(3):
                    next(f)

                # process header
                names = next(f)
                names = names.strip('#').replace('OUT_', '').split()

    kwargs['names'] = names

    if parse_dates:
        time_cols = ['YEAR', 'MONTH', 'DAY']
        if 'SECONDS' in names:
            time_cols.append('SECONDS')
        kwargs['parse_dates'] = {'datetime': time_cols}
        kwargs['index_col'] = 0

    df = pd.read_table(filepath, **kwargs)

    if datetime_index is not None:
        df.index = datetime_index

    return df
# -------------------------------------------------------------------- #

# -------------------------------------------------------------------- #
def prepare_restart_run_periods(restart_dict, state_basedir, statesec):
    ''' For restart tests, read full running period and splitting dates into datetime objects
    Parameters
    ----------
    restart_dict: <class 'configobj.Section'>
        A section of the config file for exact restart test setup
    state_basedir: <str>
        Basedirectory of output state files.
        For restart tests, state files will be output as <state_basedir>/<run_start_date>_<run_end_date>/<state_file>
    statesec: <int>
        STATESEC option in global parameter file

    Returns
    ----------
    run_periods: OrderedDict
        A list of running periods, including the full-period run, and all splitted 
        runs in order. Each element of run_period is a dictionary with keys:
            start_date
            end_date
            init_state  # None, or full path of the initial state file
                        # e.g., '/path/19490101_19490105/states_19490105_82800'
        
    '''

    #--- Read in full running period ---#
    start_date = datetime.datetime(restart_dict['start_date'][0],
                                   restart_dict['start_date'][1],
                                   restart_dict['start_date'][2])
    end_date = datetime.datetime(restart_dict['end_date'][0],
                                 restart_dict['end_date'][1],
                                 restart_dict['end_date'][2])
    #--- Identify each of the splitted running period ---#
    list_split_dates = []
    for i, split in enumerate(restart_dict['split_dates']):
        list_split_dates.append(datetime.datetime(restart_dict['split_dates'][split][0],
                                                  restart_dict['split_dates'][split][1],
                                                  restart_dict['split_dates'][split][2]))
    #--- Prepare running periods ---#
    # run_periods is a list of running periods, including the full-period run, 
    # and all splitted runs in order. Each element of run_period is a dictionary 
    # with keys:
    #       start_date
    #       end_date
    #       init_state  # None, or full path of the initial state file
    #                   # e.g., '/path/19490101_19490105/states_19490105_82800'
    run_periods = []
    # Append the full run
    d = dict(start_date=start_date, end_date=end_date)
    d['init_state'] = None
    run_periods.append(d)
    # First splitted running period - start_date to first split date
    d = dict(start_date=start_date, end_date=list_split_dates[0])
    d['init_state'] = None
    run_periods.append(d)
    # Loop over each of the rest splitted periods
    for i in range(len(list_split_dates)-1):
        d = dict(start_date=list_split_dates[i] + dt.timedelta(days=1),
                 end_date=list_split_dates[i+1])
        d['init_state'] = os.path.join(
                state_basedir,
                '{}_{}'.format(run_periods[-1]['start_date'].strftime("%Y%m%d"),
                               run_periods[-1]['end_date'].strftime("%Y%m%d")),
                '{}{}_{}'.format('states_',
                                 run_periods[-1]['end_date'].strftime("%Y%m%d"),
                                 statesec))
    run_periods.append(d)
    # Last splitted running period - last split date to end_date
    d = dict(start_date=list_split_dates[len(list_split_dates)-1] + 
                        datetime.timedelta(days=1),
             end_date=end_date)
    d['init_state'] = os.path.join(
            state_basedir,
            '{}_{}'.format(run_periods[-1]['start_date'].strftime("%Y%m%d"),
                           run_periods[-1]['end_date'].strftime("%Y%m%d")),
            '{}{}_{}'.format('states_',
                             run_periods[-1]['end_date'].strftime("%Y%m%d"),
                             statesec))
    run_periods.append(d)

    return run_periods
# -------------------------------------------------------------------- #

# -------------------------------------------------------------------- #
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
        if line_list==[]:
            continue
        key = line_list[0]
        if key==param_name:
            return line_list[1]
# -------------------------------------------------------------------- #

# -------------------------------------------------------------------- #
def setup_subdirs_restart_test(result_basedir, state_basedir, run_periods):
    ''' Set up subdirectories for multiple runs for restart testing
    Parameters
    ----------
    result_basedir: <str>
        Base directory of output fluxes results; running periods are subdirectories under the base directory
    state_basedir: <str>
        Base directory of output state results; running periods are subdirectories under the base directory
    run_periods: <list>
        A list of running periods. Return from prepare_restart_run_periods()

    Returns
    ----------
    None
    '''

    for j, run_period in enumerate(run_periods):
        # Set up subdirectories for results and states
        run_start_date = run_period['start_date']
        run_end_date = run_period['end_date']
        result_dir = os.path.join(result_basedir,
                                  '{}_{}'.format(run_start_date.strftime("%Y%m%d"),
                                                 run_end_date.strftime("%Y%m%d")))
        state_dir = os.path.join(state_basedir,
                                 '{}_{}'.format(run_start_date.strftime("%Y%m%d"),
                                                run_end_date.strftime("%Y%m%d")))
        os.makedirs(result_dir, exist_ok=True)
        os.makedirs(state_dir, exist_ok=True)

# -------------------------------------------------------------------- #

# -------------------------------------------------------------------- #
def setup_subdirs_and_fill_in_global_param_restart_test(s, run_periods, driver, result_basedir, state_basedir, test_data_dir):
    ''' Fill in global parameter options for multiple runs for restart testing

    Parameters
    ----------
    s: <string.Template>
        Template of the global param file to be filled in
    run_periods: <list>
        A list of running periods. Return from prepare_restart_run_periods()
    driver: <str>
        'classic' or 'image'
    result_basedir: <str>
        Base directory of output fluxes results; running periods are subdirectories under the base directory
    state_basedir: <str>
        Base directory of output state results; running periods are subdirectories under the base directory
    test_data_dir: <str>
        Base directory of test data

    Returns
    ----------
    '''

    list_global_param = []
    for j, run_period in enumerate(run_periods):
        # Set up subdirectories for results and states
        run_start_date = run_period['start_date']
        run_end_date = run_period['end_date']
        result_dir = os.path.join(result_basedir,
                                  '{}_{}'.format(run_start_date.strftime("%Y%m%d"),
                                                 run_end_date.strftime("%Y%m%d")))
        state_dir = os.path.join(state_basedir,
                                 '{}_{}'.format(run_start_date.strftime("%Y%m%d"),
                                                run_end_date.strftime("%Y%m%d")))
        os.makedirs(result_dir, exist_ok=True)
        os.makedirs(state_dir, exist_ok=True)
        # Determine initial state
        if run_period['init_state']==None: # if no initial state
            init_state='#INIT_STATE'
        else: # else, the initial state is the last time step
            init_state='INIT_STATE {}'.format(run_period['init_state'])
            if driver=='image': # In image driver, the name of
                                # the state file is 'basepath.*'
                                # instead of 'basepath_*', and
                                # ends with ".nc"
                init_state = init_state.replace("states_", "states.") + '.nc'

        # Fill in global parameter options
        list_global_param.append(s.safe_substitute(
                test_data_dir=test_data_dir,
                result_dir=result_dir,
                state_dir=state_dir,
                startyear=run_start_date.year,
                startmonth=run_start_date.month,
                startday=run_start_date.day,
                endyear=run_end_date.year,
                endmonth=run_end_date.month,
                endday=run_end_date.day,
                init_state=init_state,
                stateyear=run_end_date.year,
                statemonth=run_end_date.month,
                stateday=run_end_date.day))
    return(list_global_param)
# -------------------------------------------------------------------- #

# -------------------------------------------------------------------- #
def check_exact_restart_fluxes(result_basedir, driver, run_periods):
    ''' Checks whether all the fluxes are the same w/ or w/o restart

    Parameters
    ----------
    result_basedir: <str>
        Base directory of output fluxes results; running periods are subdirectories under the base directory
    driver: <str>
        'classic' or 'image'
    run_periods: <list>
        A list of running periods. Return from prepare_restart_run_periods()

    Returns
    ----------
    None
    '''

    #--- Extract full run period ---#
    run_full_start_date = run_periods[0]['start_date']
    run_full_end_date = run_periods[0]['end_date']
    #--- Read full run fluxes ---#
    if driver=='classic':
        result_dir = os.path.join(
            result_basedir,
            '{}_{}'.format(run_full_start_date.strftime('%Y%m%d'),
                           run_full_end_date.strftime('%Y%m%d')))
        # Read in each of the output flux files
        dict_df_full_run = {}  # a dict of flux at each grid cell, keyed by flux basename
        for fname in glob.glob(os.path.join(result_dir, '*')):
            df = read_vic_ascii(fname, header=True)
            dict_df_full_run[os.path.basename(fname)] = df

    #--- Loop over the result of each split run period ---#
    for i, run_period in enumerate(run_periods):
        # Skip the full run
        if i==0:
            continue
        # Extract running period
        start_date = run_period['start_date']
        end_date = run_period['end_date']
        # Loop over each of the output flux files
        result_dir = os.path.join(
            result_basedir,
            '{}_{}'.format(start_date.strftime('%Y%m%d'),
                           end_date.strftime('%Y%m%d')))
        if driver=='classic':
            for flux_basename in dict_df_full_run.keys():
                # Read in flux data
                fname = os.path.join(result_dir, flux_basename)
                df = read_vic_ascii(fname, header=True)
                # Extract the same period from the full run
                df_full_run_split_period = dict_df_full_run[flux_basename].truncate(
                                                         before=start_date,
                                                         after=end_date)
                # Compare split run fluxes with full run
                df_diff = df - df_full_run_split_period
                if np.absolute(df_diff).max().max() > 0:
                    raise VICTestError('Restart causes inexact flux outputs '
                                       'for running period {} - {} at grid cell {}!'.
                                format(start_date.strftime('%Y%m%d'),
                                       end_date.strftime('%Y%m%d'),
                                       flux_basename))
                else:
                    continue
    return

# -------------------------------------------------------------------- #

# -------------------------------------------------------------------- #
def check_exact_restart_states(state_basedir, driver, run_periods, statesec, state_format='ASCII'):
    ''' Checks whether all the states are the same w/ or w/o restart.
        Only test the state at the last time step.

    Parameters
    ----------
    state_basedir: <str>
        Base directory of output state results; running periods are subdirectories under the base directory
    driver: <str>
        'classic' or 'image'
    run_periods: <list>
        A list of running periods. Return from prepare_restart_run_periods()
    statesec: <int>
        STATESEC option in global parameter file
    state_format: <str>
        state file format, 'ASCII' or 'BINARY'; only need to specify when driver=='classic'

    Returns
    ----------
    None
    '''

    #--- Read the state at the end of the full run ---#
    # Extract full run period
    run_full_start_date = run_periods[0]['start_date']
    run_full_end_date = run_periods[0]['end_date']
    # Read the state file
    if driver=='classic':
        state_fname = os.path.join(
            state_basedir,
            '{}_{}'.format(run_full_start_date.strftime('%Y%m%d'),
                           run_full_end_date.strftime('%Y%m%d')),
            'states_{}_{}'.format(run_full_end_date.strftime('%Y%m%d'), statesec))
        if state_format=='ASCII':
            states_full_run = read_ascii_state(state_fname)
        elif state_format=='BINARY':
            states_full_run = read_binary_state(state_fname)

    #--- Read the state at the end of the last period of run ---#
    # Extract the last split run period
    run_last_period_start_date = run_periods[-1]['start_date']
    run_last_period_end_date = run_periods[-1]['end_date']
    # Read the state file
    if driver=='classic':
        state_fname = os.path.join(
            state_basedir,
            '{}_{}'.format(run_last_period_start_date.strftime('%Y%m%d'),
                           run_last_period_end_date.strftime('%Y%m%d')),
            'states_{}_{}'.format(run_last_period_end_date.strftime('%Y%m%d'), statesec))
        if state_format=='ASCII':
            states = read_ascii_state(state_fname)
        elif state_format=='BINARY':
            states = read_binary_state(state_fname)

    #--- Compare split run states with full run ---#
    if driver=='classic':
        # If ASCII state file, check if almost the same
        if state_format=='ASCII':
            states_diff = states - states_full_run
            if np.absolute(states_diff).max() > pow(10, -6):
                raise VICTestError('Restart causes inexact state outputs!')
            else:
                return
        # If BINARY state file, check if exactly the same
        elif state_format=='BINARY':
            if states!=states_full_run:
                raise VICTestError('Restart causes inexact state outputs!')
            else:
                return

# -------------------------------------------------------------------- #

# -------------------------------------------------------------------- #
def read_ascii_state(state_fname):
    ''' Read in ascii format state file and convert to a list of numbers

    Parameters
    ----------
    state_fname: <str>
        Path of the state file to be read

    Returns
    ----------
    states: <np.array>
        A np.array of float numbers of the state file
    '''

    with open(state_fname, 'r') as f:
        list_states = f.read().split()
        for i, item in enumerate(list_states):
            list_states[i] = float(item)
    return np.asarray(list_states)
# -------------------------------------------------------------------- #

# -------------------------------------------------------------------- #
def read_binary_state(state_fname):
    ''' Read in ascii format state file and convert to a list of numbers

    Parameters
    ----------
    state_fname: <str>
        Path of the state file to be read

    Returns
    ----------
    states: <bytes>
        The full binary state file content
    '''

    with open(state_fname, 'rb') as f:
        states = f.read()
    return states

# -------------------------------------------------------------------- #

# -------------------------------------------------------------------- #
if __name__ == '__main__':
    main()
# -------------------------------------------------------------------- #
