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
        infile = os.path.join(test_dir, 'samples',
                              test_dict['global_parameter_file'])

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
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # replace global options from config file
        if 'options' in test_dict:
            replacements = test_dict['options']
        else:
            replacements = OrderedDict()
        global_param = replace_global_values(global_param, replacements)
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # write global parameter file
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
            returncode = vic_exe.run(test_global_file, logdir=dirs['logs'])
            test_complete = True

            expectation = test_dict['expected_retval']
            if returncode != int(test_dict['expected_retval']):
                raise VICReturnCodeError('VIC return code ({0}) '
                                         'does not match expected '
                                         '({1})'.format(returncode,
                                                        expectation))

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
if __name__ == '__main__':
    main()
# -------------------------------------------------------------------- #
