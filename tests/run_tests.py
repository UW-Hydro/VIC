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
import warnings

import pytest

from tonic.models.vic.vic import VIC
from tonic.io import read_config, read_configobj
from tonic.testing import VICTestError
from test_utils import (setup_test_dirs, print_test_dict,
                        replace_global_values, drop_tests, pop_run_kwargs,
                        check_returncode, process_error,
                        test_classic_driver_all_complete,
                        test_classic_driver_no_output_file_nans,
                        find_global_param_value,
                        check_multistream,
                        plot_science_tests)
from test_image_driver import test_image_driver_no_output_file_nans,
                        check_multistream_classic)
from test_image_driver import (test_image_driver_no_output_file_nans,
                               check_multistream_image)
from test_restart import (prepare_restart_run_periods,
                          setup_subdirs_and_fill_in_global_param_restart_test,
                          check_exact_restart_fluxes,
                          check_exact_restart_states)

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


def main():
    '''
    Run VIC tests
    '''

    # dates and times
    starttime = datetime.datetime.now()
    ymd = starttime.strftime('%Y%m%d')

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
                        default=os.path.join(test_dir, 'system/system_tests.cfg'))
    parser.add_argument('--science', type=str,
                        help='science tests configuration file',
                        default=os.path.join(test_dir, 'science/science.cfg'))
    parser.add_argument('--examples', type=str,
                        help='examples tests configuration file',
                        default=os.path.join(test_dir, 'examples/examples.cfg'))
    parser.add_argument('--release', type=str,
                        help='release tests configuration file',
                        default=os.path.join(test_dir, 'release/release.cfg'))
    parser.add_argument('--vic_exe', type=str,
                        help='VIC executable to test',
                        default=os.path.join(
                            test_dir, '../vic/drivers/classic/vic_classic.exe'))
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
    parser.add_argument('--science_test_data_dir', type=str,
                        help='directory to find science test data',
                        default='./test_runs')
                        default=os.path.join(test_dir, '../samples/VIC_sample_data'))

    args = parser.parse_args()

    # Define test directories
    data_dir = args.data_dir
    out_dir = os.path.expandvars(args.output_dir)
    os.makedirs(out_dir, exist_ok=True)

    # check to make sure science test data directory exists
    science_test_data_dir = args.science_test_data_dir
    if not os.path.exists(science_test_data_dir):
        raise VICTestError("directory for science test data does not exist or has not been defined")

    # Validate input directories
    if not (len(args.tests) == 1 and args.tests[0] == 'unit'):
        for d in [data_dir, test_dir]:
            if not os.path.exists(d):
                raise VICTestError('Directory: {0} does not exist'.format(d))

    # Print welcome information
    print(description)
    print('\nStarting tests now...Start Time: {0}\n'.format(starttime))
    print('Running Test Set: {0}'.format(', '.join(args.tests)))

    # Setup VIC executable
    if not (len(args.tests) == 1 and args.tests[0] == 'unit'):
        vic_exe = VIC(args.vic_exe)
        print('VIC version information:\n\n{0}'.format(vic_exe.version.decode()))

    # run test sets
    # unit
    if any(i in ['all', 'unit'] for i in args.tests):
        test_results['unit'] = run_unit_tests(test_dir)

    # system
    if any(i in ['all', 'system'] for i in args.tests):
        test_results['system'] = run_system(args.system, vic_exe, data_dir,
                                            os.path.join(out_dir, 'system'),
                                            args.driver)
    # science
    if any(i in ['all', 'science'] for i in args.tests):
        test_results['science'] = run_science(args.science, vic_exe,
                                              science_test_data_dir,
                                              data_dir,
                                              os.path.join(out_dir, 'science'),
                                              args.driver)
    # examples
    if any(i in ['all', 'examples'] for i in args.tests):
        test_results['examples'] = run_examples(args.examples, vic_exe,
                                                data_dir,
                                                os.path.join(
                                                    out_dir, 'examples'),
                                                args.driver)
    # release
    if any(i in ['all', 'release'] for i in args.tests):
        test_results['release'] = run_release(args.release)

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

    # end date and times
    endtime = datetime.datetime.now()
    elapsed = endtime - starttime
    print('\nFinished testing VIC. Endtime: {0}'.format(endtime))
    print('Time elapsed during testing:  {0}\n'.format(elapsed))

    # return exit code
    sys.exit(failed)


def run_unit_tests(test_dir):
    '''Run unittests in test_dir

    Parameters
    ----------
    test_dir : str
        Path to unittests

    Returns
    -------
    test_results : dict
        Test results for all tests in config_file.

    See Also
    --------
    run_system
    run_examples
    run_science
    run_release

    '''

    print('\n-'.ljust(OUTPUT_WIDTH + 1, '-'))
    print('Running Unit Tests')
    print('-'.ljust(OUTPUT_WIDTH, '-'))

    retcode = pytest.main(['-x',  os.path.join(test_dir, 'unit'), '--boxed'])
    return {'unittests': TestResults('unittests',
                                     test_complete=True,
                                     passed=retcode == 0,
                                     comment='see stdout from pytest',
                                     returncode=retcode)}

    print('\n-'.ljust(OUTPUT_WIDTH + 1, '-'))
    print('Finished unit tests.')
    print('-'.ljust(OUTPUT_WIDTH, '-'))


def run_system(config_file, vic_exe, test_data_dir, out_dir, driver):
    '''Run system tests from config file

    Parameters
    ----------
    config_file : str
        Configuration file for system tests.
    vic_exe : VIC (object)
        VIC executable object (see tonic documentation).
    test_data_dir : str
        Path to test data sets.
    out_dir : str
        Path to output location
    driver : {'classic', 'image'}
        Driver to run tests on.

    Returns
    -------
    test_results : dict
        Test results for all tests in config_file.

    See Also
    --------
    run_unit_tests
    run_examples
    run_science
    run_release
    '''

    # Print test set welcome
    print('\n-'.ljust(OUTPUT_WIDTH + 1, '-'))
    print('Running System Tests')
    print('-'.ljust(OUTPUT_WIDTH, '-'))

    # Get setup
    config = read_configobj(config_file)

    # drop invalid driver tests
    config = drop_tests(config, driver)

    test_results = OrderedDict()

    # Run individual system tests
    for i, (testname, test_dict) in enumerate(config.items()):

        # print out status info
        print('Running test {0}/{1}: {2}'.format(i + 1, len(config.items()),
                                                 testname))

        # Setup directories for test
        dirs = setup_test_dirs(testname, out_dir,
                               mkdirs=['results', 'state', 'logs', 'plots'])

        # read template global parameter file
        infile = os.path.join(test_dir, 'system',
                              test_dict['global_parameter_file'])

        with open(infile, 'r') as global_file:
            global_param = global_file.read()

        # If restart test, prepare running periods
        # (1) Find STATESEC option (and STATE_FORMAT option for later use)
        statesec = find_global_param_value(global_param, 'STATESEC')
        if driver == 'classic':
            state_format = find_global_param_value(global_param,
                                                   'STATE_FORMAT')
        # (2) Prepare running periods and initial state file info for restart
        # test
        if 'exact_restart' in test_dict['check']:
            run_periods = prepare_restart_run_periods(
                                test_dict['restart'],
                                dirs['state'], statesec)

        # create template string
        s = string.Template(global_param)

        # fill in global parameter options
        # --- if restart test, multiple runs --- #
        if 'exact_restart' in test_dict['check']:
            # Set up subdirectories and fill in global parameter options
            # for restart testing
            list_global_param =\
                setup_subdirs_and_fill_in_global_param_restart_test(
                    s, run_periods, driver, dirs['results'], dirs['state'],
                    test_data_dir)
        # else, single run
        else:
            global_param = s.safe_substitute(test_data_dir=test_data_dir,
                                             result_dir=dirs['results'],
                                             state_dir=dirs['state'])

        # replace global options from config file
        # --- extract global options to be substitute --- #
        if 'options' in test_dict:
            replacements = test_dict['options']
        else:
            replacements = OrderedDict()
        # --- if STATE_FORMAT is specified, then the specified value (instead
        # of the one in the global template file) --- #
        if 'STATE_FORMAT' in replacements:
            state_format = replacements['STATE_FORMAT']
        # --- replace global options --- #
        if 'exact_restart' in test_dict['check']:
            for j, gp in enumerate(list_global_param):
                # save a copy of replacements for the next global file
                replacements_cp = replacements.copy()
                # replace global options for this global file
                list_global_param[j] = replace_global_values(gp, replacements)
                replacements = replacements_cp
        else:
            global_param = replace_global_values(global_param, replacements)

        # write global parameter file
        if 'exact_restart' in test_dict['check']:
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
            test_global_file = os.path.join(
                                dirs['test'],
                                '{0}_globalparam.txt'.format(testname))
            with open(test_global_file, mode='w') as f:
                for line in global_param:
                    f.write(line)

        # Get optional kwargs for run executable
        run_kwargs = pop_run_kwargs(test_dict)

        # run VIC
        test_complete = False
        test_passed = False
        test_comment = ''
        error_message = ''

        try:
            if 'exact_restart' in test_dict['check']:
                for j, test_global_file in enumerate(list_test_global_file):
                    returncode = vic_exe.run(test_global_file,
                                             logdir=dirs['logs'])
                    # Check return code
                    check_returncode(vic_exe,
                                     test_dict.pop('expected_retval', 0))
            else:
                returncode = vic_exe.run(test_global_file, logdir=dirs['logs'],
                                         **run_kwargs)
                # Check return code
                check_returncode(vic_exe,
                                 test_dict.pop('expected_retval', 0))

            test_complete = True

            # check output files (different tests depending on driver)
            if 'check' in test_dict:

                # Check that the simulation completed for all grid cells
                if 'complete' in test_dict['check']:
                    fnames = glob.glob(os.path.join(dirs['results'], '*'))
                    if driver == 'classic':
                        test_classic_driver_all_complete(fnames)

                # check for nans in all example files
                if 'output_file_nans' in test_dict['check']:
                    fnames = glob.glob(os.path.join(dirs['results'], '*'))
                    if driver == 'classic':
                        test_classic_driver_no_output_file_nans(fnames)
                    elif driver == 'image':
                        domain_file = os.path.join(test_data_dir,
                                                   test_dict['domain_file'])
                        test_image_driver_no_output_file_nans(fnames,
                                                              domain_file)
                    else:
                        raise ValueError('unknown driver')

                # check for exact restarts
                if 'exact_restart' in test_dict['check']:
                    check_exact_restart_fluxes(dirs['results'], driver,
                                               run_periods)
                    if driver == 'classic':
                        check_exact_restart_states(dirs['state'], driver,
                                                   run_periods, statesec,
                                                   state_format)
                    elif driver == 'image':
                        check_exact_restart_states(dirs['state'], driver,
                                                   run_periods, statesec)

                if 'multistream' in test_dict['check']:
                    fnames = glob.glob(os.path.join(dirs['results'], '*'))
                    if driver == 'classic':
                        check_multistream_classic(fnames)
                    elif driver == 'image':
                        warnings.warn('Skipping multistream image driver test')
                        # TODO: check_multistream_image(fnames)

            # if we got this far, the test passed.
            test_passed = True

        # Handle errors
        except Exception as e:
            test_comment, error_message = process_error(e, vic_exe)

        # record the test results
        test_results[testname] = TestResults(testname,
                                             test_complete=test_complete,
                                             passed=test_passed,
                                             comment=test_comment,
                                             error_message=error_message,
                                             returncode=returncode)

    print('-'.ljust(OUTPUT_WIDTH, '-'))
    print('Finished testing system tests.')
    print('-'.ljust(OUTPUT_WIDTH, '-'))

    return test_results


def run_science(config_file, vic_exe, science_test_data_dir,
                test_data_dir, out_dir, driver):
    '''Run science tests from config file

    Parameters
    ----------
    config_file : str
        Configuration file for science tests.
    vic_exe : VIC (object)
        VIC executable object (see tonic documentation).
    science_test_data_dir: str
        Path to science test data sets (archived VIC runs and observations)
    test_data_dir : str
        Path to test data sets.
    out_dir : str
        Path to output location
    driver : {'classic', 'image'}
        Driver to run tests on.

    Returns
    -------
    test_results : dict
        Test results for all tests in config_file.

    See Also
    --------
    run_unit_tests
    run_examples
    run_system
    run_release
    '''

    # Print test set welcome
    print('\n-'.ljust(OUTPUT_WIDTH + 1, '-'))
    print('Running Science Tests')
    print('-'.ljust(OUTPUT_WIDTH, '-'))

    # Get setup
    config = read_config(config_file)

    # drop invalid driver tests
    config = drop_tests(config, driver)

    test_results = OrderedDict()

    # Run individual tests
    for i, (testname, test_dict) in enumerate(config.items()):

        # print out status info
        print('Running test {0}/{1}: {2}'.format(i + 1, len(config.items()),
                                                 testname))

        # Setup directories for test
        dirs = setup_test_dirs(testname, out_dir,
                               mkdirs=['results', 'state', 'logs', 'plots'])

        # read template global parameter file
        infile = os.path.join(test_dir, 'science', test_dict['global_parameter_file'])
        print(infile)

        with open(infile, 'r') as global_file:
            global_param = global_file.read()

        # create template string
        s = string.Template(global_param)

        # fill in global parameter options
        global_param = s.safe_substitute(test_data_dir=science_test_data_dir,
                                         result_dir=dirs['results'],
                                         state_dir=dirs['state'],
                                         testname=testname,
                                         test_root=test_dir)

        test_global_file = os.path.join(dirs['test'],
                                        '{0}_globalparam.txt'.format(testname))

        # write global parameter file
        with open(test_global_file, 'w') as f:
            f.write(global_param)

        # Get optional kwargs for run executable
        run_kwargs = pop_run_kwargs(test_dict)

        # run VIC
        test_complete = False
        test_passed = False
        test_comment = ''
        error_message = ''

        try:
            # Run the VIC simulation
            returncode = vic_exe.run(test_global_file, logdir=dirs['logs'],
                                     **run_kwargs)
            test_complete = True

            # Check return code
            check_returncode(vic_exe)

            # check output files (different tests depending on driver)
            if test_dict['check']:
                fnames = glob.glob(os.path.join(dirs['results'], '*'))

                # Check that the simulation completed for all grid cells
                if 'complete' in test_dict['check'] and driver == 'classic':
                    test_classic_driver_all_complete(fnames)

                # check for nans in all example files
                if 'output_file_nans' in test_dict['check']:
                    if driver == 'classic':
                        test_classic_driver_no_output_file_nans(fnames)
                    elif driver == 'image':
                        domain_file = os.path.join(test_data_dir,
                                                   test_dict['domain_file'])
                        test_image_driver_no_output_file_nans(fnames,
                                                              domain_file)
                    else:
                        raise ValueError('unknown driver')

            # if we got this far, the test passed.
            test_passed = True

            # plot science test results
            plot_science_tests(test_dict['driver'],
                                testname,
                                science_test_data_dir,
                                dirs['results'],
                                dirs['plots'],
                                test_dict['vic4.2.d'],
                                test_dict['vic5.0.0'],
                                test_dict['observations_path'],
                                test_dict['plots'])

        # Handle errors
        except Exception as e:
            test_comment, error_message = process_error(e, vic_exe)

        # record the test results
        test_results[testname] = TestResults(testname,
                                             test_complete=test_complete,
                                             passed=test_passed,
                                             comment=test_comment,
                                             error_message=error_message,
                                             returncode=returncode)

    print('-'.ljust(OUTPUT_WIDTH, '-'))
    print('Finished testing science tests.')
    print('-'.ljust(OUTPUT_WIDTH, '-'))

    return test_results


def run_examples(config_file, vic_exe, test_data_dir, out_dir, driver):
    '''Run examples tests from config file

    Parameters
    ----------
    config_file : str
        Configuration file for example tests.
    vic_exe : VIC (object)
        VIC executable object (see tonic documentation).
    test_data_dir : str
        Path to test data sets.
    out_dir : str
        Path to output location
    driver : {'classic', 'image'}
        Driver to run tests on.

    Returns
    -------
    test_results : dict
        Test results for all tests in config_file.

    See Also
    --------
    run_unit_tests
    run_system
    run_science
    run_release
    '''

    # Print test set welcome
    print('\n-'.ljust(OUTPUT_WIDTH + 1, '-'))
    print('Running Examples')
    print('-'.ljust(OUTPUT_WIDTH, '-'))

    # Get setup
    config = read_config(config_file)

    # drop invalid driver tests
    config = drop_tests(config, driver)

    test_results = OrderedDict()

    # Run individual examples
    for i, (testname, test_dict) in enumerate(config.items()):

        # print out status info
        print('Running test {0}/{1}: {2}'.format(i + 1, len(config.items()),
                                                 testname))

        # Setup directories for test
        dirs = setup_test_dirs(testname, out_dir,
                               mkdirs=['results', 'state', 'logs', 'plots'])

        # read template global parameter file
        infile = os.path.join(test_dir, 'examples',
                              test_dict['global_parameter_file'])

        with open(infile, 'r') as global_file:
            global_param = global_file.read()

        # create template string
        s = string.Template(global_param)

        # fill in global parameter options
        global_param = s.safe_substitute(test_data_dir=test_data_dir,
                                         result_dir=dirs['results'],
                                         state_dir=dirs['state'],
                                         testname=testname,
                                         test_root=test_dir)

        test_global_file = os.path.join(dirs['test'],
                                        '{0}_globalparam.txt'.format(testname))

        # write global parameter file
        with open(test_global_file, 'w') as f:
            f.write(global_param)

        # Get optional kwargs for run executable
        run_kwargs = pop_run_kwargs(test_dict)

        # run VIC
        test_complete = False
        test_passed = False
        test_comment = ''
        error_message = ''

        try:
            # Run the VIC simulation
            returncode = vic_exe.run(test_global_file, logdir=dirs['logs'],
                                     **run_kwargs)
            test_complete = True

            # Check return code
            check_returncode(vic_exe)

            # check output files (different tests depending on driver)
            if test_dict['check']:
                fnames = glob.glob(os.path.join(dirs['results'], '*'))

                # Check that the simulation completed for all grid cells
                if 'complete' in test_dict['check'] and driver == 'classic':
                    test_classic_driver_all_complete(fnames)

                # check for nans in all example files
                if 'output_file_nans' in test_dict['check']:
                    if driver == 'classic':
                        test_classic_driver_no_output_file_nans(fnames)
                    elif driver == 'image':
                        domain_file = os.path.join(test_data_dir,
                                                   test_dict['domain_file'])
                        test_image_driver_no_output_file_nans(fnames,
                                                              domain_file)
                    else:
                        raise ValueError('unknown driver')

            # if we got this far, the test passed.
            test_passed = True

        # Handle errors
        except Exception as e:
            test_comment, error_message = process_error(e, vic_exe)

        # record the test results
        test_results[testname] = TestResults(testname,
                                             test_complete=test_complete,
                                             passed=test_passed,
                                             comment=test_comment,
                                             error_message=error_message,
                                             returncode=returncode)

    # Print examples footer
    print('-'.ljust(OUTPUT_WIDTH, '-'))
    print('Finished testing examples.')
    print('-'.ljust(OUTPUT_WIDTH, '-'))

    return test_results


def run_release(config_file):
    '''Run release from config file

    NOT IMPLEMENTED
    '''
    return OrderedDict()


if __name__ == '__main__':
    main()
