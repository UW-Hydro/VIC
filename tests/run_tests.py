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

from tonic.models.vic.vic import VIC, default_vic_valgrind_suppressions_path
from tonic.io import read_config, read_configobj
from tonic.testing import VICTestError
from test_utils import (
    setup_test_dirs, print_test_dict,
    replace_global_values, drop_tests, pop_run_kwargs,
    check_returncode, process_error,
    test_classic_driver_all_complete,
    test_classic_driver_no_output_file_nans,
    find_global_param_value,
    check_multistream_classic,
    setup_subdirs_and_fill_in_global_param_driver_match_test,
    check_drivers_match_fluxes,
    plot_science_tests)
from test_image_driver import (test_image_driver_no_output_file_nans,
                               setup_subdirs_and_fill_in_global_param_mpi_test,
                               check_mpi_fluxes, check_mpi_states)
from test_restart import (prepare_restart_run_periods,
                          setup_subdirs_and_fill_in_global_param_restart_test,
                          check_exact_restart_fluxes,
                          check_exact_restart_states)

test_dir = os.path.dirname(os.path.abspath(__file__))

# Set path to valgrind supressions file if not already set.
if 'VIC_VALGRIND_SUPPRESSIONS' not in os.environ:
    sup_file = os.path.join(test_dir, default_vic_valgrind_suppressions_path)
    if os.path.isfile(sup_file):
        os.environ["VIC_VALGRIND_SUPPRESSIONS"] = sup_file

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
                        default=os.path.join(test_dir,
                                             'system/system_tests.cfg'))
    parser.add_argument('--science', type=str,
                        help='science tests configuration file',
                        default=os.path.join(test_dir, 'science/science.cfg'))
    parser.add_argument('--examples', type=str,
                        help='examples tests configuration file',
                        default=os.path.join(test_dir,
                                             'examples/examples.cfg'))
    parser.add_argument('--release', type=str,
                        help='release tests configuration file',
                        default=os.path.join(test_dir, 'release/release.cfg'))
    parser.add_argument('--classic', type=str,
                        help='classic driver executable to test')
    parser.add_argument('--image', type=str,
                        help='image driver executable to test')
    parser.add_argument('--output_dir', type=str,
                        help='directory to write test output to',
                        default='$WORKDIR/VIC_tests_{0}'.format(ymd))
    parser.add_argument('--data_dir', type=str,
                        help='directory to find test data',
                        default='./samples/VIC_sample_data')
    parser.add_argument('--science_test_data_dir', type=str,
                        help='directory to find science test data',
                        default='./samples/VIC_sample_data')
    parser.add_argument('--nproc', type=int,
                        help='number of processors to use for science tests',
                        default=1)

    args = parser.parse_args()

    # Define test directories
    data_dir = args.data_dir
    out_dir = os.path.expandvars(args.output_dir)
    os.makedirs(out_dir, exist_ok=True)

    # check to make sure science test data directory exists
    science_test_data_dir = args.science_test_data_dir
    if 'science' in args.tests and not os.path.exists(science_test_data_dir):
        raise VICTestError('directory for science test data does not exist or '
                           'has not been defined')

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
    # --- if not only unit test --- #
    if not (len(args.tests) == 1 and args.tests[0] == 'unit'):
        dict_drivers = {}
        if args.classic:
            dict_drivers['classic'] = VIC(args.classic)
            print('VIC classic version information:\n\n{0}'.format(
                dict_drivers['classic'].version.decode()))
        if args.image:
            dict_drivers['image'] = VIC(args.image)
            print('VIC image version information:\n\n{0}'.format(
                dict_drivers['image'].version.decode()))

    # run test sets
    # unit
    if any(i in ['all', 'unit'] for i in args.tests):
        test_results['unit'] = run_unit_tests(test_dir)

    # system
    if any(i in ['all', 'system'] for i in args.tests):
        test_results['system'] = run_system(args.system, dict_drivers,
                                            data_dir,
                                            os.path.join(out_dir,
                                                         'system'))

    # science
    if any(i in ['all', 'science'] for i in args.tests):
        test_results['science'] = run_science(
            args.science, dict_drivers['classic'],
            science_test_data_dir,
            data_dir,
            os.path.join(out_dir, 'science'),
            'classic',
            args.nproc)
    # examples
    if any(i in ['all', 'examples'] for i in args.tests):
        if len(dict_drivers) == 1:  # if only one driver
            driver = list(dict_drivers.keys())[0]
            vic_exe = dict_drivers[driver]
            test_results['examples'] = run_examples(args.examples, vic_exe,
                                                    data_dir,
                                                    os.path.join(
                                                        out_dir, 'examples'),
                                                    driver)
        else:
            raise ValueError('example test only supports single driver')
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

    retcode = pytest.main(['-x', os.path.join(test_dir, 'unit'), '--boxed'])
    return {'unittests': TestResults('unittests',
                                     test_complete=True,
                                     passed=retcode == 0,
                                     comment='see stdout from pytest',
                                     returncode=retcode)}

    print('\n-'.ljust(OUTPUT_WIDTH + 1, '-'))
    print('Finished unit tests.')
    print('-'.ljust(OUTPUT_WIDTH, '-'))


def run_system(config_file, dict_drivers, test_data_dir, out_dir):
    '''Run system tests from config file

    Parameters
    ----------
    config_file : str
        Configuration file for system tests.
    dict_drivers : dict
        Keys: driver names {'classic', 'image'}
        Content: corresponding VIC executable object (see tonic documentation)
    test_data_dir : str
        Path to test data sets.
    out_dir : str
        Path to output location

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

    # Process driver info
    if len(dict_drivers) == 1:  # if single driver
        driver = list(dict_drivers.keys())[0]
        vic_exe = dict_drivers[driver]

    # Drop invalid driver tests
    if len(dict_drivers) == 1:  # if single driver
        config = drop_tests(config, driver)
    else:  # if multiple drivers
        config = drop_tests(config, list(dict_drivers))

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
        dict_global_param = {}
        # --- if single driver --- #
        if len(dict_drivers) == 1:
            infile = os.path.join(test_dir, 'system',
                                  test_dict['global_parameter_file'])
            with open(infile, 'r') as global_file:
                dict_global_param[driver] = global_file.read()
        # --- if multiple drivers --- #
        else:
            for j, dr in enumerate(test_dict['driver']):
                infile = os.path.join(test_dir, 'system',
                                      test_dict['global_parameter_file'][j])
                with open(infile, 'r') as global_file:
                    dict_global_param[dr] = global_file.read()

        # If restart test, prepare running periods
        if 'exact_restart' in test_dict['check']:
            if len(dict_drivers) > 1:
                raise ValueError('Only support single driver for restart'
                                 'tests!')
            global_param = dict_global_param[driver]
            # () Find STATE_FORMAT option for later use
            if driver == 'classic':
                state_format = find_global_param_value(global_param,
                                                       'STATE_FORMAT')
            # (2) Prepare running periods and initial state file info for
            # restart test
            run_periods = prepare_restart_run_periods(
                test_dict['restart'],
                dirs['state'])

        # If mpi test, prepare a list of number of processors to be run
        elif 'mpi' in test_dict['check']:
            if len(dict_drivers) > 1:
                raise ValueError('Only support single driver for MPI'
                                 'tests!')
            if not isinstance(test_dict['mpi']['n_proc'], list):
                raise ValueError('Need at least two values in n_proc to run'
                                 'mpi test!')
            list_n_proc = test_dict['mpi']['n_proc']

        # create template string
        dict_s = {}
        for dr, global_param in dict_global_param.items():
            dict_s[dr] = string.Template(global_param)

        # fill in global parameter options
        # --- if restart test, multiple runs --- #
        if 'exact_restart' in test_dict['check']:
            s = dict_s[driver]
            # Set up subdirectories and fill in global parameter options
            # for restart testing
            list_global_param =\
                setup_subdirs_and_fill_in_global_param_restart_test(
                    s, run_periods, driver, dirs['results'], dirs['state'],
                    test_data_dir)
        # --- if mpi test, multiple runs --- #
        elif 'mpi' in test_dict['check']:
            s = dict_s[driver]
            # Set up subdirectories and output directories in global file for
            # multiprocessor testing
            list_global_param = \
                setup_subdirs_and_fill_in_global_param_mpi_test(
                    s, list_n_proc, dirs['results'], dirs['state'],
                    test_data_dir)
        # --- if driver-match test, one run for each driver --- #
        elif 'driver_match' in test_dict['check']:
            # Set up subdirectories and output directories in global file for
            # driver-match testing
            dict_global_param = \
                setup_subdirs_and_fill_in_global_param_driver_match_test(
                    dict_s, dirs['results'], dirs['state'], test_data_dir)
        # --- else, single run --- #
        else:
            if len(dict_drivers) > 1:
                raise RuntimeError('Only support single driver for test'
                                   '{}!'.format(testname))
            s = dict_s[driver]
            global_param = s.safe_substitute(test_data_dir=test_data_dir,
                                             result_dir=dirs['results'],
                                             state_dir=dirs['state'])

        # replace global options from config file
        # --- extract global options to be substitute --- #
        if 'options' in test_dict:
            replacements = test_dict['options']
        else:
            replacements = OrderedDict()
        # --- replace global options --- #
        # For the purpose of exact restart, if STATE_FORMAT is specified,
        # then record the specified value (instead of the one in the global
        # template file)
        if 'exact_restart' in test_dict['check']:
            if 'STATE_FORMAT' in replacements:
                state_format = replacements['STATE_FORMAT']
        if 'exact_restart' in test_dict['check'] or\
           'mpi' in test_dict['check']:  # if multiple runs
            for j, gp in enumerate(list_global_param):
                # save a copy of replacements for the next global file
                replacements_cp = replacements.copy()
                # replace global options for this global file
                list_global_param[j] = replace_global_values(gp, replacements)
                replacements = replacements_cp
        elif 'driver_match' in test_dict['check']:  # if cross-driver runs
            for dr, gp in dict_global_param.items():
                # save a copy of replacements for the next global file
                replacements_cp = replacements.copy()
                # replace global options for this global file
                dict_global_param[dr] = replace_global_values(gp,
                                                              replacements)
                replacements = replacements_cp
        else:  # if single run
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
        elif 'mpi' in test_dict['check']:
            list_test_global_file = []
            for j, gp in enumerate(list_global_param):
                test_global_file = os.path.join(
                    dirs['test'],
                    '{}_globalparam_processors_{}.txt'.format(
                        testname, list_n_proc[j]))
                list_test_global_file.append(test_global_file)
                with open(test_global_file, mode='w') as f:
                    for line in gp:
                        f.write(line)
        elif 'driver_match' in test_dict['check']:
            dict_test_global_file = {}
            for dr, gp in dict_global_param.items():
                test_global_file = os.path.join(
                    dirs['test'],
                    '{}_globalparam_{}.txt'.format(
                        testname, dr))
                dict_test_global_file[dr] = test_global_file
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
                                             logdir=dirs['logs'],
                                             **run_kwargs)
                    # Check return code
                    check_returncode(vic_exe,
                                     test_dict.pop('expected_retval', 0))
            elif 'mpi' in test_dict['check']:
                for j, test_global_file in enumerate(list_test_global_file):
                    # Overwrite mpi_proc in option kwargs
                    n_proc = list_n_proc[j]
                    if n_proc == 1:
                        run_kwargs['mpi_proc'] = None
                    else:
                        run_kwargs['mpi_proc'] = list_n_proc[j]
                    # Run VIC
                    returncode = vic_exe.run(test_global_file,
                                             logdir=dirs['logs'],
                                             **run_kwargs)
                    # Check return code
                    check_returncode(vic_exe,
                                     test_dict.pop('expected_retval', 0))
            elif 'driver_match' in test_dict['check']:
                for dr in dict_test_global_file.keys():
                    # Reset mpi_proc in option kwargs to None for classic
                    # driver run
                    if dr == 'classic':
                        run_kwargs_classic = run_kwargs
                        run_kwargs_classic['mpi_proc'] = None
                        returncode = dict_drivers[dr].run(
                            dict_test_global_file[dr],
                            logdir=dirs['logs'],
                            **run_kwargs_classic)
                    else:
                        returncode = dict_drivers[dr].run(
                            dict_test_global_file[dr],
                            logdir=dirs['logs'],
                            **run_kwargs)
                    # Check return code
                    check_returncode(dict_drivers[dr],
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
                    if len(dict_drivers) > 1:
                        raise RuntimeError('Only support single driver for '
                                           'complete check')
                    fnames = glob.glob(os.path.join(dirs['results'], '*'))
                    if driver == 'classic':
                        test_classic_driver_all_complete(fnames)
                    else:
                        raise RuntimeError('complete check only supports '
                                           'classic driver')

                # check for nans in all example files
                if 'output_file_nans' in test_dict['check']:
                    if len(dict_drivers) > 1:
                        raise RuntimeError('Only support single driver for '
                                           'output_file_nans check')
                    fnames = glob.glob(os.path.join(dirs['results'], '*'))
                    if driver == 'classic':
                        test_classic_driver_no_output_file_nans(fnames)
                    elif driver == 'image':
                        domain_file = os.path.join(
                            test_data_dir,
                            test_dict['domain_file'])
                        test_image_driver_no_output_file_nans(
                            fnames,
                            domain_file)
                    else:
                        raise ValueError('unknown driver')

                # check for exact restarts
                if 'exact_restart' in test_dict['check']:
                    check_exact_restart_fluxes(dirs['results'], driver,
                                               run_periods)
                    if driver == 'classic':
                        check_exact_restart_states(dirs['state'], driver,
                                                   run_periods,
                                                   state_format)
                    elif driver == 'image':
                        check_exact_restart_states(dirs['state'], driver,
                                                   run_periods)
                    else:
                        raise ValueError('unknown driver')

                # check for multistream output
                if 'multistream' in test_dict['check']:
                    if len(dict_drivers) > 1:
                        raise ValueError('Only support single driver for '
                                         'multistream check')
                    fnames = glob.glob(os.path.join(dirs['results'], '*'))
                    if driver == 'classic':
                        check_multistream_classic(fnames)
                    elif driver == 'image':
                        warnings.warn('Skipping multistream image driver test')
                        # TODO: check_multistream_image(fnames)

                # check for mpi multiprocessor results
                if 'mpi' in test_dict['check']:
                    check_mpi_fluxes(dirs['results'], list_n_proc)
                    check_mpi_states(dirs['state'], list_n_proc)

                # check that results from different drivers match
                if 'driver_match' in test_dict['check']:
                    check_drivers_match_fluxes(list(dict_drivers.keys()),
                                               dirs['results'])

            # if we got this far, the test passed.
            test_passed = True

        # Handle errors
        except Exception as e:
            for dr, exe in dict_drivers.items():
                test_comment, error_message = process_error(e, exe)

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
                test_data_dir, out_dir, driver, nproc):
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
    nproc : int
        Number of processors to use for science tests

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
    config = read_configobj(config_file)

    # drop invalid driver tests
    config = drop_tests(config, driver)

    test_results = OrderedDict()

    # Run individual tests
    for i, (test_type, test_dict) in enumerate(config.items()):

        # print out status info
        print('Running test {0}/{1}: {2}'.format(i + 1, len(config.items()),
                                                 test_type))

        # Setup directories for test
        dirs = setup_test_dirs(test_type, out_dir,
                               mkdirs=['results', 'state', 'logs', 'plots'])

        # read template global parameter file
        infile = os.path.join(test_dir, 'science',
                              test_dict['global_parameter_file'])

        with open(infile, 'r') as global_file:
            global_param = global_file.read()

        # create template string
        s = string.Template(global_param)

        # fill in global parameter options
        global_param = s.safe_substitute(test_data_dir=science_test_data_dir,
                                         test_dir=test_dir,
                                         result_dir=dirs['results'],
                                         state_dir=dirs['state'],
                                         testname=test_type,
                                         test_root=test_dir)

        test_global_file = os.path.join(
            dirs['test'], '{0}_globalparam.txt'.format(test_type))

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

            # plot science test results
            plot_science_tests(test_dict['driver'],
                               test_type,
                               science_test_data_dir,
                               dirs['results'],
                               dirs['plots'],
                               test_dict['plots'],
                               test_dict['compare_data'],
                               nproc=nproc)

            # if we got this far, the test passed.
            test_passed = True

        # Handle errors
        except Exception as e:
            test_comment, error_message = process_error(e, vic_exe)

        # record the test results
        test_results[test_type] = TestResults(test_type,
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
