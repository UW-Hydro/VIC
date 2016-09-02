#!/usr/bin/env python
'''VIC testing command line interface'''

from __future__ import print_function
import os
import argparse
from collections import namedtuple
import psutil
import string
import subprocess
from subprocess import check_call
import datetime
import getpass
import socket
import time

import numpy as np

from tonic.models.vic.vic import VIC

host_config = namedtuple('host_config',
                         ('profile', 'template', 'submit', 'mpiexec'))


def log2_range(m):
    '''
    make an array of integers that increase by 2^n with maximum value of m
    '''
    n = int(np.floor(np.log2(m))) + 1
    return np.exp2(np.arange(n)).astype(np.int)


table_header = '''----------------- START VIC SCALING PROFILE -----------------

Date                      : $date
Machine                   : $hostname
User                      : $user
VIC Test Git Version      : $git_version
VIC Executable            : $vic_exe
VIC Global Parameter File : $vic_global

VIC Executable Version Info
---------------------------
$vic_version
Cores | Time (Seconds)
----------------------
'''


hosts = {
    'local': host_config(profile=[dict(np=np) for np in
                                  log2_range(psutil.cpu_count())],
                         submit=None,
                         template=None,
                         mpiexec=os.getenv('MPIEXEC', 'mpiexec')),
    'hydra': host_config(profile=[dict(np=np) for np in log2_range(64)],
                         submit='qsub', mpiexec='mpiexec',
                         template='''#!/bin/bash
#
#$ -N VIC_scaling_test_$np
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -m be
#$ -pe orte $np

# Qsub template for UW's Hydra Cluster
# Scheduler: SGE
# Valid values for np 1-64
if [ "$np" -gt "64" ]
    echo "$np exceeds maximum number of processes on Hydra"
    exit 1
fi

START=$(date +%s)
$mpiexec -np $np $vic_exe -g $vic_global
END=$(date +%s)
DIFF=$(echo "$END - $START" | bc)
printf "%5s | %f" $np $DIFF >> $timing_table_file'''),
    'topaz': host_config(profile=[dict(select=1, mpiprocs=1),
                                  dict(select=1, mpiprocs=3),
                                  dict(select=1, mpiprocs=9),
                                  dict(select=1, mpiprocs=18),
                                  dict(select=1, mpiprocs=36),
                                  dict(select=2, mpiprocs=36),
                                  dict(select=3, mpiprocs=36),
                                  dict(select=4, mpiprocs=36),
                                  dict(select=5, mpiprocs=36),
                                  dict(select=6, mpiprocs=36),
                                  dict(select=8, mpiprocs=36),
                                  dict(select=10, mpiprocs=36),
                                  dict(select=12, mpiprocs=36)],
                         submit='qsub', mpiexec='mpiexec_mpt',
                         template='''#!/bin/bash

#!/bin/bash
#PBS -N VIC_scaling_test_$i
#PBS -q standard
#PBS -A NPSCA07935242
#PBS -l application=VIC
#PBS -l select=$select:ncpus=36:mpiprocs=$mpiprocs
#PBS -l walltime=06:00:00
#PBS -j oe

# Qsub template for ERDC TOPAZ
# Scheduler: PBS

module load usp-netcdf/intel-15.0.3/4.3.3.1

START=$(date +%s)
mpiexec_mpt -np ${BC_MPI_TASKS_ALLOC} $vic_exe -g $vic_global
END=$(date +%s)
DIFF=$(echo "$END - $START" | bc)
printf "%5s | %f\n" ${BC_MPI_TASKS_ALLOC} $DIFF >> $timing_table_file''')}

OUT_WIDTH = 100

description = '''
                            VIC Test Suite
-------------------------------------------------------------------------------
This is the VIC Profiling Test Suite. There are 2 main test types:

    1. Gprof Profiling: This test will generate a profiling call graph using
        gprof. This test requires building your VIC executable with the
        flags `-pg`.
    2. Scaling: This test will generate a MPI scaling timing table.
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


def main():
    ''' '''
    # dates and times
    starttime = datetime.datetime.now()
    ymd = starttime.strftime('%Y%m%d')

    # parse command line options
    parser = argparse.ArgumentParser(description=description, epilog=epilog,
                                     formatter_class=CustomFormatter)

    parser.add_argument('vic_exe', type=str,
                        help='VIC executable to test')
    parser.add_argument('--kind', type=str,
                        help='Specify which type of test should be run',
                        choices=['scaling', 'profile'],
                        default='scaling')
    parser.add_argument('--host', type=str,
                        help='Host machine to run test on, if not specified, '
                             'test will be run locally',
                        choices=list(hosts.keys()),
                        default='local')
    parser.add_argument('--global_param', '-g', type=str,
                        help='global parameter file to test')
    parser.add_argument('--timing', '-t', type=str,
                        default='vic_timing_{}.txt'.format(ymd),
                        help='path to timing file')
    parser.add_argument('--clean', action='store_true',
                        help='Clean up run files')
    parser.add_argument('--test', action='store_true',
                        help='Test the setup but do not run VIC')

    args = parser.parse_args()

    if args.global_param is None:
        raise ValueError('Global Parameter option is required')

    if args.kind == 'scaling':
        run_scaling(args)
    elif args.kind == 'profile':
        run_profiling(args)
    else:
        raise ValueError('Unknown test kind %s' % args.kind)


def run_profiling(args):
    '''wrapper function for profiling tests'''
    cmd = './profiling/run_gprof.bash -e {vic_exe} -g {vic_global}'.format(
        vic_exe=args.vic_exe, vic_global=args.global_param)
    check_call(cmd, shell=True)


def run_scaling(args):
    '''wrapper function for scaling tests'''
    config = hosts[args.host]
    vic_exe = VIC(args.vic_exe)

    # write timing file header
    header = string.Template(table_header)
    header_kwargs = get_header_info(args.vic_exe, args.global_param)
    header = header.safe_substitute(**header_kwargs)
    with open(args.timing, 'w') as f:
        f.write(header)

    for i, kwargs in enumerate(config.profile):
        if config.template:
            # run on a cluster of some kind
            # start by printing the template
            print('-'.ljust(OUT_WIDTH, '-'))
            print('{host} template'.format(
                host=args.host).center(OUT_WIDTH))
            print('-'.ljust(OUT_WIDTH, '-'))
            print(config.template)
            print('-'.ljust(OUT_WIDTH, '-'))
            template = string.Template(config.template)

            run_string = template.safe_substitute(
                vic_exe=args.vic_exe, vic_global=args.global_param,
                timing_table_file=args.timing, i=i, **kwargs)
            run_file = 'vic_{host}_{i}.sh'.format(host=args.host, i=i)
            with open(run_file, 'w') as f:
                f.write(run_string)

            cmd = '{submit} {run_file}'.format(submit=config.submit,
                                               run_file=run_file)
            print(cmd)
            if not args.test:
                check_call(cmd, shell=True)

            if args.clean:
                os.remove(run_file)
        else:
            # run locally
            n = kwargs['np']
            print('Running {} with {} processors'.format(args.vic_exe, n))
            if not args.test:
                start = time.time()
                vic_exe.run(args.global_param, mpi_proc=int(n))
                end = time.time()
                diff = end - start
                with open(args.timing, 'a') as f:
                    f.write('%5s | %.2f\n' % (n, diff))

    print('See %s for scaling table' % args.timing)


def get_header_info(vic_exe, vic_global):
    '''get info for timing table headers'''
    header_kwargs = {}
    header_kwargs['date'] = datetime.datetime.now()
    header_kwargs['hostname'] = socket.gethostname()
    header_kwargs['user'] = getpass.getuser()
    header_kwargs['git_version'] = subprocess.check_output(
        ['git', 'describe', '--abbrev=4',
         '--dirty', '--always', '--tags']).decode()
    header_kwargs['vic_exe'] = vic_exe
    header_kwargs['vic_global'] = vic_global
    try:
        header_kwargs['vic_version'] = subprocess.check_output(
            [vic_exe, '-v']).decode()
    except subprocess.CalledProcessError:
        pass
    return header_kwargs


if __name__ == '__main__':
    main()
