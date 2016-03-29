import os
from sys import platform
import subprocess
import pytest


class CESMTestException(Exception):
    pass

if os.getenv('TESTID', 'notset') != 'cesm':
    wrong_driver = True
else:
    wrong_driver = False
    test_dir = os.path.dirname(__file__)
    shared_object = os.path.join(test_dir, os.pardir, os.pardir,
                                 'vic', 'drivers', 'cesm', 'lndlib.a')
    if os.path.isfile(shared_object):
        lib_built = True
    else:
        raise CESMTestException('lndlib.a has not been built!')


def nm(shared_object, *args):
    '''wrapper of the unix nm utility'''

    if not os.path.isfile(shared_object):
        raise FileNotFoundError('%s is not a file' % shared_object)

    cmd = ['nm', '-g']
    cmd.extend(args)
    cmd.append(shared_object)
    output = subprocess.check_output(cmd).decode()

    # parse the output
    objs = [line.split(' ')[-1] for line in output.split('\n')]

    # the nm utility on osx includes a prepended underscore on all the objects
    if platform == 'darwin':
        for i, obj in enumerate(objs):
            if len(obj) > 1 and obj[0] == '_':
                objs[i] = obj[1:]
    return objs


@pytest.mark.skipif(wrong_driver, reason='Only run for CESM driver builds')
def test_vic_lib_present():
    assert lib_built


@pytest.mark.skipif(wrong_driver, reason='Only run for CESM driver builds')
def test_fortran_symbols_present():
    lib_symbols = nm(shared_object)

    vic_fortran_symbols = ['initialize_log',
                           'initialize_vic_cesm_mpi',
                           'vic_cesm_init',
                           'vic_cesm_run',
                           'vic_cesm_final']

    for symbol in vic_fortran_symbols:
        assert symbol in lib_symbols
