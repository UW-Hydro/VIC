#!/usr/bin/env python
from __future__ import print_function
import os
import re
import sys
import sysconfig
import warnings
import glob
import subprocess
from tempfile import mkstemp
import shutil

if sys.version_info[0] > 2:
    raise ValueError('ctypes bindings currently unavailable for Python 3')

from setuptools import setup
from setuptools.extension import Extension
from setuptools.command.develop import develop
from setuptools.command.install import install

try:
    import ctypesgencore as ctypesgen
except ImportError:
    ctypesgen = False
    raise
try:
    import autopep8
except:
    autopep8 = False

MAJOR = 5
MINOR = 0
MICRO = 0
ISRELEASED = False
VERSION = '%d.%d.%d' % (MAJOR, MINOR, MICRO)
QUALIFIER = ''

FULLVERSION = VERSION
write_version = False

vic_root_rel_path = os.path.join(os.pardir, os.pardir)
vic_root_abs_path = os.path.abspath(vic_root_rel_path)


# -------------------------------------------------------------------- #
def write_version_py(filename=None):
    cnt = """\
version = '%s'
short_version = '%s'
"""
    if not filename:
        filename = os.path.join(
            os.path.dirname(__file__), 'vic', 'version.py')

    a = open(filename, 'w')
    try:
        a.write(cnt % (FULLVERSION, VERSION))
    finally:
        a.close()

    return
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
def replace(file_path, pattern, subst):
    '''replace lines in file that match pattern with subst'''
    # Create temp file
    fh, abs_path = mkstemp()
    with open(abs_path, 'w') as new_file:
        with open(file_path) as old_file:
            for line in old_file:
                new_file.write(line.replace(pattern, subst))
    os.close(fh)
    # Remove original file
    os.remove(file_path)
    # Move new file
    shutil.move(abs_path, file_path)

    return
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
# Get version string
if ISRELEASED:
    FULLVERSION += QUALIFIER
else:
    FULLVERSION += '.dev'

    pipe = None
    for cmd in ['git', 'git.cmd']:
        try:
            pipe = subprocess.Popen(
                [cmd, "describe", "--always", "--match", "v[0-9]*"],
                stdout=subprocess.PIPE)
            (so, serr) = pipe.communicate()
            if pipe.returncode == 0:
                break
        except:
            pass

    if pipe is None or pipe.returncode != 0:
        # no git, or not in git dir
        if os.path.exists('vic/version.py'):
            warnings.warn("WARNING: Couldn't get git revision, using existing "
                          "vic/version.py")
            write_version = False
        else:
            warnings.warn("WARNING: Couldn't get git revision, using generic "
                          "version string")
    else:
        # have git, in git dir, but may have used a shallow clone
        # (travis does this)
        rev = so.strip()
        # makes distutils blow up on Python 2.7
        if sys.version_info[0] >= 3:
            rev = rev.decode('ascii')

        if not rev.startswith('v') and re.match("[a-zA-Z0-9]{7,9}", rev):
            # partial clone, manually construct version string
            # this is the format before we started using git-describe
            # to get an ordering on dev version strings.
            rev = "v%s.dev-%s" % (VERSION, rev)

        # Strip leading v from tags format "vx.y.z" to get th version string
        FULLVERSION = rev.lstrip('v')
# -------------------------------------------------------------------- #

# -------------------------------------------------------------------- #
if write_version:
    write_version_py()
# -------------------------------------------------------------------- #

sources = []
sources.extend(glob.glob(os.path.join(vic_root_abs_path, 'vic_run', 'src',
                                      '*c')))
sources.extend(glob.glob(os.path.join(vic_root_abs_path, 'drivers', 'shared',
                                      'src', '*c')))
sources.extend(glob.glob(os.path.join(vic_root_abs_path, 'drivers', 'python',
                                      'src', '*c')))

includes = []
includes.append(os.path.join(vic_root_abs_path, 'vic_run', 'include', ''))
includes.append(os.path.join(vic_root_abs_path,
                             'drivers', 'shared', 'include', ''))
includes.append(os.path.join(vic_root_abs_path,
                             'drivers', 'python', 'include', ''))

ext_name = 'vic_core'
# platform safe path to extension
# ext_obj = os.path.join(os.pardir, ext_name + sysconfig.get_config_var('SO'))
ext_obj = ext_name + sysconfig.get_config_var('SO')
ext_module = Extension(ext_name,
                       sources=sources,
                       include_dirs=includes,
                       extra_compile_args=['-std=c99'])


def vic_ctypes_build(command_subclass):
    """A decorator for classes subclassing one of the setuptools commands.

    It modifies the run() method so that it builds the ctypes bindings.
    """
    orig_run = command_subclass.run

    def modified_run(self):
        if ctypesgen:
            ctypes_binding_fname = '_vic_run_lib.py'
            ctypes_binding_fpath = os.path.join('vic', ctypes_binding_fname)
            print('running ctypesgen to generate ctypes bindings to vic_run '
                  'and drivers/shared: {0}'.format(ctypes_binding_fname))

            headers = []
            for inc in includes:
                headers += glob.glob(inc + '*h')

            args = ['ctypesgen.py']
            args += headers
            args += ['-o', ctypes_binding_fpath]
            args += ['-l', ext_obj]

            subprocess.check_output(args)

            # Put location independent path to shared object into ctypes
            # bindings
            replace(ctypes_binding_fpath,
                    '_libs["{0}"] = load_library("{0}")'.format(ext_obj),
                    '_libs["{0}"] = load_library(os.path.join(os.path.dirname('
                    '__file__), os.pardir, "{0}"))'.format(ext_obj))

            if autopep8:
                print('running autopep8 to reformat ctypesgen '
                      'bindings: {0}'.format(ctypes_binding_fname))
                args = ['autopep8', '-i',
                        os.path.join('vic', ctypes_binding_fname),
                        '-a', '-a', '-a']
                subprocess.check_output(args)

            so_name_file = os.path.join('vic', '_vic_run_lib_names.py')
            print('writing private namelist: {0}'.format(so_name_file))
            with open(so_name_file, 'w') as f:
                f.write('"""{0}"""\n'.format(so_name_file))
                f.write('# Autogenerated name list for ctypes object names\n')
                f.write('# This file should not need to be modified\n')
                f.write('vic_core = "{0}"\n'.format(ext_obj))
        else:
            warnings.warn(
                'Using cached version of vic/_vic_run_lib.py ctypes bindings')

        orig_run(self)

    command_subclass.run = modified_run
    return command_subclass


@vic_ctypes_build
class VicDevelopCommand(develop):
    pass


@vic_ctypes_build
class VicInstallCommand(install):
    pass

# -------------------------------------------------------------------- #
# Run Setup
setup(name='vic',
      cmdclass={'install': VicInstallCommand, 'develop': VicDevelopCommand},
      version=FULLVERSION,
      description='The Variable Infiltration Capacity model',
      author='Joe Hamman',
      author_email='jhamman1@uw.edu',
      tests_require=['pytest'],
      url='https://github.com/UW-Hydro/VIC',
      packages=['vic'],
      py_modules=['vic.vic', 'vic.driver'],
      ext_modules=[ext_module])
# -------------------------------------------------------------------- #
