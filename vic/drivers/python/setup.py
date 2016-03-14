#!/usr/bin/env python
from __future__ import print_function
import os
import re
import sys
import sysconfig
import warnings
import glob
import subprocess
import shutil

from setuptools import setup, find_packages, Command
from setuptools.extension import Extension

# Set the log level
# To turn off warning statements, set LOG_LVL >= 30
# | Level     | Numeric value    |
# |---------  |---------------   |
# | ERROR     | Always Active    |
# | WARNING   | < 30             |
# | INFO      | < 20             |
# | DEBUG     | < 10             |
log_level = 0

MAJOR = 5
MINOR = 0
MICRO = 0
ISRELEASED = False
VERSION = '%d.%d.%d' % (MAJOR, MINOR, MICRO)
QUALIFIER = ''

FULLVERSION = VERSION
write_version = False

vic_root_rel_path = os.path.join(os.pardir, os.pardir, os.pardir)
setup_py_path = os.path.abspath(os.path.dirname(__file__))
vic_root_abs_path = os.path.abspath(
    os.path.join(setup_py_path, os.pardir, os.pardir, os.pardir))
start_dir = os.path.abspath(os.getcwd())
setup_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(setup_dir)


class CleanCommand(Command):
    """Custom clean command to tidy up the project root."""
    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        for path in ['./build', './dist', './__pycache__', './vic.egg-info']:
            try:
                print('removing %s' % path)
                shutil.rmtree(path)
            except:
                pass

        files = ['vic/_vic.py']
        files.extend(glob.glob('vic/*pyc'))
        files.extend(glob.glob('vic_core*'))

        for filename in files:
            try:
                print('removing %s' % filename)
                os.remove(filename)
            except:
                pass


# -------------------------------------------------------------------- #
def write_version_py(filename=None):
    version_text = """\
version = '{0}'
short_version = '{1}'
"""
    if not filename:
        filename = os.path.join(setup_dir, 'vic', 'version.py')

    with open(filename, 'w') as f:
        f.write(version_text.format(FULLVERSION, VERSION))
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
sources.extend(glob.glob(os.path.join(vic_root_abs_path, 'vic', 'vic_run',
                                      'src', '*c')))
sources.extend(glob.glob(os.path.join(vic_root_abs_path, 'vic', 'drivers',
                                      'shared_all', 'src', '*c')))
sources.extend(glob.glob(os.path.join(vic_root_abs_path, 'vic', 'drivers',
                                      'python', 'src', '*c')))

includes = []
includes.append(os.path.join(vic_root_abs_path, 'vic', 'vic_run',
                             'include', ''))
includes.append(os.path.join(vic_root_abs_path, 'vic',
                             'drivers', 'shared_all', 'include', ''))
includes.append(os.path.join(vic_root_abs_path, 'vic',
                             'drivers', 'python', 'include', ''))

ext_name = 'vic_core'
# platform safe path to extension
ext_obj = ext_name + sysconfig.get_config_var('SO')
ext_module = Extension(ext_name,
                       sources=sources,
                       include_dirs=includes,
                       extra_compile_args=['-std=c99',
                                           '-DLOG_LVL={0}'.format(log_level)])

# -------------------------------------------------------------------- #
# Run Setup
setup(name='vic',
      version=FULLVERSION,
      description='The Variable Infiltration Capacity model',
      author='Joe Hamman',
      author_email='jhamman1@uw.edu',
      cmdclass={'clean': CleanCommand},
      setup_requires=["cffi>=1.0.0"],
      install_requires=["cffi>=1.0.0"],
      tests_require=['pytest'],
      url='https://github.com/UW-Hydro/VIC',
      py_modules=["vic"],
      cffi_modules=["vic_build.py:ffi"],
      packages=find_packages(),
      ext_modules=[ext_module])
# -------------------------------------------------------------------- #

os.chdir(start_dir)
