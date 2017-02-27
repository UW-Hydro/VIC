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
from datetime import datetime

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
MICRO = 1
ISRELEASED = True
VERSION = '%d.%d.%d' % (MAJOR, MINOR, MICRO)
QUALIFIER = ''

FULLVERSION = VERSION
write_version = False
write_headers = True

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
                shutil.rmtree(path)
                print('removed %s' % path)
            except FileNotFoundError:
                pass

        files = ['vic/_vic.py', 'vic_headers.py']
        files.extend(glob.glob('vic/*pyc'))
        files.extend(glob.glob('vic_core*'))

        for filename in files:
            try:
                os.remove(filename)
                print('removed %s' % filename)
            except FileNotFoundError:
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
def maybe_eval_between_brackets(s):
    '''run eval on strings between brackets e.g. [20 + 3] becomes [23]'''
    matches = re.findall("\[(.*?)\]", s)
    for match in matches:
        if match:
            try:
                s = s.replace(match, str(eval(match)))
            except (NameError, SyntaxError):
                pass
    return s
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
def make_cffi_headers():
    '''Process the C headers such that CFFI can interpret them'''

    omissions = ['va_list',
                 'error_print_atmos_energy_bal',
                 'error_print_atmos_moist_bal',
                 'error_print_canopy_energy_bal',
                 'error_print_solve_T_profile',
                 'error_print_surf_energy_bal',
                 'ErrorPrintIcePackEnergyBalance',
                 'ErrorPrintSnowPackEnergyBalance',
                 'func_atmos_energy_bal',
                 'func_atmos_moist_bal',
                 'func_canopy_energy_bal',
                 'func_surf_energy_bal',
                 'IceEnergyBalance',
                 'root_brent',
                 'SnowPackEnergyBalance',
                 'soil_thermal_eqn',
                 'zwtvmoist_zwt',
                 'zwtvmoist_moist']

    args = ['gcc', '-std=c99', '-E',
            '-P', os.path.join(vic_root_abs_path, 'vic', 'drivers',
                               'python', 'src', 'globals.c'),
            '-I%s' % os.path.join(vic_root_abs_path, 'vic', 'drivers',
                                  'python', 'include', ''),
            '-I%s' % os.path.join(vic_root_abs_path, 'vic', 'drivers',
                                  'shared_all', 'include', ''),
            '-I%s' % os.path.join(vic_root_abs_path, 'vic',
                                  'vic_run', 'include', '')]

    proc = subprocess.Popen(' '.join(args),
                            shell=True,
                            stderr=subprocess.PIPE,
                            stdout=subprocess.PIPE)
    retvals = proc.communicate()

    stdout = retvals[0].decode()
    stderr = retvals[1].decode()

    if proc.returncode:
        print(' '.join(args))
        raise RuntimeError('preprocessor experienced an error: %s' % stderr)

    with open(os.path.join(setup_dir, './vic_headers.py'), 'w') as f:
        # write python script header information
        f.write('#!/usr/bin/env python\n')
        f.write("'''\n    Preprocessed Headers for VIC Python Driver\n"
                "    Last updated %s\n'''\n\n\n" % datetime.now())
        # now we write the preprocessed headers, skipping the system headers
        f.write("headers = '''\n")
        skip_headers = True
        for line in stdout.split('\n'):
            # Note: This check for LOG_DEST is here because the subprocess call
            # above includes system headers which we don't want in headers.py.
            # This could be improved by adding a more sophisticated determination
            # of when we've passed the system headers (always at the begining of
            # the standard out stream.
            if 'LOG_DEST' in line:
                skip_headers = False

            # Skip some lines if they contain functions cffi cannot use
            for omit in omissions:
                if omit in line:
                    skip_omit = True
                    break
                else:
                    skip_omit = False
            if not skip_headers and not skip_omit:
                # Evaluate strings that are not completely evaluated by the
                # preprocessor.
                line = maybe_eval_between_brackets(line)
                f.write(line)
                f.write('\n')
        f.write("'''\n")

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

if write_headers:
    make_cffi_headers()
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
