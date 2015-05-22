#!/usr/bin/env python
import os
import re
import sys
import warnings
import glob
try:
    from setuptools import setup
    from setuptools.extension import Extension
except:
    from distutils.core import setup
    from distutils.extension import Extension

MAJOR = 5
MINOR = 0
MICRO = 0
ISRELEASED = False
VERSION = '%d.%d.%d' % (MAJOR, MINOR, MICRO)
QUALIFIER = ''

FULLVERSION = VERSION
write_version = False


vic_root_path = os.path.abspath('../../')


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
# Get version string
if ISRELEASED:
    FULLVERSION += QUALIFIER
else:
    import subprocess
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
sources.extend(glob.glob(os.path.join(vic_root_path, 'vic_run/src/*c')))
sources.extend(glob.glob(os.path.join(vic_root_path, 'drivers/python/src/*c')))

includes = []
includes.append(os.path.join(vic_root_path, 'vic_run/include'))

ext_modules = [Extension('vic_run',
                         sources=sources,
                         include_dirs=includes,
                         extra_compile_args=['-std=c99'])]

# -------------------------------------------------------------------- #
# Run Setup
setup(name='vic',
      version=FULLVERSION,
      description='The Variable Infiltration Capacity model',
      author='Joe Hamman',
      author_email='jhamman1@uw.edu',
      tests_require=['pytest >= 2.5.2'],
      url='https://github.com/UW-Hydro/VIC',
      test_suite='pytest.collector',
      packages=['vic'],
      py_modules=['vic'],
      ext_modules=ext_modules)
# -------------------------------------------------------------------- #
