'''Wrapper for vic_def.h

Generated with:
./ctypesgen.py /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_log.h /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_physical_constants.h /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h -o vic_run.py -l ./build/lib.macosx-10.5-x86_64-2.7/vic_run.so

Do not modify this file.
'''

__docformat__ = 'restructuredtext'

# Begin preamble

import ctypes
import os
import sys
from ctypes import *

_int_types = (c_int16, c_int32)
if hasattr(ctypes, 'c_int64'):
    # Some builds of ctypes apparently do not have c_int64
    # defined; it's a pretty good bet that these builds do not
    # have 64-bit pointers.
    _int_types += (c_int64,)
for t in _int_types:
    if sizeof(t) == sizeof(c_size_t):
        c_ptrdiff_t = t
del t
del _int_types


class c_void(Structure):
    # c_void_p is a buggy return type, converting to int, so
    # POINTER(None) == c_void_p is actually written as
    # POINTER(c_void), so it can be treated as a real pointer.
    _fields_ = [('dummy', c_int)]


def POINTER(obj):
    p = ctypes.POINTER(obj)

    # Convert None to a real NULL pointer to work around bugs
    # in how ctypes handles None on 64-bit platforms
    if not isinstance(p.from_param, classmethod):
        def from_param(cls, x):
            if x is None:
                return cls()
            else:
                return x
        p.from_param = classmethod(from_param)

    return p


class UserString:

    def __init__(self, seq):
        if isinstance(seq, basestring):
            self.data = seq
        elif isinstance(seq, UserString):
            self.data = seq.data[:]
        else:
            self.data = str(seq)

    def __str__(self): return str(self.data)

    def __repr__(self): return repr(self.data)

    def __int__(self): return int(self.data)

    def __long__(self): return long(self.data)

    def __float__(self): return float(self.data)

    def __complex__(self): return complex(self.data)

    def __hash__(self): return hash(self.data)

    def __cmp__(self, string):
        if isinstance(string, UserString):
            return cmp(self.data, string.data)
        else:
            return cmp(self.data, string)

    def __contains__(self, char):
        return char in self.data

    def __len__(self): return len(self.data)

    def __getitem__(self, index): return self.__class__(self.data[index])

    def __getslice__(self, start, end):
        start = max(start, 0)
        end = max(end, 0)
        return self.__class__(self.data[start:end])

    def __add__(self, other):
        if isinstance(other, UserString):
            return self.__class__(self.data + other.data)
        elif isinstance(other, basestring):
            return self.__class__(self.data + other)
        else:
            return self.__class__(self.data + str(other))

    def __radd__(self, other):
        if isinstance(other, basestring):
            return self.__class__(other + self.data)
        else:
            return self.__class__(str(other) + self.data)

    def __mul__(self, n):
        return self.__class__(self.data * n)
    __rmul__ = __mul__

    def __mod__(self, args):
        return self.__class__(self.data % args)

    # the following methods are defined in alphabetical order:
    def capitalize(self): return self.__class__(self.data.capitalize())

    def center(self, width, *args):
        return self.__class__(self.data.center(width, *args))

    def count(self, sub, start=0, end=sys.maxsize):
        return self.data.count(sub, start, end)

    def decode(self, encoding=None, errors=None):  # XXX improve this?
        if encoding:
            if errors:
                return self.__class__(self.data.decode(encoding, errors))
            else:
                return self.__class__(self.data.decode(encoding))
        else:
            return self.__class__(self.data.decode())

    def encode(self, encoding=None, errors=None):  # XXX improve this?
        if encoding:
            if errors:
                return self.__class__(self.data.encode(encoding, errors))
            else:
                return self.__class__(self.data.encode(encoding))
        else:
            return self.__class__(self.data.encode())

    def endswith(self, suffix, start=0, end=sys.maxsize):
        return self.data.endswith(suffix, start, end)

    def expandtabs(self, tabsize=8):
        return self.__class__(self.data.expandtabs(tabsize))

    def find(self, sub, start=0, end=sys.maxsize):
        return self.data.find(sub, start, end)

    def index(self, sub, start=0, end=sys.maxsize):
        return self.data.index(sub, start, end)

    def isalpha(self): return self.data.isalpha()

    def isalnum(self): return self.data.isalnum()

    def isdecimal(self): return self.data.isdecimal()

    def isdigit(self): return self.data.isdigit()

    def islower(self): return self.data.islower()

    def isnumeric(self): return self.data.isnumeric()

    def isspace(self): return self.data.isspace()

    def istitle(self): return self.data.istitle()

    def isupper(self): return self.data.isupper()

    def join(self, seq): return self.data.join(seq)

    def ljust(self, width, *args):
        return self.__class__(self.data.ljust(width, *args))

    def lower(self): return self.__class__(self.data.lower())

    def lstrip(self, chars=None): return self.__class__(
        self.data.lstrip(chars))

    def partition(self, sep):
        return self.data.partition(sep)

    def replace(self, old, new, maxsplit=-1):
        return self.__class__(self.data.replace(old, new, maxsplit))

    def rfind(self, sub, start=0, end=sys.maxsize):
        return self.data.rfind(sub, start, end)

    def rindex(self, sub, start=0, end=sys.maxsize):
        return self.data.rindex(sub, start, end)

    def rjust(self, width, *args):
        return self.__class__(self.data.rjust(width, *args))

    def rpartition(self, sep):
        return self.data.rpartition(sep)

    def rstrip(self, chars=None): return self.__class__(
        self.data.rstrip(chars))

    def split(self, sep=None, maxsplit=-1):
        return self.data.split(sep, maxsplit)

    def rsplit(self, sep=None, maxsplit=-1):
        return self.data.rsplit(sep, maxsplit)

    def splitlines(self, keepends=0): return self.data.splitlines(keepends)

    def startswith(self, prefix, start=0, end=sys.maxsize):
        return self.data.startswith(prefix, start, end)

    def strip(self, chars=None): return self.__class__(self.data.strip(chars))

    def swapcase(self): return self.__class__(self.data.swapcase())

    def title(self): return self.__class__(self.data.title())

    def translate(self, *args):
        return self.__class__(self.data.translate(*args))

    def upper(self): return self.__class__(self.data.upper())

    def zfill(self, width): return self.__class__(self.data.zfill(width))


class MutableString(UserString):

    """mutable string objects

    Python strings are immutable objects.  This has the advantage, that
    strings may be used as dictionary keys.  If this property isn't needed
    and you insist on changing string values in place instead, you may cheat
    and use MutableString.

    But the purpose of this class is an educational one: to prevent
    people from inventing their own mutable string class derived
    from UserString and than forget thereby to remove (override) the
    __hash__ method inherited from UserString.  This would lead to
    errors that would be very hard to track down.

    A faster and better solution is to rewrite your program using lists."""

    def __init__(self, string=""):
        self.data = string

    def __hash__(self):
        raise TypeError("unhashable type (it is mutable)")

    def __setitem__(self, index, sub):
        if index < 0:
            index += len(self.data)
        if index < 0 or index >= len(self.data):
            raise IndexError
        self.data = self.data[:index] + sub + self.data[index + 1:]

    def __delitem__(self, index):
        if index < 0:
            index += len(self.data)
        if index < 0 or index >= len(self.data):
            raise IndexError
        self.data = self.data[:index] + self.data[index + 1:]

    def __setslice__(self, start, end, sub):
        start = max(start, 0)
        end = max(end, 0)
        if isinstance(sub, UserString):
            self.data = self.data[:start] + sub.data + self.data[end:]
        elif isinstance(sub, basestring):
            self.data = self.data[:start] + sub + self.data[end:]
        else:
            self.data = self.data[:start] + str(sub) + self.data[end:]

    def __delslice__(self, start, end):
        start = max(start, 0)
        end = max(end, 0)
        self.data = self.data[:start] + self.data[end:]

    def immutable(self):
        return UserString(self.data)

    def __iadd__(self, other):
        if isinstance(other, UserString):
            self.data += other.data
        elif isinstance(other, basestring):
            self.data += other
        else:
            self.data += str(other)
        return self

    def __imul__(self, n):
        self.data *= n
        return self


class String(MutableString, Union):

    _fields_ = [('raw', POINTER(c_char)),
                ('data', c_char_p)]

    def __init__(self, obj=""):
        if isinstance(obj, (str, unicode, UserString)):
            self.data = str(obj)
        else:
            self.raw = obj

    def __len__(self):
        return self.data and len(self.data) or 0

    def from_param(cls, obj):
        # Convert None or 0
        if obj is None or obj == 0:
            return cls(POINTER(c_char)())

        # Convert from String
        elif isinstance(obj, String):
            return obj

        # Convert from str
        elif isinstance(obj, str):
            return cls(obj)

        # Convert from c_char_p
        elif isinstance(obj, c_char_p):
            return obj

        # Convert from POINTER(c_char)
        elif isinstance(obj, POINTER(c_char)):
            return obj

        # Convert from raw pointer
        elif isinstance(obj, int):
            return cls(cast(obj, POINTER(c_char)))

        # Convert from object
        else:
            return String.from_param(obj._as_parameter_)
    from_param = classmethod(from_param)


def ReturnString(obj):
    return String.from_param(obj)

# As of ctypes 1.0, ctypes does not support custom error-checking
# functions on callbacks, nor does it support custom datatypes on
# callbacks, so we must ensure that all callbacks return
# primitive datatypes.
#
# Non-primitive return values wrapped with UNCHECKED won't be
# typechecked, and will be converted to c_void_p.


def UNCHECKED(type):
    if (hasattr(type, "_type_") and isinstance(type._type_, str)
            and type._type_ != "P"):
        return type
    else:
        return c_void_p

# ctypes doesn't have direct support for variadic functions, so we have to write
# our own wrapper class


class _variadic_function(object):

    def __init__(self, func, restype, argtypes):
        self.func = func
        self.func.restype = restype
        self.argtypes = argtypes

    def _as_parameter_(self):
        # So we can pass this variadic function as a function pointer
        return self.func

    def __call__(self, *args):
        fixed_args = []
        i = 0
        for argtype in self.argtypes:
            # Typecheck what we can
            fixed_args.append(argtype.from_param(args[i]))
            i += 1
        return self.func(*fixed_args + list(args[i:]))


# End preamble

_libs = {}
_libdirs = []

# Begin loader

# ----------------------------------------------------------------------------
# Copyright (c) 2008 David James
# Copyright (c) 2006-2008 Alex Holkner
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
#
#  * Redistributions of source code must retain the above copyright
#    notice, this list of conditions and the following disclaimer.
#  * Redistributions in binary form must reproduce the above copyright
#    notice, this list of conditions and the following disclaimer in
#    the documentation and/or other materials provided with the
#    distribution.
#  * Neither the name of pyglet nor the names of its
#    contributors may be used to endorse or promote products
#    derived from this software without specific prior written
#    permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
# FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
# COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
# LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
# ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
# ----------------------------------------------------------------------------

import os.path
import re
import sys
import glob
import ctypes
import ctypes.util


def _environ_path(name):
    if name in os.environ:
        return os.environ[name].split(":")
    else:
        return []


class LibraryLoader(object):

    def __init__(self):
        self.other_dirs = []

    def load_library(self, libname):
        """Given the name of a library, load it."""
        paths = self.getpaths(libname)

        for path in paths:
            if os.path.exists(path):
                return self.load(path)

        raise ImportError("%s not found." % libname)

    def load(self, path):
        """Given a path to a library, load it."""
        try:
            # Darwin requires dlopen to be called with mode RTLD_GLOBAL instead
            # of the default RTLD_LOCAL.  Without this, you end up with
            # libraries not being loadable, resulting in "Symbol not found"
            # errors
            if sys.platform == 'darwin':
                return ctypes.CDLL(path, ctypes.RTLD_GLOBAL)
            else:
                return ctypes.cdll.LoadLibrary(path)
        except OSError as e:
            raise ImportError(e)

    def getpaths(self, libname):
        """Return a list of paths where the library might be found."""
        if os.path.isabs(libname):
            yield libname

        else:
            for path in self.getplatformpaths(libname):
                yield path

            path = ctypes.util.find_library(libname)
            if path:
                yield path

    def getplatformpaths(self, libname):
        return []

# Darwin (Mac OS X)


class DarwinLibraryLoader(LibraryLoader):
    name_formats = ["lib%s.dylib", "lib%s.so", "lib%s.bundle", "%s.dylib",
                    "%s.so", "%s.bundle", "%s"]

    def getplatformpaths(self, libname):
        if os.path.pathsep in libname:
            names = [libname]
        else:
            names = [format % libname for format in self.name_formats]

        for dir in self.getdirs(libname):
            for name in names:
                yield os.path.join(dir, name)

    def getdirs(self, libname):
        '''Implements the dylib search as specified in Apple documentation:

        http://developer.apple.com/documentation/DeveloperTools/Conceptual/
            DynamicLibraries/Articles/DynamicLibraryUsageGuidelines.html

        Before commencing the standard search, the method first checks
        the bundle's ``Frameworks`` directory if the application is running
        within a bundle (OS X .app).
        '''

        dyld_fallback_library_path = _environ_path(
            "DYLD_FALLBACK_LIBRARY_PATH")
        if not dyld_fallback_library_path:
            dyld_fallback_library_path = [os.path.expanduser('~/lib'),
                                          '/usr/local/lib', '/usr/lib']

        dirs = []

        if '/' in libname:
            dirs.extend(_environ_path("DYLD_LIBRARY_PATH"))
        else:
            dirs.extend(_environ_path("LD_LIBRARY_PATH"))
            dirs.extend(_environ_path("DYLD_LIBRARY_PATH"))

        dirs.extend(self.other_dirs)
        dirs.append(".")

        if hasattr(sys, 'frozen') and sys.frozen == 'macosx_app':
            dirs.append(os.path.join(
                os.environ['RESOURCEPATH'],
                '..',
                'Frameworks'))

        dirs.extend(dyld_fallback_library_path)

        return dirs

# Posix


class PosixLibraryLoader(LibraryLoader):
    _ld_so_cache = None

    def _create_ld_so_cache(self):
        # Recreate search path followed by ld.so.  This is going to be
        # slow to build, and incorrect (ld.so uses ld.so.cache, which may
        # not be up-to-date).  Used only as fallback for distros without
        # /sbin/ldconfig.
        #
        # We assume the DT_RPATH and DT_RUNPATH binary sections are omitted.

        directories = []
        for name in ("LD_LIBRARY_PATH",
                     "SHLIB_PATH",  # HPUX
                     "LIBPATH",  # OS/2, AIX
                     "LIBRARY_PATH",  # BE/OS
                     ):
            if name in os.environ:
                directories.extend(os.environ[name].split(os.pathsep))
        directories.extend(self.other_dirs)
        directories.append(".")

        try:
            directories.extend([dir.strip()
                                for dir in open('/etc/ld.so.conf')])
        except IOError:
            pass

        directories.extend(['/lib', '/usr/lib', '/lib64', '/usr/lib64'])

        cache = {}
        lib_re = re.compile(r'lib(.*)\.s[ol]')
        ext_re = re.compile(r'\.s[ol]$')
        for dir in directories:
            try:
                for path in glob.glob("%s/*.s[ol]*" % dir):
                    file = os.path.basename(path)

                    # Index by filename
                    if file not in cache:
                        cache[file] = path

                    # Index by library name
                    match = lib_re.match(file)
                    if match:
                        library = match.group(1)
                        if library not in cache:
                            cache[library] = path
            except OSError:
                pass

        self._ld_so_cache = cache

    def getplatformpaths(self, libname):
        if self._ld_so_cache is None:
            self._create_ld_so_cache()

        result = self._ld_so_cache.get(libname)
        if result:
            yield result

        path = ctypes.util.find_library(libname)
        if path:
            yield os.path.join("/lib", path)

# Windows


class _WindowsLibrary(object):

    def __init__(self, path):
        self.cdll = ctypes.cdll.LoadLibrary(path)
        self.windll = ctypes.windll.LoadLibrary(path)

    def __getattr__(self, name):
        try:
            return getattr(self.cdll, name)
        except AttributeError:
            try:
                return getattr(self.windll, name)
            except AttributeError:
                raise


class WindowsLibraryLoader(LibraryLoader):
    name_formats = ["%s.dll", "lib%s.dll"]

    def load(self, path):
        return _WindowsLibrary(path)

    def getplatformpaths(self, libname):
        if os.path.sep not in libname:
            for name in self.name_formats:
                path = ctypes.util.find_library(name % libname)
                if path:
                    yield path

# Platform switching

# If your value of sys.platform does not appear in this dict, please contact
# the Ctypesgen maintainers.

loaderclass = {
    "darwin": DarwinLibraryLoader,
    "cygwin": WindowsLibraryLoader,
    "win32": WindowsLibraryLoader
}

loader = loaderclass.get(sys.platform, PosixLibraryLoader)()


def add_library_search_dirs(other_dirs):
    loader.other_dirs = other_dirs

load_library = loader.load_library

del loaderclass

# End loader

add_library_search_dirs([])

# Begin libraries

_libs["./build/lib.macosx-10.5-x86_64-2.7/vic_run.so"] = load_library(
    "./build/lib.macosx-10.5-x86_64-2.7/vic_run.so")

# 1 libraries
# End libraries

# No modules

__int64_t = c_longlong  # /usr/include/i386/_types.h: 46

__darwin_off_t = __int64_t  # /usr/include/sys/_types.h: 71

fpos_t = __darwin_off_t  # /usr/include/stdio.h: 77

# /usr/include/stdio.h: 88


class struct___sbuf(Structure):
    pass

struct___sbuf.__slots__ = [
    '_base',
    '_size',
]
struct___sbuf._fields_ = [
    ('_base', POINTER(c_ubyte)),
    ('_size', c_int),
]

# /usr/include/stdio.h: 94


class struct___sFILEX(Structure):
    pass

# /usr/include/stdio.h: 153


class struct___sFILE(Structure):
    pass

struct___sFILE.__slots__ = [
    '_p',
    '_r',
    '_w',
    '_flags',
    '_file',
    '_bf',
    '_lbfsize',
    '_cookie',
    '_close',
    '_read',
    '_seek',
    '_write',
    '_ub',
    '_extra',
    '_ur',
    '_ubuf',
    '_nbuf',
    '_lb',
    '_blksize',
    '_offset',
]
struct___sFILE._fields_ = [
    ('_p', POINTER(c_ubyte)),
    ('_r', c_int),
    ('_w', c_int),
    ('_flags', c_short),
    ('_file', c_short),
    ('_bf', struct___sbuf),
    ('_lbfsize', c_int),
    ('_cookie', POINTER(None)),
    ('_close', CFUNCTYPE(UNCHECKED(c_int), POINTER(None))),
    ('_read', CFUNCTYPE(UNCHECKED(c_int), POINTER(None), String, c_int)),
    ('_seek', CFUNCTYPE(UNCHECKED(fpos_t), POINTER(None), fpos_t, c_int)),
    ('_write', CFUNCTYPE(UNCHECKED(c_int), POINTER(None), String, c_int)),
    ('_ub', struct___sbuf),
    ('_extra', POINTER(struct___sFILEX)),
    ('_ur', c_int),
    ('_ubuf', c_ubyte * 3),
    ('_nbuf', c_ubyte * 1),
    ('_lb', struct___sbuf),
    ('_blksize', c_int),
    ('_offset', fpos_t),
]

FILE = struct___sFILE  # /usr/include/stdio.h: 153

# /usr/include/string.h: 81
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'strerror'):
    strerror = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].strerror
    strerror.restype = ReturnString
    strerror.argtypes = [c_int]

# /usr/include/sys/errno.h: 80
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        '__error'):
    __error = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].__error
    __error.restype = POINTER(c_int)
    __error.argtypes = []

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_log.h: 48
try:
    LOG_DEST = (
        POINTER(FILE)).in_dll(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'LOG_DEST')
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_log.h: 57
for _lib in _libs.values():
    if hasattr(_lib, 'finalize_logging'):
        finalize_logging = _lib.finalize_logging
        finalize_logging.restype = None
        finalize_logging.argtypes = []
        break

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_log.h: 58
for _lib in _libs.values():
    if hasattr(_lib, 'get_current_datetime'):
        get_current_datetime = _lib.get_current_datetime
        get_current_datetime.restype = None
        get_current_datetime.argtypes = [String]
        break

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_log.h: 59
for _lib in _libs.values():
    if hasattr(_lib, 'get_logname'):
        get_logname = _lib.get_logname
        get_logname.restype = None
        get_logname.argtypes = [String, c_int, String]
        break

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_log.h: 60
for _lib in _libs.values():
    if hasattr(_lib, 'initialize_log'):
        initialize_log = _lib.initialize_log
        initialize_log.restype = None
        initialize_log.argtypes = []
        break

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_log.h: 61
for _lib in _libs.values():
    if hasattr(_lib, 'setup_logging'):
        setup_logging = _lib.setup_logging
        setup_logging.restype = None
        setup_logging.argtypes = [c_int]
        break

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 77
try:
    ref_veg_over = (
        POINTER(c_uint8)).in_dll(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'ref_veg_over')
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 78
try:
    ref_veg_rarc = (
        POINTER(c_double)).in_dll(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'ref_veg_rarc')
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 79
try:
    ref_veg_rmin = (
        POINTER(c_double)).in_dll(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'ref_veg_rmin')
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 80
try:
    ref_veg_lai = (
        POINTER(c_double)).in_dll(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'ref_veg_lai')
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 81
try:
    ref_veg_albedo = (
        POINTER(c_double)).in_dll(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'ref_veg_albedo')
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 82
try:
    ref_veg_vegcover = (
        POINTER(c_double)).in_dll(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'ref_veg_vegcover')
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 83
try:
    ref_veg_rough = (
        POINTER(c_double)).in_dll(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'ref_veg_rough')
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 84
try:
    ref_veg_displ = (
        POINTER(c_double)).in_dll(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'ref_veg_displ')
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 85
try:
    ref_veg_wind_h = (
        POINTER(c_double)).in_dll(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'ref_veg_wind_h')
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 86
try:
    ref_veg_RGL = (
        POINTER(c_double)).in_dll(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'ref_veg_RGL')
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 87
try:
    ref_veg_rad_atten = (
        POINTER(c_double)).in_dll(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'ref_veg_rad_atten')
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 88
try:
    ref_veg_wind_atten = (
        POINTER(c_double)).in_dll(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'ref_veg_wind_atten')
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 89
try:
    ref_veg_trunk_ratio = (
        POINTER(c_double)).in_dll(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'ref_veg_trunk_ratio')
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 90
try:
    ref_veg_ref_crop = (
        POINTER(c_uint8)).in_dll(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'ref_veg_ref_crop')
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 105
try:
    NR = (c_size_t).in_dll(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'NR')
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 107
try:
    NF = (c_size_t).in_dll(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'NF')
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 113
enum_anon_8 = c_int

ASCII = 0  # /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 113

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 113
BINARY = (ASCII + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 122
enum_anon_9 = c_int

LITTLE = 0  # /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 122

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 122
BIG = (LITTLE + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 131
enum_anon_10 = c_int

DENS_BRAS = 0  # /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 131

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 131
DENS_SNTHRM = (DENS_BRAS + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 140
enum_anon_11 = c_int

ARNO = 0  # /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 140

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 140
NIJSSEN2001 = (ARNO + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 149
enum_anon_12 = c_int

AR_406 = 0  # /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 149

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 149
AR_406_LS = (AR_406 + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 149
AR_406_FULL = (AR_406_LS + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 149
AR_410 = (AR_406_FULL + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 160
enum_anon_13 = c_int

GF_406 = 0  # /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 160

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 160
GF_410 = (GF_406 + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 169
enum_anon_14 = c_int

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 169
VP_ITER_NONE = 0

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 169
VP_ITER_ALWAYS = (VP_ITER_NONE + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 169
VP_ITER_ANNUAL = (VP_ITER_ALWAYS + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 169
VP_ITER_CONVERGE = (VP_ITER_ANNUAL + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 180
enum_anon_15 = c_int

LW_TVA = 0  # /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 180

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 180
LW_ANDERSON = (LW_TVA + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 180
LW_BRUTSAERT = (LW_ANDERSON + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 180
LW_SATTERLUND = (LW_BRUTSAERT + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 180
LW_IDSO = (LW_SATTERLUND + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 180
LW_PRATA = (LW_IDSO + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 193
enum_anon_16 = c_int

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 193
LW_CLOUD_BRAS = 0

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 193
LW_CLOUD_DEARDORFF = (LW_CLOUD_BRAS + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 202
enum_anon_17 = c_int

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 202
FROM_VEGLIB = 0

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 202
FROM_VEGPARAM = (FROM_VEGLIB + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 211
enum_anon_18 = c_int

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 211
LAI_FROM_VEGLIB = 0

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 211
LAI_FROM_VEGPARAM = (LAI_FROM_VEGLIB + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 220
enum_anon_19 = c_int

RC_JARVIS = 0  # /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 220

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 220
RC_PHOTO = (RC_JARVIS + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 229
enum_anon_20 = c_int

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 229
PS_FARQUHAR = 0

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 229
PS_MONTEITH = (PS_FARQUHAR + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 238
enum_anon_21 = c_int

PHOTO_C3 = 0  # /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 238

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 238
PHOTO_C4 = (PHOTO_C3 + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 247
enum_anon_22 = c_int

AIR_TEMP = 0  # /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 247

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 247
ALBEDO = (AIR_TEMP + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 247
CATM = (ALBEDO + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 247
CHANNEL_IN = (CATM + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 247
CRAINF = (CHANNEL_IN + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 247
CSNOWF = (CRAINF + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 247
DENSITY = (CSNOWF + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 247
FDIR = (DENSITY + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 247
LAI_IN = (FDIR + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 247
LONGWAVE = (LAI_IN + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 247
LSRAINF = (LONGWAVE + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 247
LSSNOWF = (LSRAINF + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 247
PAR = (LSSNOWF + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 247
PREC = (PAR + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 247
PRESSURE = (PREC + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 247
QAIR = (PRESSURE + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 247
RAINF = (QAIR + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 247
REL_HUMID = (RAINF + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 247
SHORTWAVE = (REL_HUMID + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 247
SNOWF = (SHORTWAVE + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 247
TMAX = (SNOWF + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 247
TMIN = (TMAX + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 247
TSKC = (TMIN + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 247
VEGCOVER = (TSKC + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 247
VP = (VEGCOVER + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 247
WIND = (VP + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 247
WIND_E = (WIND + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 247
WIND_N = (WIND_E + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 247
SKIP = (WIND_N + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 247
N_FORCING_TYPES = (SKIP + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
enum_anon_23 = c_int

OUT_ASAT = 0  # /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_LAKE_AREA_FRAC = (OUT_ASAT + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_LAKE_DEPTH = (OUT_LAKE_AREA_FRAC + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_LAKE_ICE = (OUT_LAKE_DEPTH + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_LAKE_ICE_FRACT = (OUT_LAKE_ICE + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_LAKE_ICE_HEIGHT = (OUT_LAKE_ICE_FRACT + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_LAKE_MOIST = (OUT_LAKE_ICE_HEIGHT + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_LAKE_SURF_AREA = (OUT_LAKE_MOIST + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_LAKE_SWE = (OUT_LAKE_SURF_AREA + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_LAKE_SWE_V = (OUT_LAKE_SWE + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_LAKE_VOLUME = (OUT_LAKE_SWE_V + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_ROOTMOIST = (OUT_LAKE_VOLUME + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_SMFROZFRAC = (OUT_ROOTMOIST + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_SMLIQFRAC = (OUT_SMFROZFRAC + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_SNOW_CANOPY = (OUT_SMLIQFRAC + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_SNOW_COVER = (OUT_SNOW_CANOPY + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_SNOW_DEPTH = (OUT_SNOW_COVER + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_SOIL_ICE = (OUT_SNOW_DEPTH + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_SOIL_LIQ = (OUT_SOIL_ICE + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_SOIL_MOIST = (OUT_SOIL_LIQ + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_SOIL_WET = (OUT_SOIL_MOIST + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_SURFSTOR = (OUT_SOIL_WET + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_SURF_FROST_FRAC = (OUT_SURFSTOR + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_SWE = (OUT_SURF_FROST_FRAC + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_WDEW = (OUT_SWE + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_ZWT = (OUT_WDEW + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_ZWT_LUMPED = (OUT_ZWT + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_BASEFLOW = (OUT_ZWT_LUMPED + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_DELINTERCEPT = (OUT_BASEFLOW + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_DELSOILMOIST = (OUT_DELINTERCEPT + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_DELSURFSTOR = (OUT_DELSOILMOIST + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_DELSWE = (OUT_DELSURFSTOR + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_EVAP = (OUT_DELSWE + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_EVAP_BARE = (OUT_EVAP + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_EVAP_CANOP = (OUT_EVAP_BARE + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_INFLOW = (OUT_EVAP_CANOP + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_LAKE_BF_IN = (OUT_INFLOW + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_LAKE_BF_IN_V = (OUT_LAKE_BF_IN + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_LAKE_BF_OUT = (OUT_LAKE_BF_IN_V + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_LAKE_BF_OUT_V = (OUT_LAKE_BF_OUT + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_LAKE_CHAN_IN = (OUT_LAKE_BF_OUT_V + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_LAKE_CHAN_IN_V = (OUT_LAKE_CHAN_IN + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_LAKE_CHAN_OUT = (OUT_LAKE_CHAN_IN_V + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_LAKE_CHAN_OUT_V = (OUT_LAKE_CHAN_OUT + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_LAKE_DSTOR = (OUT_LAKE_CHAN_OUT_V + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_LAKE_DSTOR_V = (OUT_LAKE_DSTOR + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_LAKE_DSWE = (OUT_LAKE_DSTOR_V + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_LAKE_DSWE_V = (OUT_LAKE_DSWE + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_LAKE_EVAP = (OUT_LAKE_DSWE_V + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_LAKE_EVAP_V = (OUT_LAKE_EVAP + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_LAKE_PREC_V = (OUT_LAKE_EVAP_V + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_LAKE_RCHRG = (OUT_LAKE_PREC_V + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_LAKE_RCHRG_V = (OUT_LAKE_RCHRG + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_LAKE_RO_IN = (OUT_LAKE_RCHRG_V + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_LAKE_RO_IN_V = (OUT_LAKE_RO_IN + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_LAKE_VAPFLX = (OUT_LAKE_RO_IN_V + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_LAKE_VAPFLX_V = (OUT_LAKE_VAPFLX + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_PET_SATSOIL = (OUT_LAKE_VAPFLX_V + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_PET_H2OSURF = (OUT_PET_SATSOIL + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_PET_SHORT = (OUT_PET_H2OSURF + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_PET_TALL = (OUT_PET_SHORT + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_PET_NATVEG = (OUT_PET_TALL + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_PET_VEGNOCR = (OUT_PET_NATVEG + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_PREC = (OUT_PET_VEGNOCR + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_RAINF = (OUT_PREC + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_REFREEZE = (OUT_RAINF + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_RUNOFF = (OUT_REFREEZE + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_SNOW_MELT = (OUT_RUNOFF + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_SNOWF = (OUT_SNOW_MELT + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_SUB_BLOWING = (OUT_SNOWF + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_SUB_CANOP = (OUT_SUB_BLOWING + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_SUB_SNOW = (OUT_SUB_CANOP + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_SUB_SURFACE = (OUT_SUB_SNOW + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_TRANSP_VEG = (OUT_SUB_SURFACE + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_WATER_ERROR = (OUT_TRANSP_VEG + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_ALBEDO = (OUT_WATER_ERROR + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_BARESOILT = (OUT_ALBEDO + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_FDEPTH = (OUT_BARESOILT + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_LAKE_ICE_TEMP = (OUT_FDEPTH + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_LAKE_SURF_TEMP = (OUT_LAKE_ICE_TEMP + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_RAD_TEMP = (OUT_LAKE_SURF_TEMP + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_SALBEDO = (OUT_RAD_TEMP + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_SNOW_PACK_TEMP = (OUT_SALBEDO + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_SNOW_SURF_TEMP = (OUT_SNOW_PACK_TEMP + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_SNOWT_FBFLAG = (OUT_SNOW_SURF_TEMP + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_SOIL_TEMP = (OUT_SNOWT_FBFLAG + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_SOIL_TNODE = (OUT_SOIL_TEMP + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_SOIL_TNODE_WL = (OUT_SOIL_TNODE + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_SOILT_FBFLAG = (OUT_SOIL_TNODE_WL + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_SURF_TEMP = (OUT_SOILT_FBFLAG + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_SURFT_FBFLAG = (OUT_SURF_TEMP + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_TCAN_FBFLAG = (OUT_SURFT_FBFLAG + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_TDEPTH = (OUT_TCAN_FBFLAG + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_TFOL_FBFLAG = (OUT_TDEPTH + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_VEGT = (OUT_TFOL_FBFLAG + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_ADV_SENS = (OUT_VEGT + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_ADVECTION = (OUT_ADV_SENS + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_DELTACC = (OUT_ADVECTION + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_DELTAH = (OUT_DELTACC + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_ENERGY_ERROR = (OUT_DELTAH + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_FUSION = (OUT_ENERGY_ERROR + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_GRND_FLUX = (OUT_FUSION + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_IN_LONG = (OUT_GRND_FLUX + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_LATENT = (OUT_IN_LONG + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_LATENT_SUB = (OUT_LATENT + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_MELT_ENERGY = (OUT_LATENT_SUB + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_NET_LONG = (OUT_MELT_ENERGY + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_NET_SHORT = (OUT_NET_LONG + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_R_NET = (OUT_NET_SHORT + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_RFRZ_ENERGY = (OUT_R_NET + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_SENSIBLE = (OUT_RFRZ_ENERGY + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_SNOW_FLUX = (OUT_SENSIBLE + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_AERO_COND = (OUT_SNOW_FLUX + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_AERO_COND1 = (OUT_AERO_COND + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_AERO_COND2 = (OUT_AERO_COND1 + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_AERO_RESIST = (OUT_AERO_COND2 + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_AERO_RESIST1 = (OUT_AERO_RESIST + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_AERO_RESIST2 = (OUT_AERO_RESIST1 + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_AIR_TEMP = (OUT_AERO_RESIST2 + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_CATM = (OUT_AIR_TEMP + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_COSZEN = (OUT_CATM + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_DENSITY = (OUT_COSZEN + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_FDIR = (OUT_DENSITY + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_LAI = (OUT_FDIR + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_LONGWAVE = (OUT_LAI + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_PAR = (OUT_LONGWAVE + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_PRESSURE = (OUT_PAR + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_QAIR = (OUT_PRESSURE + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_REL_HUMID = (OUT_QAIR + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_SHORTWAVE = (OUT_REL_HUMID + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_SURF_COND = (OUT_SHORTWAVE + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_TSKC = (OUT_SURF_COND + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_VEGCOVER = (OUT_TSKC + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_VP = (OUT_VEGCOVER + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_VPD = (OUT_VP + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_WIND = (OUT_VPD + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_ADV_SENS_BAND = (OUT_WIND + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_ADVECTION_BAND = (OUT_ADV_SENS_BAND + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_ALBEDO_BAND = (OUT_ADVECTION_BAND + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_DELTACC_BAND = (OUT_ALBEDO_BAND + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_GRND_FLUX_BAND = (OUT_DELTACC_BAND + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_IN_LONG_BAND = (OUT_GRND_FLUX_BAND + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_LATENT_BAND = (OUT_IN_LONG_BAND + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_LATENT_SUB_BAND = (OUT_LATENT_BAND + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_MELT_ENERGY_BAND = (OUT_LATENT_SUB_BAND + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_NET_LONG_BAND = (OUT_MELT_ENERGY_BAND + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_NET_SHORT_BAND = (OUT_NET_LONG_BAND + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_RFRZ_ENERGY_BAND = (OUT_NET_SHORT_BAND + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_SENSIBLE_BAND = (OUT_RFRZ_ENERGY_BAND + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_SNOW_CANOPY_BAND = (OUT_SENSIBLE_BAND + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_SNOW_COVER_BAND = (OUT_SNOW_CANOPY_BAND + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_SNOW_DEPTH_BAND = (OUT_SNOW_COVER_BAND + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_SNOW_FLUX_BAND = (OUT_SNOW_DEPTH_BAND + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_SNOW_MELT_BAND = (OUT_SNOW_FLUX_BAND + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_SNOW_PACKT_BAND = (OUT_SNOW_MELT_BAND + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_SNOW_SURFT_BAND = (OUT_SNOW_PACKT_BAND + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_SWE_BAND = (OUT_SNOW_SURFT_BAND + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_APAR = (OUT_SWE_BAND + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_GPP = (OUT_APAR + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_RAUT = (OUT_GPP + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_NPP = (OUT_RAUT + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_LITTERFALL = (OUT_NPP + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_RHET = (OUT_LITTERFALL + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_NEE = (OUT_RHET + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_CLITTER = (OUT_NEE + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_CINTER = (OUT_CLITTER + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
OUT_CSLOW = (OUT_CINTER + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 286
N_OUTVAR_TYPES = (OUT_CSLOW + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 470
enum_anon_24 = c_int

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 470
OUT_TYPE_DEFAULT = 0

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 470
OUT_TYPE_CHAR = (OUT_TYPE_DEFAULT + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 470
OUT_TYPE_SINT = (OUT_TYPE_CHAR + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 470
OUT_TYPE_USINT = (OUT_TYPE_SINT + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 470
OUT_TYPE_INT = (OUT_TYPE_USINT + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 470
OUT_TYPE_FLOAT = (OUT_TYPE_INT + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 470
OUT_TYPE_DOUBLE = (OUT_TYPE_FLOAT + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 484
enum_anon_25 = c_int

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 484
AGG_TYPE_AVG = 0

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 484
AGG_TYPE_BEG = (AGG_TYPE_AVG + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 484
AGG_TYPE_END = (AGG_TYPE_BEG + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 484
AGG_TYPE_MAX = (AGG_TYPE_END + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 484
AGG_TYPE_MIN = (AGG_TYPE_MAX + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 484
AGG_TYPE_SUM = (AGG_TYPE_MIN + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 498
enum_anon_26 = c_int

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 498
DISP_VERSION = 0

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 498
DISP_COMPILE_TIME = (DISP_VERSION + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 498
DISP_ALL = (DISP_COMPILE_TIME + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 508
enum_calendars = c_int

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 508
CALENDAR_STANDARD = 0

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 508
CALENDAR_GREGORIAN = (CALENDAR_STANDARD + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 508
CALENDAR_PROLEPTIC_GREGORIAN = (CALENDAR_GREGORIAN + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 508
CALENDAR_NOLEAP = (CALENDAR_PROLEPTIC_GREGORIAN + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 508
CALENDAR_365_DAY = (CALENDAR_NOLEAP + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 508
CALENDAR_360_DAY = (CALENDAR_365_DAY + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 508
CALENDAR_JULIAN = (CALENDAR_360_DAY + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 508
CALENDAR_ALL_LEAP = (CALENDAR_JULIAN + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 508
CALENDAR_366_DAY = (CALENDAR_ALL_LEAP + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 524
enum_time_units = c_int

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 524
TIME_UNITS_SECONDS = 0

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 524
TIME_UNITS_MINUTES = (TIME_UNITS_SECONDS + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 524
TIME_UNITS_HOURS = (TIME_UNITS_MINUTES + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 524
TIME_UNITS_DAYS = (TIME_UNITS_HOURS + 1)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 550


class struct_anon_27(Structure):
    pass

struct_anon_27.__slots__ = [
    'forcing',
    'globalparam',
    'constants',
    'domain',
    'init_state',
    'lakeparam',
    'snowband',
    'soilparam',
    'statefile',
    'veglib',
    'vegparam',
    'logfile',
]
struct_anon_27._fields_ = [
    ('forcing', POINTER(FILE) * 2),
    ('globalparam', POINTER(FILE)),
    ('constants', POINTER(FILE)),
    ('domain', POINTER(FILE)),
    ('init_state', POINTER(FILE)),
    ('lakeparam', POINTER(FILE)),
    ('snowband', POINTER(FILE)),
    ('soilparam', POINTER(FILE)),
    ('statefile', POINTER(FILE)),
    ('veglib', POINTER(FILE)),
    ('vegparam', POINTER(FILE)),
    ('logfile', POINTER(FILE)),
]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 550
filep_struct = struct_anon_27

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 570


class struct_anon_28(Structure):
    pass

struct_anon_28.__slots__ = [
    'forcing',
    'f_path_pfx',
    '_global',
    'domain',
    'constants',
    'init_state',
    'lakeparam',
    'result_dir',
    'snowband',
    'soil',
    'statefile',
    'veg',
    'veglib',
    'log_path',
]
struct_anon_28._fields_ = [
    ('forcing', (c_char * 2048) * 2),
    ('f_path_pfx', (c_char * 2048) * 2),
    ('_global', c_char * 2048),
    ('domain', c_char * 2048),
    ('constants', c_char * 2048),
    ('init_state', c_char * 2048),
    ('lakeparam', c_char * 2048),
    ('result_dir', c_char * 2048),
    ('snowband', c_char * 2048),
    ('soil', c_char * 2048),
    ('statefile', c_char * 2048),
    ('veg', c_char * 2048),
    ('veglib', c_char * 2048),
    ('log_path', c_char * 2048),
]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 570
filenames_struct = struct_anon_28

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 724


class struct_anon_29(Structure):
    pass

struct_anon_29.__slots__ = [
    'AboveTreelineVeg',
    'AERO_RESIST_CANSNOW',
    'BLOWING',
    'BLOWING_VAR_THRESHOLD',
    'BLOWING_CALC_PROB',
    'BLOWING_SIMPLE',
    'BLOWING_FETCH',
    'BLOWING_SPATIAL_WIND',
    'CARBON',
    'CLOSE_ENERGY',
    'COMPUTE_TREELINE',
    'CONTINUEONERROR',
    'CORRPREC',
    'EQUAL_AREA',
    'EXP_TRANS',
    'FROZEN_SOIL',
    'FULL_ENERGY',
    'GRND_FLUX_TYPE',
    'IMPLICIT',
    'JULY_TAVG_SUPPLIED',
    'LAKES',
    'LW_CLOUD',
    'LW_TYPE',
    'MTCLIM_SWE_CORR',
    'Ncanopy',
    'Nfrost',
    'Nlakenode',
    'Nlayer',
    'Nnode',
    'NOFLUX',
    'NVEGTYPES',
    'PLAPSE',
    'RC_MODE',
    'ROOT_ZONES',
    'QUICK_FLUX',
    'QUICK_SOLVE',
    'SHARE_LAYER_MOIST',
    'SNOW_DENSITY',
    'SNOW_BAND',
    'SPATIAL_FROST',
    'SPATIAL_SNOW',
    'TFALLBACK',
    'VP_INTERP',
    'VP_ITER',
    'ALMA_INPUT',
    'BASEFLOW',
    'GRID_DECIMAL',
    'VEGLIB_PHOTO',
    'VEGLIB_VEGCOVER',
    'VEGPARAM_ALB',
    'VEGPARAM_LAI',
    'VEGPARAM_VEGCOVER',
    'ALB_SRC',
    'LAI_SRC',
    'VEGCOVER_SRC',
    'LAKE_PROFILE',
    'ORGANIC_FRACT',
    'BINARY_STATE_FILE',
    'INIT_STATE',
    'SAVE_STATE',
    'ALMA_OUTPUT',
    'BINARY_OUTPUT',
    'COMPRESS',
    'MOISTFRACT',
    'Noutfiles',
    'OUTPUT_FORCE',
    'PRT_HEADER',
    'PRT_SNOW_BAND',
]
struct_anon_29._fields_ = [
    ('AboveTreelineVeg', c_short),
    ('AERO_RESIST_CANSNOW', c_uint),
    ('BLOWING', c_uint8),
    ('BLOWING_VAR_THRESHOLD', c_uint8),
    ('BLOWING_CALC_PROB', c_uint8),
    ('BLOWING_SIMPLE', c_uint8),
    ('BLOWING_FETCH', c_uint8),
    ('BLOWING_SPATIAL_WIND', c_uint8),
    ('CARBON', c_uint8),
    ('CLOSE_ENERGY', c_uint8),
    ('COMPUTE_TREELINE', c_uint8),
    ('CONTINUEONERROR', c_uint8),
    ('CORRPREC', c_uint8),
    ('EQUAL_AREA', c_uint8),
    ('EXP_TRANS', c_uint8),
    ('FROZEN_SOIL', c_uint8),
    ('FULL_ENERGY', c_uint8),
    ('GRND_FLUX_TYPE', c_uint),
    ('IMPLICIT', c_uint8),
    ('JULY_TAVG_SUPPLIED', c_uint8),
    ('LAKES', c_uint8),
    ('LW_CLOUD', c_uint),
    ('LW_TYPE', c_uint),
    ('MTCLIM_SWE_CORR', c_uint8),
    ('Ncanopy', c_size_t),
    ('Nfrost', c_size_t),
    ('Nlakenode', c_size_t),
    ('Nlayer', c_size_t),
    ('Nnode', c_size_t),
    ('NOFLUX', c_uint8),
    ('NVEGTYPES', c_size_t),
    ('PLAPSE', c_uint8),
    ('RC_MODE', c_uint),
    ('ROOT_ZONES', c_size_t),
    ('QUICK_FLUX', c_uint8),
    ('QUICK_SOLVE', c_uint8),
    ('SHARE_LAYER_MOIST', c_uint8),
    ('SNOW_DENSITY', c_uint),
    ('SNOW_BAND', c_size_t),
    ('SPATIAL_FROST', c_uint8),
    ('SPATIAL_SNOW', c_uint8),
    ('TFALLBACK', c_uint8),
    ('VP_INTERP', c_uint8),
    ('VP_ITER', c_ushort),
    ('ALMA_INPUT', c_uint8),
    ('BASEFLOW', c_uint8),
    ('GRID_DECIMAL', c_uint),
    ('VEGLIB_PHOTO', c_uint8),
    ('VEGLIB_VEGCOVER', c_uint8),
    ('VEGPARAM_ALB', c_uint8),
    ('VEGPARAM_LAI', c_uint8),
    ('VEGPARAM_VEGCOVER', c_uint8),
    ('ALB_SRC', c_uint),
    ('LAI_SRC', c_uint),
    ('VEGCOVER_SRC', c_uint),
    ('LAKE_PROFILE', c_uint8),
    ('ORGANIC_FRACT', c_uint8),
    ('BINARY_STATE_FILE', c_uint8),
    ('INIT_STATE', c_uint8),
    ('SAVE_STATE', c_uint8),
    ('ALMA_OUTPUT', c_uint8),
    ('BINARY_OUTPUT', c_uint8),
    ('COMPRESS', c_uint8),
    ('MOISTFRACT', c_uint8),
    ('Noutfiles', c_size_t),
    ('OUTPUT_FORCE', c_uint8),
    ('PRT_HEADER', c_uint8),
    ('PRT_SNOW_BAND', c_uint8),
]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 724
option_struct = struct_anon_29

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 772


class struct_anon_30(Structure):
    pass

struct_anon_30.__slots__ = [
    'measure_h',
    'wind_h',
    'resolution',
    'dt',
    'snow_dt',
    'runoff_dt',
    'atmos_dt',
    'out_dt',
    'model_steps_per_day',
    'snow_steps_per_day',
    'runoff_steps_per_day',
    'atmos_steps_per_day',
    'output_steps_per_day',
    'endday',
    'endmonth',
    'endyear',
    'forceday',
    'forcesec',
    'forcemonth',
    'forceoffset',
    'forceskip',
    'forceyear',
    'nrecs',
    'skipyear',
    'startday',
    'startsec',
    'startmonth',
    'startyear',
    'stateday',
    'statemonth',
    'stateyear',
    'calendar',
    'time_units',
    'time_origin_num',
]
struct_anon_30._fields_ = [
    ('measure_h', c_double),
    ('wind_h', c_double),
    ('resolution', c_double),
    ('dt', c_double),
    ('snow_dt', c_double),
    ('runoff_dt', c_double),
    ('atmos_dt', c_double),
    ('out_dt', c_double),
    ('model_steps_per_day', c_size_t),
    ('snow_steps_per_day', c_size_t),
    ('runoff_steps_per_day', c_size_t),
    ('atmos_steps_per_day', c_size_t),
    ('output_steps_per_day', c_size_t),
    ('endday', c_uint),
    ('endmonth', c_uint),
    ('endyear', c_uint),
    ('forceday', c_uint * 2),
    ('forcesec', c_uint * 2),
    ('forcemonth', c_uint * 2),
    ('forceoffset', c_uint * 2),
    ('forceskip', c_uint * 2),
    ('forceyear', c_uint * 2),
    ('nrecs', c_size_t),
    ('skipyear', c_uint),
    ('startday', c_uint),
    ('startsec', c_uint),
    ('startmonth', c_uint),
    ('startyear', c_uint),
    ('stateday', c_uint),
    ('statemonth', c_uint),
    ('stateyear', c_uint),
    ('calendar', c_uint),
    ('time_units', c_uint),
    ('time_origin_num', c_double),
]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 772
global_param_struct = struct_anon_30

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 982


class struct_anon_31(Structure):
    pass

struct_anon_31.__slots__ = [
    'LAPSE_RATE',
    'GAUGE_HEIGHT',
    'WIND_SPEED_DEFAULT',
    'WIND_SPEED_MIN',
    'HUGE_RESIST',
    'ALBEDO_BARE_SOIL',
    'ALBEDO_H20_SURF',
    'EMISS_GRND',
    'EMISS_VEG',
    'EMISS_ICE',
    'EMISS_SNOW',
    'EMISS_H2O',
    'SOIL_RESID_MOIST',
    'SOIL_SLAB_MOIST_FRACT',
    'VEG_LAI_SNOW_MULTIPLIER',
    'VEG_MIN_INTERCEPTION_STORAGE',
    'VEG_LAI_WATER_FACTOR',
    'CANOPY_CLOSURE',
    'CANOPY_RSMAX',
    'CANOPY_VPDMINFACTOR',
    'MTCLIM_TDAYCOEF',
    'MTCLIM_SOLAR_CONSTANT',
    'MTCLIM_SNOW_TCRIT',
    'MTCLIM_SNOW_TRATE',
    'MTCLIM_TBASE',
    'MTCLIM_ABASE',
    'MTCLIM_C',
    'MTCLIM_B0',
    'MTCLIM_B1',
    'MTCLIM_B2',
    'MTCLIM_RAIN_SCALAR',
    'MTCLIM_DIF_ALB',
    'MTCLIM_SC_INT',
    'MTCLIM_SC_SLOPE',
    'MTCLIM_SRADDT',
    'MTCLIM_SW_PREC_THRESH',
    'LAKE_TMELT',
    'LAKE_MAX_SURFACE',
    'LAKE_BETA',
    'LAKE_FRACMIN',
    'LAKE_FRACLIM',
    'LAKE_DM',
    'LAKE_SNOWCRIT',
    'LAKE_ZWATER',
    'LAKE_ZSNOW',
    'LAKE_RHOSNOW',
    'LAKE_CONDI',
    'LAKE_CONDS',
    'LAKE_LAMISW',
    'LAKE_LAMILW',
    'LAKE_LAMSSW',
    'LAKE_LAMSLW',
    'LAKE_LAMWSW',
    'LAKE_LAMWLW',
    'LAKE_A1',
    'LAKE_A2',
    'LAKE_QWTAU',
    'LAKE_MAX_ITER',
    'SVP_A',
    'SVP_B',
    'SVP_C',
    'CARBON_CATMCURRENT',
    'CARBON_SW2PAR',
    'PHOTO_OMEGA',
    'PHOTO_LAIMAX',
    'PHOTO_LAILIMIT',
    'PHOTO_LAIMIN',
    'PHOTO_EPAR',
    'PHOTO_FCMAX',
    'PHOTO_FCMIN',
    'PHOTO_ZENITHMIN',
    'PHOTO_ZENITHMINPAR',
    'PHOTO_ALBSOIPARMIN',
    'PHOTO_MINMAXETRANS',
    'PHOTO_MINSTOMCOND',
    'PHOTO_FCI1C3',
    'PHOTO_FCI1C4',
    'PHOTO_OX',
    'PHOTO_KC',
    'PHOTO_KO',
    'PHOTO_EC',
    'PHOTO_EO',
    'PHOTO_EV',
    'PHOTO_ER',
    'PHOTO_ALC3',
    'PHOTO_FRDC3',
    'PHOTO_EK',
    'PHOTO_ALC4',
    'PHOTO_FRDC4',
    'PHOTO_THETA',
    'PHOTO_FRLEAF',
    'PHOTO_FRGROWTH',
    'SRESP_E0_LT',
    'SRESP_T0_LT',
    'SRESP_WMINFM',
    'SRESP_WMAXFM',
    'SRESP_WOPTFM',
    'SRESP_RHSAT',
    'SRESP_RFACTOR',
    'SRESP_TAULITTER',
    'SRESP_TAUINTER',
    'SRESP_TAUSLOW',
    'SRESP_FAIR',
    'SRESP_FINTER',
    'SNOW_MAX_SURFACE_SWE',
    'SNOW_LIQUID_WATER_CAPACITY',
    'SNOW_NEW_SNOW_DENSITY',
    'SNOW_DENS_DMLIMIT',
    'SNOW_DENS_MAX_CHANGE',
    'SNOW_DENS_ETA0',
    'SNOW_DENS_C1',
    'SNOW_DENS_C2',
    'SNOW_DENS_C5',
    'SNOW_DENS_C6',
    'SNOW_DENS_F',
    'SNOW_MIN_SWQ_EB_THRES',
    'SNOW_A1',
    'SNOW_A2',
    'SNOW_L1',
    'SNOW_L2',
    'SNOW_NEW_SNOW_ALB',
    'SNOW_ALB_ACCUM_A',
    'SNOW_ALB_ACCUM_B',
    'SNOW_ALB_THAW_A',
    'SNOW_ALB_THAW_B',
    'SNOW_TRACESNOW',
    'SNOW_CONDUCT',
    'SNOW_MAX_SNOW_TEMP',
    'SNOW_MIN_RAIN_TEMP',
    'BLOWING_KA',
    'BLOWING_CSALT',
    'BLOWING_UTHRESH',
    'BLOWING_KIN_VIS',
    'BLOWING_MAX_ITER',
    'BLOWING_K',
    'BLOWING_SETTLING',
    'BLOWING_NUMINCS',
    'TREELINE_TEMPERATURE',
    'SNOW_DT',
    'SURF_DT',
    'SOIL_DT',
    'CANOPY_DT',
    'CANOPY_VP',
    'TOL_GRND',
    'TOL_OVER',
    'FROZEN_MAXITER',
    'NEWT_RAPH_MAXTRIAL',
    'NEWT_RAPH_TOLX',
    'NEWT_RAPH_TOLF',
    'NEWT_RAPH_R_MAX',
    'NEWT_RAPH_R_MIN',
    'NEWT_RAPH_RELAX1',
    'NEWT_RAPH_RELAX2',
    'NEWT_RAPH_RELAX3',
    'NEWT_RAPH_EPS2',
    'ROOT_BRENT_MAXTRIES',
    'ROOT_BRENT_MAXITER',
    'ROOT_BRENT_TSTEP',
    'ROOT_BRENT_T',
]
struct_anon_31._fields_ = [
    ('LAPSE_RATE', c_double),
    ('GAUGE_HEIGHT', c_double),
    ('WIND_SPEED_DEFAULT', c_double),
    ('WIND_SPEED_MIN', c_double),
    ('HUGE_RESIST', c_double),
    ('ALBEDO_BARE_SOIL', c_double),
    ('ALBEDO_H20_SURF', c_double),
    ('EMISS_GRND', c_double),
    ('EMISS_VEG', c_double),
    ('EMISS_ICE', c_double),
    ('EMISS_SNOW', c_double),
    ('EMISS_H2O', c_double),
    ('SOIL_RESID_MOIST', c_double),
    ('SOIL_SLAB_MOIST_FRACT', c_double),
    ('VEG_LAI_SNOW_MULTIPLIER', c_double),
    ('VEG_MIN_INTERCEPTION_STORAGE', c_double),
    ('VEG_LAI_WATER_FACTOR', c_double),
    ('CANOPY_CLOSURE', c_double),
    ('CANOPY_RSMAX', c_double),
    ('CANOPY_VPDMINFACTOR', c_double),
    ('MTCLIM_TDAYCOEF', c_double),
    ('MTCLIM_SOLAR_CONSTANT', c_double),
    ('MTCLIM_SNOW_TCRIT', c_double),
    ('MTCLIM_SNOW_TRATE', c_double),
    ('MTCLIM_TBASE', c_double),
    ('MTCLIM_ABASE', c_double),
    ('MTCLIM_C', c_double),
    ('MTCLIM_B0', c_double),
    ('MTCLIM_B1', c_double),
    ('MTCLIM_B2', c_double),
    ('MTCLIM_RAIN_SCALAR', c_double),
    ('MTCLIM_DIF_ALB', c_double),
    ('MTCLIM_SC_INT', c_double),
    ('MTCLIM_SC_SLOPE', c_double),
    ('MTCLIM_SRADDT', c_double),
    ('MTCLIM_SW_PREC_THRESH', c_double),
    ('LAKE_TMELT', c_double),
    ('LAKE_MAX_SURFACE', c_double),
    ('LAKE_BETA', c_double),
    ('LAKE_FRACMIN', c_double),
    ('LAKE_FRACLIM', c_double),
    ('LAKE_DM', c_double),
    ('LAKE_SNOWCRIT', c_double),
    ('LAKE_ZWATER', c_double),
    ('LAKE_ZSNOW', c_double),
    ('LAKE_RHOSNOW', c_double),
    ('LAKE_CONDI', c_double),
    ('LAKE_CONDS', c_double),
    ('LAKE_LAMISW', c_double),
    ('LAKE_LAMILW', c_double),
    ('LAKE_LAMSSW', c_double),
    ('LAKE_LAMSLW', c_double),
    ('LAKE_LAMWSW', c_double),
    ('LAKE_LAMWLW', c_double),
    ('LAKE_A1', c_double),
    ('LAKE_A2', c_double),
    ('LAKE_QWTAU', c_double),
    ('LAKE_MAX_ITER', c_int),
    ('SVP_A', c_double),
    ('SVP_B', c_double),
    ('SVP_C', c_double),
    ('CARBON_CATMCURRENT', c_double),
    ('CARBON_SW2PAR', c_double),
    ('PHOTO_OMEGA', c_double),
    ('PHOTO_LAIMAX', c_double),
    ('PHOTO_LAILIMIT', c_double),
    ('PHOTO_LAIMIN', c_double),
    ('PHOTO_EPAR', c_double),
    ('PHOTO_FCMAX', c_double),
    ('PHOTO_FCMIN', c_double),
    ('PHOTO_ZENITHMIN', c_double),
    ('PHOTO_ZENITHMINPAR', c_double),
    ('PHOTO_ALBSOIPARMIN', c_double),
    ('PHOTO_MINMAXETRANS', c_double),
    ('PHOTO_MINSTOMCOND', c_double),
    ('PHOTO_FCI1C3', c_double),
    ('PHOTO_FCI1C4', c_double),
    ('PHOTO_OX', c_double),
    ('PHOTO_KC', c_double),
    ('PHOTO_KO', c_double),
    ('PHOTO_EC', c_double),
    ('PHOTO_EO', c_double),
    ('PHOTO_EV', c_double),
    ('PHOTO_ER', c_double),
    ('PHOTO_ALC3', c_double),
    ('PHOTO_FRDC3', c_double),
    ('PHOTO_EK', c_double),
    ('PHOTO_ALC4', c_double),
    ('PHOTO_FRDC4', c_double),
    ('PHOTO_THETA', c_double),
    ('PHOTO_FRLEAF', c_double),
    ('PHOTO_FRGROWTH', c_double),
    ('SRESP_E0_LT', c_double),
    ('SRESP_T0_LT', c_double),
    ('SRESP_WMINFM', c_double),
    ('SRESP_WMAXFM', c_double),
    ('SRESP_WOPTFM', c_double),
    ('SRESP_RHSAT', c_double),
    ('SRESP_RFACTOR', c_double),
    ('SRESP_TAULITTER', c_double),
    ('SRESP_TAUINTER', c_double),
    ('SRESP_TAUSLOW', c_double),
    ('SRESP_FAIR', c_double),
    ('SRESP_FINTER', c_double),
    ('SNOW_MAX_SURFACE_SWE', c_double),
    ('SNOW_LIQUID_WATER_CAPACITY', c_double),
    ('SNOW_NEW_SNOW_DENSITY', c_double),
    ('SNOW_DENS_DMLIMIT', c_double),
    ('SNOW_DENS_MAX_CHANGE', c_double),
    ('SNOW_DENS_ETA0', c_double),
    ('SNOW_DENS_C1', c_double),
    ('SNOW_DENS_C2', c_double),
    ('SNOW_DENS_C5', c_double),
    ('SNOW_DENS_C6', c_double),
    ('SNOW_DENS_F', c_double),
    ('SNOW_MIN_SWQ_EB_THRES', c_double),
    ('SNOW_A1', c_double),
    ('SNOW_A2', c_double),
    ('SNOW_L1', c_double),
    ('SNOW_L2', c_double),
    ('SNOW_NEW_SNOW_ALB', c_double),
    ('SNOW_ALB_ACCUM_A', c_double),
    ('SNOW_ALB_ACCUM_B', c_double),
    ('SNOW_ALB_THAW_A', c_double),
    ('SNOW_ALB_THAW_B', c_double),
    ('SNOW_TRACESNOW', c_double),
    ('SNOW_CONDUCT', c_double),
    ('SNOW_MAX_SNOW_TEMP', c_double),
    ('SNOW_MIN_RAIN_TEMP', c_double),
    ('BLOWING_KA', c_double),
    ('BLOWING_CSALT', c_double),
    ('BLOWING_UTHRESH', c_double),
    ('BLOWING_KIN_VIS', c_double),
    ('BLOWING_MAX_ITER', c_int),
    ('BLOWING_K', c_int),
    ('BLOWING_SETTLING', c_double),
    ('BLOWING_NUMINCS', c_int),
    ('TREELINE_TEMPERATURE', c_double),
    ('SNOW_DT', c_double),
    ('SURF_DT', c_double),
    ('SOIL_DT', c_double),
    ('CANOPY_DT', c_double),
    ('CANOPY_VP', c_double),
    ('TOL_GRND', c_double),
    ('TOL_OVER', c_double),
    ('FROZEN_MAXITER', c_int),
    ('NEWT_RAPH_MAXTRIAL', c_int),
    ('NEWT_RAPH_TOLX', c_double),
    ('NEWT_RAPH_TOLF', c_double),
    ('NEWT_RAPH_R_MAX', c_double),
    ('NEWT_RAPH_R_MIN', c_double),
    ('NEWT_RAPH_RELAX1', c_double),
    ('NEWT_RAPH_RELAX2', c_double),
    ('NEWT_RAPH_RELAX3', c_double),
    ('NEWT_RAPH_EPS2', c_double),
    ('ROOT_BRENT_MAXTRIES', c_int),
    ('ROOT_BRENT_MAXITER', c_int),
    ('ROOT_BRENT_TSTEP', c_double),
    ('ROOT_BRENT_T', c_double),
]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 982
parameters_struct = struct_anon_31

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 1057


class struct_anon_32(Structure):
    pass

struct_anon_32.__slots__ = [
    'FS_ACTIVE',
    'Ds',
    'Dsmax',
    'Ksat',
    'Wcr',
    'Wpwp',
    'Ws',
    'AlbedoPar',
    'alpha',
    'annual_prec',
    'avg_temp',
    'avgJulyAirTemp',
    'b_infilt',
    'beta',
    'bubble',
    'bubble_node',
    'bulk_density',
    'bulk_dens_min',
    'bulk_dens_org',
    'c',
    'depth',
    'dp',
    'dz_node',
    'Zsum_node',
    'expt',
    'expt_node',
    'frost_fract',
    'frost_slope',
    'gamma',
    'init_moist',
    'max_infil',
    'max_moist',
    'max_moist_node',
    'max_snow_distrib_slope',
    'phi_s',
    'porosity',
    'quartz',
    'organic',
    'resid_moist',
    'rough',
    'snow_rough',
    'soil_density',
    'soil_dens_min',
    'soil_dens_org',
    'BandElev',
    'AreaFract',
    'Pfactor',
    'Tfactor',
    'AboveTreeLine',
    'elevation',
    'lat',
    'lng',
    'cell_area',
    'time_zone_lng',
    'gridcel',
    'zwtvmoist_zwt',
    'zwtvmoist_moist',
    'slope',
    'aspect',
    'ehoriz',
    'whoriz',
]
struct_anon_32._fields_ = [
    ('FS_ACTIVE', c_uint8),
    ('Ds', c_double),
    ('Dsmax', c_double),
    ('Ksat', c_double * 3),
    ('Wcr', c_double * 3),
    ('Wpwp', c_double * 3),
    ('Ws', c_double),
    ('AlbedoPar', c_double),
    ('alpha', c_double * 50),
    ('annual_prec', c_double),
    ('avg_temp', c_double),
    ('avgJulyAirTemp', c_double),
    ('b_infilt', c_double),
    ('beta', c_double * 50),
    ('bubble', c_double * 3),
    ('bubble_node', c_double * 50),
    ('bulk_density', c_double * 3),
    ('bulk_dens_min', c_double * 3),
    ('bulk_dens_org', c_double * 3),
    ('c', c_double),
    ('depth', c_double * 3),
    ('dp', c_double),
    ('dz_node', c_double * 50),
    ('Zsum_node', c_double * 50),
    ('expt', c_double * 3),
    ('expt_node', c_double * 50),
    ('frost_fract', c_double * 10),
    ('frost_slope', c_double),
    ('gamma', c_double * 50),
    ('init_moist', c_double * 3),
    ('max_infil', c_double),
    ('max_moist', c_double * 3),
    ('max_moist_node', c_double * 50),
    ('max_snow_distrib_slope', c_double),
    ('phi_s', c_double * 3),
    ('porosity', c_double * 3),
    ('quartz', c_double * 3),
    ('organic', c_double * 3),
    ('resid_moist', c_double * 3),
    ('rough', c_double),
    ('snow_rough', c_double),
    ('soil_density', c_double * 3),
    ('soil_dens_min', c_double * 3),
    ('soil_dens_org', c_double * 3),
    ('BandElev', POINTER(c_double)),
    ('AreaFract', POINTER(c_double)),
    ('Pfactor', POINTER(c_double)),
    ('Tfactor', POINTER(c_double)),
    ('AboveTreeLine', POINTER(c_uint8)),
    ('elevation', c_double),
    ('lat', c_double),
    ('lng', c_double),
    ('cell_area', c_double),
    ('time_zone_lng', c_double),
    ('gridcel', c_uint),
    ('zwtvmoist_zwt', (c_double * 11) * (3 + 2)),
    ('zwtvmoist_moist', (c_double * 11) * (3 + 2)),
    ('slope', c_double),
    ('aspect', c_double),
    ('ehoriz', c_double),
    ('whoriz', c_double),
]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 1057
soil_con_struct = struct_anon_32

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 1076


class struct_anon_33(Structure):
    pass

struct_anon_33.__slots__ = [
    'Cv',
    'Cv_sum',
    'root',
    'zone_depth',
    'zone_fract',
    'veg_class',
    'vegetat_type_num',
    'sigma_slope',
    'lag_one',
    'fetch',
    'LAKE',
    'CanopLayerBnd',
]
struct_anon_33._fields_ = [
    ('Cv', c_double),
    ('Cv_sum', c_double),
    ('root', c_double * 3),
    ('zone_depth', POINTER(c_double)),
    ('zone_fract', POINTER(c_double)),
    ('veg_class', c_int),
    ('vegetat_type_num', c_size_t),
    ('sigma_slope', c_double),
    ('lag_one', c_double),
    ('fetch', c_double),
    ('LAKE', c_int),
    ('CanopLayerBnd', POINTER(c_double)),
]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 1076
veg_con_struct = struct_anon_33

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 1118


class struct_anon_34(Structure):
    pass

struct_anon_34.__slots__ = [
    'overstory',
    'LAI',
    'vegcover',
    'Wdmax',
    'albedo',
    'displacement',
    'emissivity',
    'NVegLibTypes',
    'rad_atten',
    'rarc',
    'rmin',
    'roughness',
    'trunk_ratio',
    'wind_atten',
    'wind_h',
    'RGL',
    'veg_class',
    'Ctype',
    'MaxCarboxRate',
    'MaxETransport',
    'CO2Specificity',
    'LightUseEff',
    'NscaleFlag',
    'Wnpp_inhib',
    'NPPfactor_sat',
]
struct_anon_34._fields_ = [
    ('overstory', c_uint8),
    ('LAI', c_double * 12),
    ('vegcover', c_double * 12),
    ('Wdmax', c_double * 12),
    ('albedo', c_double * 12),
    ('displacement', c_double * 12),
    ('emissivity', c_double * 12),
    ('NVegLibTypes', c_size_t),
    ('rad_atten', c_double),
    ('rarc', c_double),
    ('rmin', c_double),
    ('roughness', c_double * 12),
    ('trunk_ratio', c_double),
    ('wind_atten', c_double),
    ('wind_h', c_double),
    ('RGL', c_double),
    ('veg_class', c_uint),
    ('Ctype', c_char),
    ('MaxCarboxRate', c_double),
    ('MaxETransport', c_double),
    ('CO2Specificity', c_double),
    ('LightUseEff', c_double),
    ('NscaleFlag', c_uint8),
    ('Wnpp_inhib', c_double),
    ('NPPfactor_sat', c_double),
]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 1118
veg_lib_struct = struct_anon_34

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 1129


class struct_anon_35(Structure):
    pass

struct_anon_35.__slots__ = [
    'albedo',
    'LAI',
    'vegcover',
]
struct_anon_35._fields_ = [
    ('albedo', POINTER(c_double)),
    ('LAI', POINTER(c_double)),
    ('vegcover', POINTER(c_double)),
]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 1129
veg_hist_struct = struct_anon_35

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 1161


class struct_anon_36(Structure):
    pass

struct_anon_36.__slots__ = [
    'air_temp',
    'Catm',
    'channel_in',
    'coszen',
    'density',
    'fdir',
    'longwave',
    'out_prec',
    'out_rain',
    'out_snow',
    'par',
    'prec',
    'pressure',
    'shortwave',
    'snowflag',
    'tskc',
    'vp',
    'vpd',
    'wind',
]
struct_anon_36._fields_ = [
    ('air_temp', POINTER(c_double)),
    ('Catm', POINTER(c_double)),
    ('channel_in', POINTER(c_double)),
    ('coszen', POINTER(c_double)),
    ('density', POINTER(c_double)),
    ('fdir', POINTER(c_double)),
    ('longwave', POINTER(c_double)),
    ('out_prec', c_double),
    ('out_rain', c_double),
    ('out_snow', c_double),
    ('par', POINTER(c_double)),
    ('prec', POINTER(c_double)),
    ('pressure', POINTER(c_double)),
    ('shortwave', POINTER(c_double)),
    ('snowflag', POINTER(c_uint8)),
    ('tskc', POINTER(c_double)),
    ('vp', POINTER(c_double)),
    ('vpd', POINTER(c_double)),
    ('wind', POINTER(c_double)),
]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 1161
atmos_data_struct = struct_anon_36

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 1173


class struct_anon_37(Structure):
    pass

struct_anon_37.__slots__ = [
    'day',
    'day_in_year',
    'month',
    'year',
    'dayseconds',
]
struct_anon_37._fields_ = [
    ('day', c_uint),
    ('day_in_year', c_uint),
    ('month', c_uint),
    ('year', c_int),
    ('dayseconds', c_uint),
]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 1173
dmy_struct = struct_anon_37

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 1192


class struct_anon_38(Structure):
    pass

struct_anon_38.__slots__ = [
    'Cs',
    'T',
    'bare_evap_frac',
    'evap',
    'ice',
    'kappa',
    'moist',
    'phi',
    'zwt',
]
struct_anon_38._fields_ = [
    ('Cs', c_double),
    ('T', c_double),
    ('bare_evap_frac', c_double),
    ('evap', c_double),
    ('ice', c_double * 10),
    ('kappa', c_double),
    ('moist', c_double),
    ('phi', c_double),
    ('zwt', c_double),
]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 1192
layer_data_struct = struct_anon_38

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 1227


class struct_anon_39(Structure):
    pass

struct_anon_39.__slots__ = [
    'aero_resist',
    'asat',
    'baseflow',
    'CLitter',
    'CInter',
    'CSlow',
    'inflow',
    'pot_evap',
    'runoff',
    'layer',
    'RhLitter',
    'RhLitter2Atm',
    'RhInter',
    'RhSlow',
    'RhTot',
    'rootmoist',
    'wetness',
    'zwt',
    'zwt_lumped',
]
struct_anon_39._fields_ = [
    ('aero_resist', c_double * 2),
    ('asat', c_double),
    ('baseflow', c_double),
    ('CLitter', c_double),
    ('CInter', c_double),
    ('CSlow', c_double),
    ('inflow', c_double),
    ('pot_evap', c_double * 0),
    ('runoff', c_double),
    ('layer', layer_data_struct * 3),
    ('RhLitter', c_double),
    ('RhLitter2Atm', c_double),
    ('RhInter', c_double),
    ('RhSlow', c_double),
    ('RhTot', c_double),
    ('rootmoist', c_double),
    ('wetness', c_double),
    ('zwt', c_double),
    ('zwt_lumped', c_double),
]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 1227
cell_data_struct = struct_anon_39

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 1302


class struct_anon_40(Structure):
    pass

struct_anon_40.__slots__ = [
    'AlbedoLake',
    'AlbedoOver',
    'AlbedoUnder',
    'Cs',
    'Cs_node',
    'fdepth',
    'frozen',
    'ice',
    'kappa',
    'kappa_node',
    'moist',
    'Nfrost',
    'Nthaw',
    'T',
    'T_fbflag',
    'T_fbcount',
    'T1_index',
    'Tcanopy',
    'Tcanopy_fbflag',
    'Tcanopy_fbcount',
    'tdepth',
    'Tfoliage',
    'Tfoliage_fbflag',
    'Tfoliage_fbcount',
    'Tsurf',
    'Tsurf_fbflag',
    'Tsurf_fbcount',
    'unfrozen',
    'advected_sensible',
    'advection',
    'AtmosError',
    'AtmosLatent',
    'AtmosLatentSub',
    'AtmosSensible',
    'canopy_advection',
    'canopy_latent',
    'canopy_latent_sub',
    'canopy_refreeze',
    'canopy_sensible',
    'deltaCC',
    'deltaH',
    'error',
    'fusion',
    'grnd_flux',
    'latent',
    'latent_sub',
    'longwave',
    'LongOverIn',
    'LongUnderIn',
    'LongUnderOut',
    'melt_energy',
    'NetLongAtmos',
    'NetLongOver',
    'NetLongUnder',
    'NetShortAtmos',
    'NetShortGrnd',
    'NetShortOver',
    'NetShortUnder',
    'out_long_canopy',
    'out_long_surface',
    'refreeze_energy',
    'sensible',
    'shortwave',
    'ShortOverIn',
    'ShortUnderIn',
    'snow_flux',
]
struct_anon_40._fields_ = [
    ('AlbedoLake', c_double),
    ('AlbedoOver', c_double),
    ('AlbedoUnder', c_double),
    ('Cs', c_double * 2),
    ('Cs_node', c_double * 50),
    ('fdepth', c_double * 3),
    ('frozen', c_uint8),
    ('ice', c_double * 50),
    ('kappa', c_double * 2),
    ('kappa_node', c_double * 50),
    ('moist', c_double * 50),
    ('Nfrost', c_size_t),
    ('Nthaw', c_size_t),
    ('T', c_double * 50),
    ('T_fbflag', c_uint8 * 50),
    ('T_fbcount', c_uint * 50),
    ('T1_index', c_int),
    ('Tcanopy', c_double),
    ('Tcanopy_fbflag', c_uint8),
    ('Tcanopy_fbcount', c_uint),
    ('tdepth', c_double * 3),
    ('Tfoliage', c_double),
    ('Tfoliage_fbflag', c_uint8),
    ('Tfoliage_fbcount', c_uint),
    ('Tsurf', c_double),
    ('Tsurf_fbflag', c_uint8),
    ('Tsurf_fbcount', c_uint),
    ('unfrozen', c_double),
    ('advected_sensible', c_double),
    ('advection', c_double),
    ('AtmosError', c_double),
    ('AtmosLatent', c_double),
    ('AtmosLatentSub', c_double),
    ('AtmosSensible', c_double),
    ('canopy_advection', c_double),
    ('canopy_latent', c_double),
    ('canopy_latent_sub', c_double),
    ('canopy_refreeze', c_double),
    ('canopy_sensible', c_double),
    ('deltaCC', c_double),
    ('deltaH', c_double),
    ('error', c_double),
    ('fusion', c_double),
    ('grnd_flux', c_double),
    ('latent', c_double),
    ('latent_sub', c_double),
    ('longwave', c_double),
    ('LongOverIn', c_double),
    ('LongUnderIn', c_double),
    ('LongUnderOut', c_double),
    ('melt_energy', c_double),
    ('NetLongAtmos', c_double),
    ('NetLongOver', c_double),
    ('NetLongUnder', c_double),
    ('NetShortAtmos', c_double),
    ('NetShortGrnd', c_double),
    ('NetShortOver', c_double),
    ('NetShortUnder', c_double),
    ('out_long_canopy', c_double),
    ('out_long_surface', c_double),
    ('refreeze_energy', c_double),
    ('sensible', c_double),
    ('shortwave', c_double),
    ('ShortOverIn', c_double),
    ('ShortUnderIn', c_double),
    ('snow_flux', c_double),
]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 1302
energy_bal_struct = struct_anon_40

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 1334


class struct_anon_41(Structure):
    pass

struct_anon_41.__slots__ = [
    'albedo',
    'canopyevap',
    'LAI',
    'throughfall',
    'vegcover',
    'Wdew',
    'Wdmax',
    'NscaleFactor',
    'aPARLayer',
    'CiLayer',
    'rsLayer',
    'aPAR',
    'Ci',
    'rc',
    'NPPfactor',
    'GPP',
    'Rphoto',
    'Rdark',
    'Rmaint',
    'Rgrowth',
    'Raut',
    'NPP',
    'Litterfall',
    'AnnualNPP',
    'AnnualNPPPrev',
]
struct_anon_41._fields_ = [
    ('albedo', c_double),
    ('canopyevap', c_double),
    ('LAI', c_double),
    ('throughfall', c_double),
    ('vegcover', c_double),
    ('Wdew', c_double),
    ('Wdmax', c_double),
    ('NscaleFactor', POINTER(c_double)),
    ('aPARLayer', POINTER(c_double)),
    ('CiLayer', POINTER(c_double)),
    ('rsLayer', POINTER(c_double)),
    ('aPAR', c_double),
    ('Ci', c_double),
    ('rc', c_double),
    ('NPPfactor', c_double),
    ('GPP', c_double),
    ('Rphoto', c_double),
    ('Rdark', c_double),
    ('Rmaint', c_double),
    ('Rgrowth', c_double),
    ('Raut', c_double),
    ('NPP', c_double),
    ('Litterfall', c_double),
    ('AnnualNPP', c_double),
    ('AnnualNPPPrev', c_double),
]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 1334
veg_var_struct = struct_anon_41

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 1380


class struct_anon_42(Structure):
    pass

struct_anon_42.__slots__ = [
    'albedo',
    'canopy_albedo',
    'coldcontent',
    'coverage',
    'density',
    'depth',
    'last_snow',
    'max_snow_depth',
    'MELTING',
    'pack_temp',
    'pack_water',
    'snow',
    'snow_canopy',
    'store_coverage',
    'store_snow',
    'store_swq',
    'surf_temp',
    'surf_temp_fbcount',
    'surf_temp_fbflag',
    'surf_water',
    'swq',
    'snow_distrib_slope',
    'tmp_int_storage',
    'blowing_flux',
    'canopy_vapor_flux',
    'mass_error',
    'melt',
    'Qnet',
    'surface_flux',
    'transport',
    'vapor_flux',
]
struct_anon_42._fields_ = [
    ('albedo', c_double),
    ('canopy_albedo', c_double),
    ('coldcontent', c_double),
    ('coverage', c_double),
    ('density', c_double),
    ('depth', c_double),
    ('last_snow', c_uint),
    ('max_snow_depth', c_double),
    ('MELTING', c_uint8),
    ('pack_temp', c_double),
    ('pack_water', c_double),
    ('snow', c_uint8),
    ('snow_canopy', c_double),
    ('store_coverage', c_double),
    ('store_snow', c_uint8),
    ('store_swq', c_double),
    ('surf_temp', c_double),
    ('surf_temp_fbcount', c_uint),
    ('surf_temp_fbflag', c_uint8),
    ('surf_water', c_double),
    ('swq', c_double),
    ('snow_distrib_slope', c_double),
    ('tmp_int_storage', c_double),
    ('blowing_flux', c_double),
    ('canopy_vapor_flux', c_double),
    ('mass_error', c_double),
    ('melt', c_double),
    ('Qnet', c_double),
    ('surface_flux', c_double),
    ('transport', c_double),
    ('vapor_flux', c_double),
]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 1380
snow_data_struct = struct_anon_42

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 1403


class struct_anon_43(Structure):
    pass

struct_anon_43.__slots__ = [
    'numnod',
    'z',
    'basin',
    'Cl',
    'b',
    'maxdepth',
    'mindepth',
    'maxvolume',
    'minvolume',
    'bpercent',
    'rpercent',
    'wfrac',
    'depth_in',
    'lake_idx',
]
struct_anon_43._fields_ = [
    ('numnod', c_size_t),
    ('z', c_double * (20 + 1)),
    ('basin', c_double * (20 + 1)),
    ('Cl', c_double * (20 + 1)),
    ('b', c_double),
    ('maxdepth', c_double),
    ('mindepth', c_double),
    ('maxvolume', c_double),
    ('minvolume', c_double),
    ('bpercent', c_double),
    ('rpercent', c_double),
    ('wfrac', c_double),
    ('depth_in', c_double),
    ('lake_idx', c_int),
]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 1403
lake_con_struct = struct_anon_43

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 1456


class struct_anon_44(Structure):
    pass

struct_anon_44.__slots__ = [
    'activenod',
    'dz',
    'surfdz',
    'ldepth',
    'surface',
    'sarea',
    'sarea_save',
    'volume',
    'volume_save',
    'temp',
    'tempavg',
    'areai',
    'new_ice_area',
    'ice_water_eq',
    'hice',
    'tempi',
    'swe',
    'swe_save',
    'surf_temp',
    'pack_temp',
    'coldcontent',
    'surf_water',
    'pack_water',
    'SAlbedo',
    'sdepth',
    'aero_resist',
    'density',
    'baseflow_in',
    'baseflow_out',
    'channel_in',
    'evapw',
    'ice_throughfall',
    'prec',
    'recharge',
    'runoff_in',
    'runoff_out',
    'snowmlt',
    'vapor_flux',
    'snow',
    'energy',
    'soil',
]
struct_anon_44._fields_ = [
    ('activenod', c_uint),
    ('dz', c_double),
    ('surfdz', c_double),
    ('ldepth', c_double),
    ('surface', c_double * (20 + 1)),
    ('sarea', c_double),
    ('sarea_save', c_double),
    ('volume', c_double),
    ('volume_save', c_double),
    ('temp', c_double * 20),
    ('tempavg', c_double),
    ('areai', c_double),
    ('new_ice_area', c_double),
    ('ice_water_eq', c_double),
    ('hice', c_double),
    ('tempi', c_double),
    ('swe', c_double),
    ('swe_save', c_double),
    ('surf_temp', c_double),
    ('pack_temp', c_double),
    ('coldcontent', c_double),
    ('surf_water', c_double),
    ('pack_water', c_double),
    ('SAlbedo', c_double),
    ('sdepth', c_double),
    ('aero_resist', c_double),
    ('density', c_double * 20),
    ('baseflow_in', c_double),
    ('baseflow_out', c_double),
    ('channel_in', c_double),
    ('evapw', c_double),
    ('ice_throughfall', c_double),
    ('prec', c_double),
    ('recharge', c_double),
    ('runoff_in', c_double),
    ('runoff_out', c_double),
    ('snowmlt', c_double),
    ('vapor_flux', c_double),
    ('snow', snow_data_struct),
    ('energy', energy_bal_struct),
    ('soil', cell_data_struct),
]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 1456
lake_var_struct = struct_anon_44

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 1468


class struct_anon_45(Structure):
    pass

struct_anon_45.__slots__ = [
    'cell',
    'energy',
    'lake_var',
    'snow',
    'veg_var',
]
struct_anon_45._fields_ = [
    ('cell', POINTER(POINTER(cell_data_struct))),
    ('energy', POINTER(POINTER(energy_bal_struct))),
    ('lake_var', lake_var_struct),
    ('snow', POINTER(POINTER(snow_data_struct))),
    ('veg_var', POINTER(POINTER(veg_var_struct))),
]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 1468
all_vars_struct = struct_anon_45

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 1479


class struct_anon_46(Structure):
    pass

struct_anon_46.__slots__ = [
    'total_soil_moist',
    'surfstor',
    'swe',
    'wdew',
]
struct_anon_46._fields_ = [
    ('total_soil_moist', c_double),
    ('surfstor', c_double),
    ('swe', c_double),
    ('wdew', c_double),
]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 1479
save_data_struct = struct_anon_46

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 1505


class struct_anon_47(Structure):
    pass

struct_anon_47.__slots__ = [
    'varname',
    'write',
    'format',
    'type',
    'mult',
    'aggtype',
    'nelem',
    'data',
    'aggdata',
]
struct_anon_47._fields_ = [
    ('varname', c_char * 20),
    ('write', c_uint8),
    ('format', c_char * 10),
    ('type', c_uint),
    ('mult', c_double),
    ('aggtype', c_uint),
    ('nelem', c_uint),
    ('data', POINTER(c_double)),
    ('aggdata', POINTER(c_double)),
]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 1505
out_data_struct = struct_anon_47

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 1519


class struct_anon_48(Structure):
    pass

struct_anon_48.__slots__ = [
    'prefix',
    'filename',
    'fh',
    'nvars',
    'varid',
]
struct_anon_48._fields_ = [
    ('prefix', c_char * 20),
    ('filename', c_char * 2048),
    ('fh', POINTER(FILE)),
    ('nvars', c_size_t),
    ('varid', POINTER(c_uint)),
]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 1519
out_data_file_struct = struct_anon_48

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 1537


class struct_anon_49(Structure):
    pass

struct_anon_49.__slots__ = [
    'atmos',
    'dt',
    'energy',
    'filep',
    'rec',
    'out_data',
    'out_data_files',
    'snow',
    'soil_con',
    'veg_con',
    'veg_var',
]
struct_anon_49._fields_ = [
    ('atmos', POINTER(atmos_data_struct)),
    ('dt', c_double),
    ('energy', POINTER(energy_bal_struct)),
    ('filep', filep_struct),
    ('rec', c_size_t),
    ('out_data', POINTER(out_data_struct)),
    ('out_data_files', POINTER(out_data_file_struct)),
    ('snow', POINTER(snow_data_struct)),
    ('soil_con', soil_con_struct),
    ('veg_con', POINTER(veg_con_struct)),
    ('veg_var', POINTER(veg_var_struct)),
]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 1537
Error_struct = struct_anon_49

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 30
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'advect_carbon_storage'):
    advect_carbon_storage = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].advect_carbon_storage
    advect_carbon_storage.restype = None
    advect_carbon_storage.argtypes = [
        c_double,
        c_double,
        POINTER(lake_var_struct),
        POINTER(cell_data_struct)]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 32
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'advect_snow_storage'):
    advect_snow_storage = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].advect_snow_storage
    advect_snow_storage.restype = None
    advect_snow_storage.argtypes = [
        c_double, c_double, c_double, POINTER(snow_data_struct)]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 33
if hasattr(_libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
           'advect_soil_veg_storage'):
    advect_soil_veg_storage = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].advect_soil_veg_storage
    advect_soil_veg_storage.restype = None
    advect_soil_veg_storage.argtypes = [
        c_double,
        c_double,
        c_double,
        POINTER(c_double),
        POINTER(soil_con_struct),
        POINTER(veg_con_struct),
        POINTER(cell_data_struct),
        POINTER(veg_var_struct),
        lake_con_struct]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 37
if hasattr(_libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
           'advected_sensible_heat'):
    advected_sensible_heat = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].advected_sensible_heat
    advected_sensible_heat.restype = c_double
    advected_sensible_heat.argtypes = [
        c_double, c_double, c_double, c_double, c_double]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 38
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'alblake'):
    alblake = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].alblake
    alblake.restype = None
    alblake.argtypes = [
        c_double,
        c_double,
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        c_double,
        c_double,
        c_double,
        POINTER(c_uint),
        c_double,
        POINTER(c_uint8),
        c_uint,
        c_double]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 41
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'arno_evap'):
    arno_evap = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].arno_evap
    arno_evap.restype = c_double
    arno_evap.argtypes = [
        POINTER(layer_data_struct),
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        POINTER(c_double)]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 43
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'calc_atmos_energy_bal'):
    calc_atmos_energy_bal = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].calc_atmos_energy_bal
    calc_atmos_energy_bal.restype = c_double
    calc_atmos_energy_bal.argtypes = [
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_uint8),
        POINTER(c_uint)]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 47
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'calc_density'):
    calc_density = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].calc_density
    calc_density.restype = c_double
    calc_density.argtypes = [c_double]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 48
if hasattr(_libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
           'calc_energy_balance_error'):
    calc_energy_balance_error = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].calc_energy_balance_error
    calc_energy_balance_error.restype = c_double
    calc_energy_balance_error.argtypes = [
        c_int, c_double, c_double, c_double, c_double, c_double]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 49
if hasattr(_libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
           'calc_latent_heat_of_sublimation'):
    calc_latent_heat_of_sublimation = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].calc_latent_heat_of_sublimation
    calc_latent_heat_of_sublimation.restype = c_double
    calc_latent_heat_of_sublimation.argtypes = [c_double]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 50
if hasattr(_libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
           'calc_latent_heat_of_vaporization'):
    calc_latent_heat_of_vaporization = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].calc_latent_heat_of_vaporization
    calc_latent_heat_of_vaporization.restype = c_double
    calc_latent_heat_of_vaporization.argtypes = [c_double]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 51
if hasattr(_libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
           'calc_layer_average_thermal_props'):
    calc_layer_average_thermal_props = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].calc_layer_average_thermal_props
    calc_layer_average_thermal_props.restype = c_int
    calc_layer_average_thermal_props.argtypes = [
        POINTER(energy_bal_struct),
        POINTER(layer_data_struct),
        POINTER(soil_con_struct),
        c_size_t,
        POINTER(c_double)]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 53
if hasattr(_libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
           'calc_outgoing_longwave'):
    calc_outgoing_longwave = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].calc_outgoing_longwave
    calc_outgoing_longwave.restype = c_double
    calc_outgoing_longwave.argtypes = [c_double, c_double]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 54
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'calc_scale_height'):
    calc_scale_height = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].calc_scale_height
    calc_scale_height.restype = c_double
    calc_scale_height.argtypes = [c_double, c_double]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 55
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'calc_sensible_heat'):
    calc_sensible_heat = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].calc_sensible_heat
    calc_sensible_heat.restype = c_double
    calc_sensible_heat.argtypes = [c_double, c_double, c_double, c_double]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 57
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'calc_Nscale_factors'):
    calc_Nscale_factors = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].calc_Nscale_factors
    calc_Nscale_factors.restype = None
    calc_Nscale_factors.argtypes = [
        c_char,
        POINTER(c_double),
        c_double,
        c_double,
        c_double,
        c_double,
        c_uint,
        POINTER(c_double)]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 59
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'calc_rainonly'):
    calc_rainonly = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].calc_rainonly
    calc_rainonly.restype = c_double
    calc_rainonly.argtypes = [c_double, c_double, c_double, c_double]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 60
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'calc_rc'):
    calc_rc = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].calc_rc
    calc_rc.restype = c_double
    calc_rc.argtypes = [c_double, c_double, c_double,
                        c_double, c_double, c_double, c_double, c_char]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 61
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'calc_rc_ps'):
    calc_rc_ps = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].calc_rc_ps
    calc_rc_ps.restype = None
    calc_rc_ps.argtypes = [
        c_char,
        c_double,
        c_double,
        c_double,
        POINTER(c_double),
        c_double,
        c_double,
        POINTER(c_double),
        c_double,
        c_double,
        POINTER(c_double),
        c_double,
        c_double,
        c_double,
        POINTER(c_double),
        POINTER(c_double)]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 64
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'calc_snow_coverage'):
    calc_snow_coverage = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].calc_snow_coverage
    calc_snow_coverage.restype = c_double
    calc_snow_coverage.argtypes = [
        POINTER(c_uint8),
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        POINTER(c_double),
        c_double,
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double)]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 67
if hasattr(_libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
           'calc_soil_thermal_fluxes'):
    calc_soil_thermal_fluxes = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].calc_soil_thermal_fluxes
    calc_soil_thermal_fluxes.restype = c_int
    calc_soil_thermal_fluxes.argtypes = [
        c_int,
        POINTER(c_double),
        POINTER(c_double),
        String,
        POINTER(c_uint),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        c_int,
        c_int,
        c_int]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 71
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'calc_surf_energy_bal'):
    calc_surf_energy_bal = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].calc_surf_energy_bal
    calc_surf_energy_bal.restype = c_double
    calc_surf_energy_bal.argtypes = [
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        c_double,
        POINTER(c_double),
        POINTER(c_double),
        c_double,
        POINTER(c_double),
        POINTER(c_double),
        c_int,
        c_int,
        c_size_t,
        c_size_t,
        c_double,
        c_size_t,
        c_uint,
        c_int,
        c_int,
        c_uint,
        POINTER(c_double),
        POINTER(c_double),
        POINTER(atmos_data_struct),
        POINTER(dmy_struct),
        POINTER(energy_bal_struct),
        POINTER(layer_data_struct),
        POINTER(snow_data_struct),
        POINTER(soil_con_struct),
        POINTER(veg_var_struct)]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 83
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'calc_veg_displacement'):
    calc_veg_displacement = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].calc_veg_displacement
    calc_veg_displacement.restype = c_double
    calc_veg_displacement.argtypes = [c_double]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 84
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'calc_veg_height'):
    calc_veg_height = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].calc_veg_height
    calc_veg_height.restype = c_double
    calc_veg_height.argtypes = [c_double]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 85
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'calc_veg_roughness'):
    calc_veg_roughness = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].calc_veg_roughness
    calc_veg_roughness.restype = c_double
    calc_veg_roughness.argtypes = [c_double]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 86
if hasattr(_libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
           'calc_water_balance_error'):
    calc_water_balance_error = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].calc_water_balance_error
    calc_water_balance_error.restype = c_double
    calc_water_balance_error.argtypes = [c_int, c_double, c_double, c_double]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 87
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'CalcAerodynamic'):
    CalcAerodynamic = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].CalcAerodynamic
    CalcAerodynamic.restype = c_int
    CalcAerodynamic.argtypes = [
        c_uint8,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double)]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 89
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'CalcBlowingSnow'):
    CalcBlowingSnow = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].CalcBlowingSnow
    CalcBlowingSnow.restype = c_double
    CalcBlowingSnow.argtypes = [
        c_double,
        c_double,
        c_uint,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_int,
        c_int,
        c_double,
        c_double,
        c_double,
        POINTER(c_double)]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 92
if hasattr(_libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
           'CalcIcePackEnergyBalance'):
    _func = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].CalcIcePackEnergyBalance
    _restype = c_double
    _argtypes = [c_double]
    CalcIcePackEnergyBalance = _variadic_function(_func, _restype, _argtypes)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 93
if hasattr(_libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
           'CalcSnowPackEnergyBalance'):
    _func = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].CalcSnowPackEnergyBalance
    _restype = c_double
    _argtypes = [c_double]
    CalcSnowPackEnergyBalance = _variadic_function(_func, _restype, _argtypes)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 94
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'CalcSubFlux'):
    CalcSubFlux = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].CalcSubFlux
    CalcSubFlux.restype = c_double
    CalcSubFlux.argtypes = [
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        POINTER(c_double)]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 98
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'canopy_assimilation'):
    canopy_assimilation = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].canopy_assimilation
    canopy_assimilation.restype = None
    canopy_assimilation.argtypes = [
        c_char,
        c_double,
        c_double,
        c_double,
        POINTER(c_double),
        c_double,
        c_double,
        POINTER(c_double),
        c_double,
        c_double,
        POINTER(c_double),
        c_double,
        String,
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double)]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 102
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'canopy_evap'):
    canopy_evap = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].canopy_evap
    canopy_evap.restype = c_double
    canopy_evap.argtypes = [
        POINTER(layer_data_struct),
        POINTER(veg_var_struct),
        c_uint8,
        c_uint,
        POINTER(c_double),
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        c_double,
        c_double,
        POINTER(c_double)]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 106
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'colavg'):
    colavg = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].colavg
    colavg.restype = None
    colavg.argtypes = [POINTER(c_double), POINTER(c_double), POINTER(
        c_double), c_double, POINTER(c_double), c_int, c_double, c_double]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 108
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'collect_eb_terms'):
    collect_eb_terms = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].collect_eb_terms
    collect_eb_terms.restype = None
    collect_eb_terms.argtypes = [
        energy_bal_struct,
        snow_data_struct,
        cell_data_struct,
        POINTER(c_int),
        POINTER(c_int),
        POINTER(c_int),
        POINTER(c_int),
        POINTER(c_int),
        c_double,
        c_double,
        c_double,
        c_int,
        c_int,
        c_double,
        c_int,
        c_int,
        POINTER(c_double),
        c_double,
        POINTER(out_data_struct)]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 112
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'collect_wb_terms'):
    collect_wb_terms = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].collect_wb_terms
    collect_wb_terms.restype = None
    collect_wb_terms.argtypes = [
        cell_data_struct,
        veg_var_struct,
        snow_data_struct,
        c_double,
        c_double,
        c_double,
        c_int,
        c_double,
        c_int,
        POINTER(c_double),
        POINTER(c_double),
        POINTER(out_data_struct)]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 115
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'compute_coszen'):
    compute_coszen = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].compute_coszen
    compute_coszen.restype = c_double
    compute_coszen.argtypes = [c_double, c_double, c_double, c_uint, c_uint]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 116
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'compute_pot_evap'):
    compute_pot_evap = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].compute_pot_evap
    compute_pot_evap.restype = None
    compute_pot_evap.argtypes = [
        c_uint,
        POINTER(dmy_struct),
        c_size_t,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        POINTER(
            POINTER(c_double)),
        POINTER(c_double)]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 118
if hasattr(_libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
           'compute_runoff_and_asat'):
    compute_runoff_and_asat = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].compute_runoff_and_asat
    compute_runoff_and_asat.restype = None
    compute_runoff_and_asat.argtypes = [POINTER(soil_con_struct), POINTER(
        c_double), c_double, POINTER(c_double), POINTER(c_double)]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 120
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'compute_soil_resp'):
    compute_soil_resp = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].compute_soil_resp
    compute_soil_resp.restype = None
    compute_soil_resp.argtypes = [
        c_int,
        POINTER(c_double),
        c_double,
        c_double,
        POINTER(c_double),
        POINTER(c_double),
        c_double,
        c_double,
        c_double,
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double)]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 122
if hasattr(_libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
           'compute_soil_layer_thermal_properties'):
    compute_soil_layer_thermal_properties = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].compute_soil_layer_thermal_properties
    compute_soil_layer_thermal_properties.restype = None
    compute_soil_layer_thermal_properties.argtypes = [
        POINTER(layer_data_struct),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        c_size_t]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 126
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'compute_zwt'):
    compute_zwt = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].compute_zwt
    compute_zwt.restype = c_double
    compute_zwt.argtypes = [POINTER(soil_con_struct), c_int, c_double]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 127
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'correct_precip'):
    correct_precip = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].correct_precip
    correct_precip.restype = None
    correct_precip.argtypes = [
        POINTER(c_double), c_double, c_double, c_double, c_double]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 128
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'darkinhib'):
    darkinhib = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].darkinhib
    darkinhib.restype = c_double
    darkinhib.argtypes = [c_double]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 129
if hasattr(_libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
           'distribute_node_moisture_properties'):
    distribute_node_moisture_properties = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].distribute_node_moisture_properties
    distribute_node_moisture_properties.restype = c_int
    distribute_node_moisture_properties.argtypes = [
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        c_int,
        c_int,
        c_char]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 134
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'eddy'):
    eddy = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].eddy
    eddy.restype = None
    eddy.argtypes = [c_int, c_double, POINTER(c_double), POINTER(
        c_double), c_double, c_int, c_double, c_double]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 135
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'energycalc'):
    energycalc = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].energycalc
    energycalc.restype = None
    energycalc.argtypes = [
        POINTER(c_double),
        POINTER(c_double),
        c_int,
        c_double,
        c_double,
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double)]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 137
if hasattr(_libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
           'error_calc_atmos_energy_bal'):
    _func = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].error_calc_atmos_energy_bal
    _restype = c_double
    _argtypes = [c_double]
    error_calc_atmos_energy_bal = _variadic_function(
        _func, _restype, _argtypes)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 138
if hasattr(_libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
           'error_calc_atmos_moist_bal'):
    _func = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].error_calc_atmos_moist_bal
    _restype = c_double
    _argtypes = [c_double]
    error_calc_atmos_moist_bal = _variadic_function(_func, _restype, _argtypes)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 139
if hasattr(_libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
           'error_calc_canopy_energy_bal'):
    _func = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].error_calc_canopy_energy_bal
    _restype = c_double
    _argtypes = [c_double]
    error_calc_canopy_energy_bal = _variadic_function(
        _func, _restype, _argtypes)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 140
if hasattr(_libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
           'error_calc_surf_energy_bal'):
    _func = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].error_calc_surf_energy_bal
    _restype = c_double
    _argtypes = [c_double]
    error_calc_surf_energy_bal = _variadic_function(_func, _restype, _argtypes)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 141
if hasattr(_libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
           'error_print_atmos_energy_bal'):
    error_print_atmos_energy_bal = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].error_print_atmos_energy_bal
    error_print_atmos_energy_bal.restype = c_double
    error_print_atmos_energy_bal.argtypes = [c_double, c_void_p]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 142
if hasattr(_libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
           'error_print_atmos_moist_bal'):
    error_print_atmos_moist_bal = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].error_print_atmos_moist_bal
    error_print_atmos_moist_bal.restype = c_double
    error_print_atmos_moist_bal.argtypes = [c_double, c_void_p]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 143
if hasattr(_libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
           'error_print_canopy_energy_bal'):
    error_print_canopy_energy_bal = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].error_print_canopy_energy_bal
    error_print_canopy_energy_bal.restype = c_double
    error_print_canopy_energy_bal.argtypes = [c_double, c_void_p]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 144
if hasattr(_libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
           'error_print_solve_T_profile'):
    error_print_solve_T_profile = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].error_print_solve_T_profile
    error_print_solve_T_profile.restype = c_double
    error_print_solve_T_profile.argtypes = [c_double, c_void_p]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 145
if hasattr(_libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
           'error_print_surf_energy_bal'):
    error_print_surf_energy_bal = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].error_print_surf_energy_bal
    error_print_surf_energy_bal.restype = c_double
    error_print_surf_energy_bal.argtypes = [c_double, c_void_p]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 146
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'error_solve_T_profile'):
    _func = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].error_solve_T_profile
    _restype = c_double
    _argtypes = [c_double]
    error_solve_T_profile = _variadic_function(_func, _restype, _argtypes)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 147
if hasattr(_libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
           'ErrorIcePackEnergyBalance'):
    _func = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].ErrorIcePackEnergyBalance
    _restype = c_double
    _argtypes = [c_double]
    ErrorIcePackEnergyBalance = _variadic_function(_func, _restype, _argtypes)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 148
if hasattr(_libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
           'ErrorPrintIcePackEnergyBalance'):
    ErrorPrintIcePackEnergyBalance = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].ErrorPrintIcePackEnergyBalance
    ErrorPrintIcePackEnergyBalance.restype = c_double
    ErrorPrintIcePackEnergyBalance.argtypes = [c_double, c_void_p]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 149
if hasattr(_libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
           'ErrorPrintSnowPackEnergyBalance'):
    ErrorPrintSnowPackEnergyBalance = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].ErrorPrintSnowPackEnergyBalance
    ErrorPrintSnowPackEnergyBalance.restype = c_int
    ErrorPrintSnowPackEnergyBalance.argtypes = [c_double, c_void_p]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 150
if hasattr(_libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
           'ErrorSnowPackEnergyBalance'):
    _func = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].ErrorSnowPackEnergyBalance
    _restype = c_int
    _argtypes = [c_double]
    ErrorSnowPackEnergyBalance = _variadic_function(_func, _restype, _argtypes)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 151
if hasattr(_libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
           'estimate_layer_ice_content'):
    estimate_layer_ice_content = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].estimate_layer_ice_content
    estimate_layer_ice_content.restype = c_int
    estimate_layer_ice_content.argtypes = [
        POINTER(layer_data_struct),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        c_double,
        c_size_t,
        c_size_t,
        c_char]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 154
if hasattr(_libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
           'estimate_layer_ice_content_quick_flux'):
    estimate_layer_ice_content_quick_flux = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].estimate_layer_ice_content_quick_flux
    estimate_layer_ice_content_quick_flux.restype = c_int
    estimate_layer_ice_content_quick_flux.argtypes = [
        POINTER(layer_data_struct),
        POINTER(c_double),
        c_double,
        c_double,
        c_double,
        c_double,
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        c_double,
        c_char]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 158
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'estimate_T1'):
    estimate_T1 = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].estimate_T1
    estimate_T1.restype = c_double
    estimate_T1.argtypes = [
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 160
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'exp_interp'):
    exp_interp = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].exp_interp
    exp_interp.restype = c_double
    exp_interp.argtypes = [c_double, c_double, c_double, c_double, c_double]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 161
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'faparl'):
    faparl = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].faparl
    faparl.restype = None
    faparl.argtypes = [
        POINTER(c_double),
        c_double,
        c_double,
        c_double,
        c_double,
        POINTER(c_double),
        POINTER(c_double)]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 162
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'fda_heat_eqn'):
    _func = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].fda_heat_eqn
    _restype = None
    _argtypes = [POINTER(c_double), POINTER(c_double), c_int, c_int]
    fda_heat_eqn = _variadic_function(_func, _restype, _argtypes)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 163
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'fdjac3'):
    fdjac3 = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].fdjac3
    fdjac3.restype = None
    fdjac3.argtypes = [
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        CFUNCTYPE(
            UNCHECKED(None),
            POINTER(c_double),
            POINTER(c_double),
            c_int,
            c_int),
        c_int]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 165
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'find_0_degree_fronts'):
    find_0_degree_fronts = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].find_0_degree_fronts
    find_0_degree_fronts.restype = None
    find_0_degree_fronts.argtypes = [
        POINTER(energy_bal_struct),
        POINTER(c_double),
        POINTER(c_double),
        c_int]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 166
for _lib in _libs.values():
    if hasattr(_lib, 'find_average_layer'):
        find_average_layer = _lib.find_average_layer
        find_average_layer.restype = layer_data_struct
        find_average_layer.argtypes = [
            POINTER(layer_data_struct),
            POINTER(layer_data_struct),
            c_double,
            c_double]
        break

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 168
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'func_atmos_energy_bal'):
    func_atmos_energy_bal = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].func_atmos_energy_bal
    func_atmos_energy_bal.restype = c_double
    func_atmos_energy_bal.argtypes = [c_double, c_void_p]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 169
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'func_atmos_moist_bal'):
    func_atmos_moist_bal = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].func_atmos_moist_bal
    func_atmos_moist_bal.restype = c_double
    func_atmos_moist_bal.argtypes = [c_double, c_void_p]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 170
if hasattr(_libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
           'func_canopy_energy_bal'):
    func_canopy_energy_bal = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].func_canopy_energy_bal
    func_canopy_energy_bal.restype = c_double
    func_canopy_energy_bal.argtypes = [c_double, c_void_p]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 171
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'func_surf_energy_bal'):
    func_surf_energy_bal = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].func_surf_energy_bal
    func_surf_energy_bal.restype = c_double
    func_surf_energy_bal.argtypes = [c_double, c_void_p]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 172
try:
    funcd = (
        POINTER(
            CFUNCTYPE(
                UNCHECKED(c_double),
                c_double,
                c_double,
                c_double,
                c_double,
                c_double,
                c_double,
                c_double,
                c_double,
                c_double,
                c_double,
                c_double))).in_dll(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'funcd')
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 176
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'get_depth'):
    get_depth = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].get_depth
    get_depth.restype = c_int
    get_depth.argtypes = [lake_con_struct, c_double, POINTER(c_double)]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 177
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'get_prob'):
    get_prob = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].get_prob
    get_prob.restype = c_double
    get_prob.argtypes = [c_double, c_double, c_double, c_double]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 178
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'get_sarea'):
    get_sarea = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].get_sarea
    get_sarea.restype = c_int
    get_sarea.argtypes = [lake_con_struct, c_double, POINTER(c_double)]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 179
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'get_shear'):
    get_shear = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].get_shear
    get_shear.restype = None
    get_shear.argtypes = [
        c_double, POINTER(c_double), POINTER(c_double), c_double, c_double]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 180
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'get_thresh'):
    get_thresh = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].get_thresh
    get_thresh.restype = c_double
    get_thresh.argtypes = [c_double, c_double, c_double]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 181
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'get_volume'):
    get_volume = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].get_volume
    get_volume.restype = c_int
    get_volume.argtypes = [lake_con_struct, c_double, POINTER(c_double)]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 182
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'hiTinhib'):
    hiTinhib = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].hiTinhib
    hiTinhib.restype = c_double
    hiTinhib.argtypes = [c_double]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 183
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'initialize_lake'):
    initialize_lake = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].initialize_lake
    initialize_lake.restype = c_int
    initialize_lake.argtypes = [
        POINTER(lake_var_struct),
        lake_con_struct,
        POINTER(soil_con_struct),
        POINTER(cell_data_struct),
        c_double,
        c_int]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 185
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'ice_melt'):
    ice_melt = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].ice_melt
    ice_melt.restype = c_int
    ice_melt.argtypes = [
        c_double,
        c_double,
        POINTER(c_double),
        c_double,
        POINTER(snow_data_struct),
        POINTER(lake_var_struct),
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double)]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 190
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'IceEnergyBalance'):
    IceEnergyBalance = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].IceEnergyBalance
    IceEnergyBalance.restype = c_double
    IceEnergyBalance.argtypes = [c_double, c_void_p]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 191
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'iceform'):
    iceform = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].iceform
    iceform.restype = None
    iceform.argtypes = [
        POINTER(c_double),
        POINTER(c_double),
        c_double,
        c_double,
        POINTER(c_double),
        c_int,
        c_double,
        c_double,
        c_double,
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        c_double]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 193
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'icerad'):
    icerad = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].icerad
    icerad.restype = None
    icerad.argtypes = [c_double, c_double, c_double, POINTER(
        c_double), POINTER(c_double), POINTER(c_double)]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 194
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'lakeice'):
    lakeice = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].lakeice
    lakeice.restype = c_int
    lakeice.argtypes = [
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        POINTER(c_double),
        c_double,
        POINTER(c_double),
        POINTER(c_double),
        c_double,
        c_double]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 196
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'latent_heat_from_snow'):
    latent_heat_from_snow = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].latent_heat_from_snow
    latent_heat_from_snow.restype = None
    latent_heat_from_snow.argtypes = [
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double)]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 199
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'latsens'):
    latsens = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].latsens
    latsens.restype = None
    latsens.argtypes = [
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        POINTER(c_double),
        POINTER(c_double),
        c_double]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 201
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'linear_interp'):
    linear_interp = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].linear_interp
    linear_interp.restype = c_double
    linear_interp.argtypes = [c_double, c_double, c_double, c_double, c_double]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 202
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'lkdrag'):
    lkdrag = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].lkdrag
    lkdrag.restype = c_double
    lkdrag.argtypes = [c_double, c_double, c_double, c_double, c_double]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 203
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'MassRelease'):
    MassRelease = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].MassRelease
    MassRelease.restype = None
    MassRelease.argtypes = [
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double)]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 204
if hasattr(_libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
           'maximum_unfrozen_water'):
    maximum_unfrozen_water = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].maximum_unfrozen_water
    maximum_unfrozen_water.restype = c_double
    maximum_unfrozen_water.argtypes = [c_double, c_double, c_double, c_double]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 205
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'modify_Ksat'):
    modify_Ksat = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].modify_Ksat
    modify_Ksat.restype = c_double
    modify_Ksat.argtypes = [c_double]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 206
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'new_snow_density'):
    new_snow_density = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].new_snow_density
    new_snow_density.restype = c_double
    new_snow_density.argtypes = [c_double]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 207
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'newt_raph'):
    newt_raph = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].newt_raph
    newt_raph.restype = c_int
    newt_raph.argtypes = [CFUNCTYPE(UNCHECKED(None), POINTER(
        c_double), POINTER(c_double), c_int, c_int), POINTER(c_double), c_int]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 209
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'penman'):
    penman = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].penman
    penman.restype = c_double
    penman.argtypes = [
        c_double, c_double, c_double, c_double, c_double, c_double, c_double]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 210
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'photosynth'):
    photosynth = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].photosynth
    photosynth.restype = None
    photosynth.argtypes = [
        c_char,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        String,
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double)]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 213
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'polint'):
    polint = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].polint
    polint.restype = None
    polint.argtypes = [POINTER(c_double), POINTER(
        c_double), c_int, c_double, POINTER(c_double), POINTER(c_double)]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 214
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'prepare_full_energy'):
    prepare_full_energy = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].prepare_full_energy
    prepare_full_energy.restype = None
    prepare_full_energy.argtypes = [c_int, POINTER(all_vars_struct), POINTER(
        soil_con_struct), POINTER(c_double), POINTER(c_double)]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 216
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'put_data'):
    put_data = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].put_data
    put_data.restype = c_int
    put_data.argtypes = [
        POINTER(all_vars_struct),
        POINTER(atmos_data_struct),
        POINTER(soil_con_struct),
        POINTER(veg_con_struct),
        POINTER(veg_lib_struct),
        POINTER(lake_con_struct),
        POINTER(out_data_struct),
        POINTER(save_data_struct),
        c_int]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 219
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'qromb'):
    qromb = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].qromb
    qromb.restype = c_double
    qromb.argtypes = [
        CFUNCTYPE(
            UNCHECKED(c_double),
        ),
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 222
if hasattr(_libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
           'rescale_snow_energy_fluxes'):
    rescale_snow_energy_fluxes = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].rescale_snow_energy_fluxes
    rescale_snow_energy_fluxes.restype = None
    rescale_snow_energy_fluxes.argtypes = [
        c_double,
        c_double,
        POINTER(snow_data_struct),
        POINTER(energy_bal_struct)]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 224
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'rescale_snow_storage'):
    rescale_snow_storage = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].rescale_snow_storage
    rescale_snow_storage.restype = None
    rescale_snow_storage.argtypes = [
        c_double, c_double, POINTER(snow_data_struct)]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 225
if hasattr(_libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
           'rescale_soil_veg_fluxes'):
    rescale_soil_veg_fluxes = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].rescale_soil_veg_fluxes
    rescale_soil_veg_fluxes.restype = None
    rescale_soil_veg_fluxes.argtypes = [
        c_double, c_double, POINTER(cell_data_struct), POINTER(veg_var_struct)]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 227
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'rhoinit'):
    rhoinit = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].rhoinit
    rhoinit.restype = None
    rhoinit.argtypes = [POINTER(c_double), c_double]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 228
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'root_brent'):
    _func = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].root_brent
    _restype = c_double
    _argtypes = [c_double, c_double, String, CFUNCTYPE(
        UNCHECKED(c_double), c_double, c_void_p)]
    root_brent = _variadic_function(_func, _restype, _argtypes)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 230
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'rtnewt'):
    rtnewt = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].rtnewt
    rtnewt.restype = c_double
    rtnewt.argtypes = [c_double, c_double, c_double, c_double, c_double]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 231
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'runoff'):
    runoff = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].runoff
    runoff.restype = c_int
    runoff.argtypes = [
        POINTER(cell_data_struct),
        POINTER(energy_bal_struct),
        POINTER(soil_con_struct),
        c_double,
        POINTER(c_double),
        c_double,
        c_int]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 233
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'set_node_parameters'):
    set_node_parameters = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].set_node_parameters
    set_node_parameters.restype = None
    set_node_parameters.argtypes = [
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        c_int,
        c_int]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 236
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'shear_stress'):
    shear_stress = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].shear_stress
    shear_stress.restype = None
    shear_stress.argtypes = [
        c_double, c_double, POINTER(c_double), POINTER(c_double), c_double]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 238
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'snow_albedo'):
    snow_albedo = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].snow_albedo
    snow_albedo.restype = c_double
    snow_albedo.argtypes = [
        c_double, c_double, c_double, c_double, c_double, c_int, c_char]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 239
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'snow_density'):
    snow_density = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].snow_density
    snow_density.restype = c_double
    snow_density.argtypes = [
        POINTER(snow_data_struct), c_double, c_double, c_double, c_double]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 240
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'snow_intercept'):
    snow_intercept = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].snow_intercept
    snow_intercept.restype = c_int
    snow_intercept.argtypes = [
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_uint8),
        POINTER(c_uint),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        c_int,
        c_int,
        c_int,
        c_int,
        c_int,
        c_int,
        c_uint,
        POINTER(c_double),
        POINTER(c_double),
        POINTER(atmos_data_struct),
        POINTER(layer_data_struct),
        POINTER(soil_con_struct),
        POINTER(veg_var_struct)]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 249
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'snow_melt'):
    snow_melt = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].snow_melt
    snow_melt.restype = c_int
    snow_melt.argtypes = [
        c_double,
        c_double,
        c_double,
        c_double,
        POINTER(c_double),
        c_double,
        POINTER(c_double),
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        c_int,
        c_int,
        c_int,
        c_int,
        POINTER(snow_data_struct)]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 255
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'SnowPackEnergyBalance'):
    SnowPackEnergyBalance = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].SnowPackEnergyBalance
    SnowPackEnergyBalance.restype = c_double
    SnowPackEnergyBalance.argtypes = [c_double, c_void_p]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 256
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'soil_carbon_balance'):
    soil_carbon_balance = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].soil_carbon_balance
    soil_carbon_balance.restype = None
    soil_carbon_balance.argtypes = [POINTER(soil_con_struct), POINTER(
        energy_bal_struct), POINTER(cell_data_struct), POINTER(veg_var_struct)]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 258
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'soil_conductivity'):
    soil_conductivity = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].soil_conductivity
    soil_conductivity.restype = c_double
    soil_conductivity.argtypes = [
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 260
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'soil_thermal_eqn'):
    soil_thermal_eqn = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].soil_thermal_eqn
    soil_thermal_eqn.restype = c_double
    soil_thermal_eqn.argtypes = [c_double, c_void_p]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 261
if hasattr(_libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
           'solve_atmos_energy_bal'):
    _func = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].solve_atmos_energy_bal
    _restype = c_double
    _argtypes = [c_double]
    solve_atmos_energy_bal = _variadic_function(_func, _restype, _argtypes)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 262
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'solve_atmos_moist_bal'):
    _func = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].solve_atmos_moist_bal
    _restype = c_double
    _argtypes = [c_double]
    solve_atmos_moist_bal = _variadic_function(_func, _restype, _argtypes)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 263
if hasattr(_libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
           'solve_canopy_energy_bal'):
    _func = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].solve_canopy_energy_bal
    _restype = c_double
    _argtypes = [c_double]
    solve_canopy_energy_bal = _variadic_function(_func, _restype, _argtypes)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 264
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'solve_lake'):
    solve_lake = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].solve_lake
    solve_lake.restype = c_int
    solve_lake.argtypes = [
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        POINTER(lake_var_struct),
        soil_con_struct,
        c_double,
        c_double,
        dmy_struct,
        c_double]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 267
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'solve_snow'):
    solve_snow = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].solve_snow
    solve_snow.restype = c_double
    solve_snow.argtypes = [
        c_char,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        c_int,
        c_size_t,
        c_uint,
        c_uint,
        c_double,
        c_size_t,
        c_size_t,
        c_int,
        POINTER(c_int),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(dmy_struct),
        POINTER(atmos_data_struct),
        POINTER(energy_bal_struct),
        POINTER(layer_data_struct),
        POINTER(snow_data_struct),
        POINTER(soil_con_struct),
        POINTER(veg_var_struct)]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 277
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'solve_surf_energy_bal'):
    _func = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].solve_surf_energy_bal
    _restype = c_double
    _argtypes = [c_double]
    solve_surf_energy_bal = _variadic_function(_func, _restype, _argtypes)

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 278
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'solve_T_profile'):
    solve_T_profile = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].solve_T_profile
    solve_T_profile.restype = c_int
    solve_T_profile.argtypes = [
        POINTER(c_double),
        POINTER(c_double),
        String,
        POINTER(c_uint),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        c_double,
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        c_double,
        c_int,
        POINTER(c_int),
        c_int,
        c_int,
        c_int]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 282
if hasattr(_libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
           'solve_T_profile_implicit'):
    solve_T_profile_implicit = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].solve_T_profile_implicit
    solve_T_profile_implicit.restype = c_int
    solve_T_profile_implicit.argtypes = [
        POINTER(c_double),
        POINTER(c_double),
        String,
        POINTER(c_uint),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        c_double,
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        c_double,
        c_int,
        POINTER(c_int),
        c_int,
        c_int,
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double)]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 288
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'specheat'):
    specheat = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].specheat
    specheat.restype = c_double
    specheat.argtypes = [c_double]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 289
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'StabilityCorrection'):
    StabilityCorrection = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].StabilityCorrection
    StabilityCorrection.restype = c_double
    StabilityCorrection.argtypes = [
        c_double, c_double, c_double, c_double, c_double, c_double]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 290
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'sub_with_height'):
    sub_with_height = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].sub_with_height
    sub_with_height.restype = c_double
    sub_with_height.argtypes = [
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 293
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'surface_fluxes'):
    surface_fluxes = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].surface_fluxes
    surface_fluxes.restype = c_int
    surface_fluxes.argtypes = [
        c_uint8,
        c_double,
        c_double,
        c_double,
        c_double,
        POINTER(c_double),
        POINTER(c_double),
        POINTER(
            POINTER(c_double)),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        c_size_t,
        c_size_t,
        c_uint,
        c_double,
        c_uint,
        c_size_t,
        c_uint,
        POINTER(atmos_data_struct),
        POINTER(dmy_struct),
        POINTER(energy_bal_struct),
        POINTER(global_param_struct),
        POINTER(cell_data_struct),
        POINTER(snow_data_struct),
        POINTER(soil_con_struct),
        POINTER(veg_var_struct),
        c_double,
        c_double,
        c_double,
        POINTER(c_double)]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 301
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'svp'):
    svp = _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].svp
    svp.restype = c_double
    svp.argtypes = [c_double]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 302
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'svp_slope'):
    svp_slope = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].svp_slope
    svp_slope.restype = c_double
    svp_slope.argtypes = [c_double]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 303
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'temp_area'):
    temp_area = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].temp_area
    temp_area.restype = None
    temp_area.argtypes = [
        c_double,
        c_double,
        c_double,
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        c_double,
        POINTER(c_double),
        c_int,
        c_double,
        c_double,
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double)]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 306
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'tracer_mixer'):
    tracer_mixer = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].tracer_mixer
    tracer_mixer.restype = None
    tracer_mixer.argtypes = [POINTER(c_double), POINTER(c_int), POINTER(
        c_double), c_int, c_double, c_double, POINTER(c_double)]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 307
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'transpiration'):
    transpiration = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].transpiration
    transpiration.restype = None
    transpiration.argtypes = [
        POINTER(layer_data_struct),
        POINTER(veg_var_struct),
        c_uint,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        c_double,
        c_double,
        POINTER(c_double)]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 311
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'transport_with_height'):
    transport_with_height = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].transport_with_height
    transport_with_height.restype = c_double
    transport_with_height.argtypes = [
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 314
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'trapzd'):
    trapzd = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].trapzd
    trapzd.restype = c_double
    trapzd.argtypes = [
        CFUNCTYPE(
            UNCHECKED(c_double),
        ),
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_int]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 317
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'tridia'):
    tridia = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].tridia
    tridia.restype = None
    tridia.argtypes = [c_int, POINTER(c_double), POINTER(
        c_double), POINTER(c_double), POINTER(c_double), POINTER(c_double)]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 318
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'tridiag'):
    tridiag = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].tridiag
    tridiag.restype = None
    tridiag.argtypes = [POINTER(c_double), POINTER(
        c_double), POINTER(c_double), POINTER(c_double), c_uint]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 319
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'vic_run'):
    vic_run = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].vic_run
    vic_run.restype = c_int
    vic_run.argtypes = [
        c_int,
        POINTER(atmos_data_struct),
        POINTER(all_vars_struct),
        POINTER(dmy_struct),
        POINTER(global_param_struct),
        POINTER(lake_con_struct),
        POINTER(soil_con_struct),
        POINTER(veg_con_struct),
        POINTER(veg_lib_struct),
        POINTER(veg_hist_struct)]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 322
if hasattr(_libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
           'volumetric_heat_capacity'):
    volumetric_heat_capacity = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].volumetric_heat_capacity
    volumetric_heat_capacity.restype = c_double
    volumetric_heat_capacity.argtypes = [
        c_double, c_double, c_double, c_double]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 323
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'water_balance'):
    water_balance = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].water_balance
    water_balance.restype = c_int
    water_balance.argtypes = [
        POINTER(lake_var_struct),
        lake_con_struct,
        c_double,
        POINTER(all_vars_struct),
        c_int,
        c_int,
        c_int,
        c_double,
        soil_con_struct,
        veg_con_struct]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 325
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'water_energy_balance'):
    water_energy_balance = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].water_energy_balance
    water_energy_balance.restype = c_int
    water_energy_balance.argtypes = [
        c_int,
        POINTER(c_double),
        POINTER(c_double),
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        c_double,
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        c_double]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 331
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'water_under_ice'):
    water_under_ice = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].water_under_ice
    water_under_ice.restype = c_int
    water_under_ice.argtypes = [
        c_int,
        c_double,
        c_double,
        POINTER(c_double),
        POINTER(c_double),
        c_double,
        c_int,
        c_double,
        c_double,
        c_double,
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        c_int,
        c_double,
        c_double,
        c_double,
        POINTER(c_double)]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 334
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'wrap_compute_zwt'):
    wrap_compute_zwt = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].wrap_compute_zwt
    wrap_compute_zwt.restype = None
    wrap_compute_zwt.argtypes = [
        POINTER(soil_con_struct), POINTER(cell_data_struct)]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 335
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'write_layer'):
    write_layer = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].write_layer
    write_layer.restype = None
    write_layer.argtypes = [
        POINTER(layer_data_struct), c_int, POINTER(c_double)]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 336
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'write_vegvar'):
    write_vegvar = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].write_vegvar
    write_vegvar.restype = None
    write_vegvar.argtypes = [POINTER(veg_var_struct), c_int]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_run.h: 337
if hasattr(
        _libs['./build/lib.macosx-10.5-x86_64-2.7/vic_run.so'],
        'zero_output_list'):
    zero_output_list = _libs[
        './build/lib.macosx-10.5-x86_64-2.7/vic_run.so'].zero_output_list
    zero_output_list.restype = None
    zero_output_list.argtypes = [POINTER(out_data_struct)]

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_physical_constants.h: 30
try:
    DAYS_PER_360DAY_YEAR = 360
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_physical_constants.h: 30
try:
    DAYS_PER_YEAR = 365
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_physical_constants.h: 30
try:
    DAYS_PER_LYEAR = 366
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_physical_constants.h: 30
try:
    DAYS_PER_JYEAR = 365.25
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_physical_constants.h: 30
try:
    HOURS_PER_DAY = 24
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_physical_constants.h: 30
try:
    MONTHS_PER_YEAR = 12
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_physical_constants.h: 30
try:
    MIN_PER_HOUR = 60
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_physical_constants.h: 30
try:
    MIN_PER_DAY = (MIN_PER_HOUR * HOURS_PER_DAY)
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_physical_constants.h: 30
try:
    SEC_PER_MIN = 60
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_physical_constants.h: 30
try:
    SEC_PER_HOUR = (SEC_PER_MIN * MIN_PER_HOUR)
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_physical_constants.h: 30
try:
    SEC_PER_DAY = (SEC_PER_HOUR * HOURS_PER_DAY)
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_physical_constants.h: 32
try:
    JOULES_PER_CAL = 4.1868
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_physical_constants.h: 32
try:
    GRAMS_PER_KG = 1000
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_physical_constants.h: 32
try:
    PA_PER_KPA = 1000
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_physical_constants.h: 32
try:
    BAR_PER_KPA = 100
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_physical_constants.h: 32
try:
    RAD_PER_DEG = 0.0174532925
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_physical_constants.h: 32
try:
    M_PER_KM = 1000
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_physical_constants.h: 32
try:
    MM_PER_M = 1000
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_physical_constants.h: 32
try:
    CM_PER_M = 100
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_physical_constants.h: 32
try:
    MM_PER_CM = 10
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_physical_constants.h: 32
try:
    MOLE_PER_KMOLE = 1000
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_physical_constants.h: 32
try:
    FRACT_TO_PERCENT = 100
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_physical_constants.h: 32
try:
    PPM_to_MIXRATIO = 1e-06
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_physical_constants.h: 34
try:
    CONST_PI = 3.141592653589793
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_physical_constants.h: 36
try:
    CONST_CDAY = 86400
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_physical_constants.h: 36
try:
    CONST_SDAY = 86164
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_physical_constants.h: 36
try:
    CONST_DDAYS_PER_YEAR = 365.2425
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_physical_constants.h: 38
try:
    CONST_OMEGA = ((2.0 * CONST_PI) / CONST_SDAY)
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_physical_constants.h: 38
try:
    CONST_SECPERRAD = 13750.9871
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_physical_constants.h: 38
try:
    CONST_RADPERDAY = 0.017214
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_physical_constants.h: 38
try:
    CONST_RADPERDEG = 0.01745329
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_physical_constants.h: 38
try:
    CONST_MINDECL = (-0.4092797)
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_physical_constants.h: 38
try:
    CONST_DAYSOFF = 11.25
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_physical_constants.h: 40
try:
    CONST_REARTH = 6371220.0
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_physical_constants.h: 40
try:
    CONST_G = 9.80616
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_physical_constants.h: 40
try:
    CONST_STEBOL = 5.67e-08
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_physical_constants.h: 40
try:
    CONST_BOLTZ = 1.38065e-23
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_physical_constants.h: 40
try:
    CONST_AVOGAD = 6.02214e+26
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_physical_constants.h: 40
try:
    CONST_KARMAN = 0.4
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_physical_constants.h: 41
try:
    CONST_MWDAIR = 28.966
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_physical_constants.h: 41
try:
    CONST_MWWV = 18.016
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_physical_constants.h: 41
try:
    CONST_MWCO2 = 44.011
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_physical_constants.h: 41
try:
    CONST_MWAIR = 28.97
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_physical_constants.h: 41
try:
    CONST_MWC = 12.01
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_physical_constants.h: 42
try:
    CONST_RGAS = (CONST_AVOGAD * CONST_BOLTZ)
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_physical_constants.h: 42
try:
    CONST_RDAIR = (CONST_RGAS / CONST_MWDAIR)
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_physical_constants.h: 42
try:
    CONST_RWV = (CONST_RGAS / CONST_MWWV)
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_physical_constants.h: 42
try:
    CONST_EPS = (CONST_MWWV / CONST_MWAIR)
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_physical_constants.h: 43
try:
    CONST_TKTRIP = 273.16
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_physical_constants.h: 43
try:
    CONST_TKFRZ = 273.15
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_physical_constants.h: 44
try:
    CONST_PSTD = 101325.0
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_physical_constants.h: 44
try:
    CONST_TSTD = (CONST_TKFRZ + 15.0)
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_physical_constants.h: 45
try:
    CONST_RHODAIR = (CONST_PSTD / (CONST_RDAIR * CONST_TKFRZ))
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_physical_constants.h: 45
try:
    CONST_RHOFW = 1000.0
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_physical_constants.h: 45
try:
    CONST_RHOICE = 917.0
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_physical_constants.h: 46
try:
    CONST_CPDAIR = 1004.64
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_physical_constants.h: 46
try:
    CONST_CPMAIR = 1013.0
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_physical_constants.h: 46
try:
    CONST_CPWV = 1810.0
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_physical_constants.h: 46
try:
    CONST_CPFW = 4188.0
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_physical_constants.h: 46
try:
    CONST_CPFWICE = 4200.0
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_physical_constants.h: 46
try:
    CONST_CPICE = 2117.27
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_physical_constants.h: 47
try:
    CONST_LATICE = 333700.0
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_physical_constants.h: 47
try:
    CONST_LATVAP = 2501000.0
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_physical_constants.h: 47
try:
    CONST_LATSUB = (CONST_LATICE + CONST_LATVAP)
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_physical_constants.h: 49
try:
    CONST_SPVAL = 1e+30
except:
    pass

# /usr/include/sys/errno.h: 81
try:
    errno = ((__error())[0])
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_log.h: 67
try:
    clean_errno = (errno == 0) and 'None' or (strerror(errno))
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 42
try:
    MAXSTRING = 2048
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 42
try:
    MISSING = (-99999.0)
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 42
try:
    NODATA_VH = (-1)
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 42
try:
    ERROR = (-999)
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 44
try:
    MAX_VEG = 12
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 44
try:
    MAX_LAYERS = 3
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 44
try:
    MAX_NODES = 50
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 44
try:
    MAX_BANDS = 10
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 44
try:
    MAX_FRONTS = 3
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 44
try:
    MAX_FROST_AREAS = 10
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 44
try:
    MAX_LAKE_NODES = 20
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 44
try:
    MAX_ZWTVMOIST = 11
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 46
try:
    MINSOILDEPTH = 0.001
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 46
try:
    MIN_VEGCOVER = 0.0001
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 48
try:
    MIN_SUBDAILY_STEPS_PER_DAY = 4
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 48
try:
    MAX_SUBDAILY_STEPS_PER_DAY = 1440
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 50
try:
    N_PET_TYPES = 0
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 50
try:
    N_PET_TYPES_NON_NAT = 0
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 50
try:
    PET_SATSOIL = 0
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 50
try:
    PET_H2OSURF = 1
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 50
try:
    PET_SHORT = 2
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 50
try:
    PET_TALL = 3
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 50
try:
    N_PET_TYPES_NAT = 0
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 50
try:
    PET_NATVEG = 4
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 50
try:
    PET_VEGNOCR = 5
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 68
try:
    WET = 0
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 68
try:
    DRY = 1
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 71
try:
    RAIN = 0
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 71
try:
    SNOW = 1
except:
    pass

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 73


def min(a, b):
    return (a < b) and a or b

# /Users/jhamman/Dropbox/src/VIC/vic_run/include/vic_def.h: 73


def max(a, b):
    return (a > b) and a or b

# No inserted files
