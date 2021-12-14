"""
  @section DESCRIPTION

  cffi wrapper for vic_run
"""

from ._vic import ffi


def _load_lib(lib):
    import os
    import sysconfig
    suffix = sysconfig.get_config_var('SO')
    path = os.path.join(os.path.dirname(__file__), os.pardir,
                        '{0}{1}'.format(lib, suffix))

    return ffi.dlopen(path)


lib = _load_lib('vic_core')

# Initialize global structures
lib.initialize_log()
lib.initialize_global()
lib.initialize_options()
lib.initialize_parameters()

# TODO: wrappers for individual vic functions. For now, access to lib functions
# is made directly through the lib object
