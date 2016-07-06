"""
  @section DESCRIPTION

  cffi wrapper for vic_run

  @section LICENSE

  The Variable Infiltration Capacity (VIC) macroscale hydrological model
  Copyright (C) 2016 The Computational Hydrology Group, Department of Civil
  and Environmental Engineering, University of Washington.

  The VIC model is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2
  of the License, or (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License along with
  this program; if not, write to the Free Software Foundation, Inc.,
  51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
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
