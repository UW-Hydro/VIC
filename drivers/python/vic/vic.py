"""
  @section DESCRIPTION

  ctypes wrapper for vic_run

  This module initialize a few globals missed by ctypes gen and imports all the
  objects from the vic_core.so. This is the module that should be used for
  interacting with the ctypes objects in _vic_run_lib

  @section LICENSE

  The Variable Infiltration Capacity (VIC) macroscale hydrological model
  Copyright (C) 2014 The Land Surface Hydrology Group, Department of Civil
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

import ctypes
from ._vic_run_lib import *
from ._vic_run_lib import _libs

# Map global scalars
flag = ctypes.c_int.in_dll(_libs['vic_core.so'], 'flag')
NR = ctypes.c_size_t.in_dll(_libs['vic_core.so'], 'NR')
NF = ctypes.c_size_t.in_dll(_libs['vic_core.so'], 'NF')

# Map global structures
global_param = global_param_struct.in_dll(_libs['vic_core.so'], 'global_param')
options = option_struct.in_dll(_libs['vic_core.so'], 'options')
param = parameters_struct.in_dll(_libs['vic_core.so'], 'param')
filenames = filenames_struct.in_dll(_libs['vic_core.so'], 'filenames')
filep = filep_struct.in_dll(_libs['vic_core.so'], 'filep')
param_set = param_set_struct.in_dll(_libs['vic_core.so'], 'param_set')

# Initialize global structures
initialize_log()
initialize_global()
initialize_options()
initialize_parameters()
initialize_filenames()
initialize_fileps()
