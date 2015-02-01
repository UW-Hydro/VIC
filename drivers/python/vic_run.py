"""
  @section DESCRIPTION

  ctypes wrapper for vic_run.c

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

import os
import numpy as np
import ctypes

debug = True

SHAREDOBJECT = 'vic_run.so'
LIBPATH = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../'))

try:
    _clib = np.ctypeslib.load_library(SHAREDOBJECT, LIBPATH)
except (ImportError, OSError) as e:
    print(('Error looking for shared object {0} in {1}'.format(SHAREDOBJECT,
                                                               LIBPATH)))
    raise e

_args = [ctypes.c_int,  # rec
         ctypes.py_object,  # atmos_data_struct
         ctypes.py_object,  # all_vars_struct
         ctypes.py_object,  # dmy_struct
         ctypes.py_object,  # global_param_struct
         ctypes.py_object,  # lake_con_struct
         ctypes.py_object,  # soil_con_struct
         ctypes.py_object,  # veg_con_struct
         ctypes.py_object,  # veg_lib_struct
         ctypes.py_object]  # veg_hist_struct
_clib.vic_run.argtypes = _args
_clib.vic_run.restype = ctypes.c_int


class VICInputError(RuntimeError):
    pass


class VICRuntimeError(RuntimeError):
    pass


def vic_run(*args):
    """
    Wrapper function for vic_run().  This function controls the model core,
    it solves both the energy and water balance models.

    args:
        rec (int)
        atmos_data_struct (pointer to structure)
        all_vars_struct (pointer to structure)
        dmy_struct (pointer to structure)
        global_param_struct (pointer to structure)
        lake_con_struct (pointer to structure)
        soil_con_struct (pointer to structure)
        veg_con_struct (pointer to structure)
        veg_lib_struct (pointer to structure)
        veg_hist_struct (pointer to structure)

    return:
        error (int)
    """
    if debug and len(args) != 10:
        raise VICInputError("Wrong number of arguments provided to VIC run")

    error = _clib.vic_run(*args)

    if error:
        raise VICRuntimeError("vic_run returned an error: %s", error)

    return
