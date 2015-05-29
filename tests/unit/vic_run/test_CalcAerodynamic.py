import ctypes
import numpy as np
from vic.vic import CalcAerodynamic


# segfaulting...
# def test_calc_aerodynamic():
#     r_a = np.array([0, 0], dtype=np.float)
#     u = np.array([0, 0], dtype=np.float)
#     displacement = np.array([0, 0], dtype=np.float)
#     ref_height = np.array([0, 0], dtype=np.float)
#     roughness = np.array([0, 0], dtype=np.float)

#     assert CalcAerodynamic(True, 40., 0.2, 0.0005, 0.001, 0.5,
#                            ctypes.c_void_p(r_a.ctypes.data),
#                            ctypes.c_void_p(u.ctypes.data),
#                            ctypes.c_void_p(displacement.ctypes.data),
#                            ctypes.c_void_p(ref_height.ctypes.data),
#                            ctypes.c_void_p(roughness.ctypes.data)) == 0.
