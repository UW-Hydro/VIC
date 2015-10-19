from vic import lib as vic_lib
from vic.pycompat import pyzip
import numpy as np


def test_assert_close_float():
    xs = [1e-7, 1e10, 1e-8, 1e10, 1e-8, 1.]
    ys = [1e-8, 1.00001e10, 1e-9, 1.0001e10, 1e-9, 1.]

    for (rtol, abs_tol) in [(1e-3, 1e-2), (1e-03, 0)]:
        rs = np.isclose(xs, ys, rtol, abs_tol)
        for x, y, r in pyzip(xs, ys, rs):
            assert vic_lib.assert_close_float(x, y, rtol, abs_tol) == r


def test_assert_close_double():
    xs = [1e-7, 1e10, 1e-8, 1e10, 1e-8, 1.]
    ys = [1e-8, 1.00001e10, 1e-9, 1.0001e10, 1e-9, 1.]

    for (rtol, abs_tol) in [(1e-05, 1e-08), (1e-03, 0)]:
        rs = np.isclose(xs, ys, rtol, abs_tol)
        for x, y, r in pyzip(xs, ys, rs):
            assert vic_lib.assert_close_double(x, y, rtol, abs_tol) == r
