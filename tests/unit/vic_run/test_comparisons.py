from vic import lib as vic_lib
from vic.pycompat import iteritems, pyzip
import numpy as np
import pytest


def test_assert_close_float():
    xs = [1e10, 1e-7, 1e10, 1e-8, 1e10, 1e-8, 1.]
    ys = [1.00001e10, 1e-8, 1.00001e10, 1e-9, 1.0001e10, 1e-9, 1.]

    results = {}
    results[(1e-05, 1e-08)] = [True, False, True, True, False, True, True]
    results[(1e-03, 0)] = [True, False,  True, False,  True, False,  True]

    for tols, rs in iteritems(results):
        rtol, abs_tol = tols
        for x, y, r in pyzip(xs, ys, rs):
            assert vic_lib.assert_close_float(x, y, rtol, abs_tol) == r
            # make sure numpy agrees
            assert np.isclose(x, y, rtol=rtol, atol=abs_tol) == r

def test_assert_close_double():
    xs = [1e10, 1e-7, 1e10, 1e-8, 1e10, 1e-8, 1.]
    ys = [1.00001e10, 1e-8, 1.00001e10, 1e-9, 1.0001e10, 1e-9, 1.]

    results = {}
    results[(1e-05, 1e-08)] = [True, False, True, True, False, True, True]
    results[(1e-03, 0)] = [True, False,  True, False,  True, False,  True]

    for tols, rs in iteritems(results):
        rtol, abs_tol = tols
        for x, y, r in pyzip(xs, ys, rs):
            assert vic_lib.assert_close_float(x, y, rtol, abs_tol) == r
            # make sure numpy agrees
            assert np.isclose(x, y, rtol=rtol, atol=abs_tol) == r
