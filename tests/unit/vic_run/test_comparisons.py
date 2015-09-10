from vic import lib as vic_lib
from vic.pycompat import iteritems, pyzip
import numpy as np
import pytest


def test_assert_almost_equal():
    xs = [1.0, 2.33333, 33, 77.7, 0.000001]
    ys = [1.0, 2.33339, 34, 77.70000000001, -0.00001]

    results = {}
    results[3] = [True, True, False, True, True]
    results[10] = [True, False, False, True, False]
    results[12] = [True, False, False, False, False]

    for decimal, rs in iteritems(results):
        for x, y, r in pyzip(xs, ys, rs):
            assert vic_lib.assert_almost_equal(x, y, decimal) == r
            # make sure numpy agrees
            if not r:
                with pytest.raises(AssertionError):
                    np.testing.assert_almost_equal(x, y, decimal=decimal)
            else:
                np.testing.assert_almost_equal(x, y, decimal=decimal)


def test_assert_close():
    xs = [1e10, 1e-7, 1e10, 1e-8, 1e10, 1e-8, 1.]
    ys = [1.00001e10, 1e-8, 1.00001e10, 1e-9, 1.0001e10, 1e-9, 1.]

    results = {}
    results[(1e-05, 1e-08)] = [True, False, True, True, False, True, True]
    results[(1e-03, 0)] = [True, False,  True, False,  True, False,  True]

    for tols, rs in iteritems(results):
        rtol, abs_tol = tols
        for x, y, r in pyzip(xs, ys, rs):
            assert vic_lib.assert_close(x, y, rtol, abs_tol) == r
            # make sure numpy agrees
            assert np.isclose(x, y, rtol=rtol, atol=abs_tol) == r
