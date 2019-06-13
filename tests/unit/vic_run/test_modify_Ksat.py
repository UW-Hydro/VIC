from vic import lib as vic_lib
import numpy as np
from scipy.interpolate import interp1d


def test_linear_interp():
    n = 10
    x = np.arange(n, dtype=np.float)
    y = np.arange(n, dtype=np.float)
    new_x = np.arange(n, dtype=np.float) + 0.5
    f = interp1d(x, y, bounds_error=False, fill_value='extrapolate')
    scipy_new_y = f(new_x)
    np.testing.assert_almost_equal(scipy_new_y[:5], [0.5, 1.5, 2.5, 3.5, 4.5])
    for i in range(n - 1):
        assert vic_lib.linear_interp(new_x[i], x[i], x[i + 1],
                                     y[i], y[i + 1]) == scipy_new_y[i]
