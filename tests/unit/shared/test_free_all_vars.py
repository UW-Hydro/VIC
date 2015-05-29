import pytest
from vic.vic import make_all_vars, free_all_vars
import ctypes


@pytest.fixture()
def all_vars(scope='function'):
    return make_all_vars(12)


def test_free_all_vars(all_vars):
    assert free_all_vars(ctypes.byref(all_vars), 12) is None
