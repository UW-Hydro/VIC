import pytest
from vic.vic import ffi
from vic import lib as vic_lib


# TODO: get this to work
# @pytest.fixture()
# def all_vars(scope='function'):
#     return vic_lib.make_all_vars(12)


# def test_free_all_vars():
#     all_vars_p = ffi.new('all_vars_struct *')
#     assert vic_lib.free_all_vars(all_vars_p, 12) is None
