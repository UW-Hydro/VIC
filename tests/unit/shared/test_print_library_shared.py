import pytest
from vic import lib as vic_lib, ffi


@pytest.fixture()
def dmy(scope='function'):
    d = ffi.new('dmy_struct *')
    d[0].day = 1
    d[0].day_in_year = 1
    d[0].dayseconds = 1234
    d[0].month = 1
    d[0].year = 1984
    return d


def test_print_dmy(dmy):
    assert vic_lib.print_dmy(dmy) is None


def test_print_global_param():
    assert vic_lib.print_global_param(
        ffi.addressof(vic_lib.global_param)) is None


def test_print_option():
    assert vic_lib.print_option(ffi.addressof(vic_lib.options)) is None


def test_print_parameters():
    assert vic_lib.print_parameters(ffi.addressof(vic_lib.param)) is None
