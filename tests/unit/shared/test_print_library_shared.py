import pytest
from vic import lib as vic_lib, ffi


# TODO: fix make_all_vars struct

# @pytest.fixture()
# def all_data(scope='function'):
#     return vic_lib.make_all_vars(4)


@pytest.fixture()
def dmy(scope='function'):
    d = ffi.new('dmy_struct *')
    d[0].day = 1
    d[0].day_in_year = 1
    d[0].dayseconds = 1234
    d[0].month = 1
    d[0].year = 1984
    return d


# def test_print_all_data(all_data):
#     assert vic_lib.print_cell_data(all_data.cell, vic_lib.options.Nlayer,
#                                    vic_lib.options.Nfrost, 0) is None


def test_print_dmy(dmy):
    assert vic_lib.print_dmy(dmy) is None


# def test_print_energy_bal(all_data):
#     assert vic_lib.print_energy_bal(all_data.energy, 1, 1) is None


def test_print_filenames():
    assert vic_lib.print_filenames(ffi.addressof(vic_lib.filenames)) is None


def test_print_filep():
    assert vic_lib.print_filep(ffi.addressof(vic_lib.filep)) is None


def test_print_global_param():
    assert vic_lib.print_global_param(
        ffi.addressof(vic_lib.global_param)) is None


def test_print_option():
    assert vic_lib.print_option(ffi.addressof(vic_lib.options)) is None


def test_print_parameters():
    assert vic_lib.print_parameters(ffi.addressof(vic_lib.param)) is None


# def test_print_snow_data(all_data):
#     assert vic_lib.print_snow_data(all_data.snow) is None
