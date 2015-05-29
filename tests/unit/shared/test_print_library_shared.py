import pytest
import ctypes
from vic.vic import (make_all_vars, dmy_struct, options, filenames, filep,
                     global_param, param)
from vic.vic import (print_cell_data, print_dmy, print_energy_bal,
                     print_filenames, print_filep, print_global_param,
                     print_option, print_parameters, print_snow_data)


@pytest.fixture()
def all_data(scope='function'):
    return make_all_vars(4)


@pytest.fixture()
def dmy(scope='function'):
    d = dmy_struct()
    d.day = 1
    d.day_in_year = 1
    d.dayseconds = 1234
    d.month = 1
    d.year = 1984
    return d


def test_print_all_data(all_data):
    assert print_cell_data(all_data.cell, options.Nlayer, options.Nfrost, 0) \
        is None


def test_print_dmy(dmy):
    assert print_dmy(dmy) is None


def test_print_energy_bal(all_data):
    assert print_energy_bal(all_data.energy, 1, 1) is None


def test_print_filenames():
    assert print_filenames(ctypes.byref(filenames)) is None


def test_print_filep():
    assert print_filep(ctypes.byref(filep)) is None


def test_print_global_param():
    assert print_global_param(ctypes.byref(global_param)) is None


def test_print_option():
    assert print_option(ctypes.byref(options)) is None


def test_print_parameters():
    assert print_parameters(ctypes.byref(param)) is None


def test_print_snow_data(all_data):
    assert print_snow_data(all_data.snow) is None
