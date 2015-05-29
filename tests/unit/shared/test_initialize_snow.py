import pytest
from vic.vic import make_snow_data, initialize_snow


@pytest.fixture()
def veg_num(scope='function'):
    return 5


@pytest.fixture()
def snow(veg_num, scope='function'):
    return make_snow_data(veg_num)


# def test_initialize_snow(snow, veg_num):
#     assert initialize_snow(snow, veg_num) is None
#     assert snow[0][0].albedo == 0.


# def test_initialize_snow_bools(snow, veg_num):
#     assert initialize_snow(snow, veg_num) is None
#     assert not snow[0][0].MELTING
