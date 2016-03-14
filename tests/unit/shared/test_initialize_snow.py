# import pytest
# from vic import lib as vic_lib


# @pytest.fixture()
# def veg_num(scope='function'):
#     return 5


# @pytest.fixture()
# def snow(veg_num, scope='function'):
#     return vic_lib.make_snow_data(veg_num)


# def test_initialize_snow():
#     raise NotImplementedError('problems here in test_initialize_snow.py')
    # veg_num = 5
    # snow = vic_lib.make_snow_data(veg_num)

    # assert vic_lib.initialize_snow(snow, 5) is None
    # assert snow[0][0].albedo == 0.


# def test_initialize_snow_bools(snow, veg_num):
#     assert vic_lib.initialize_snow(snow, veg_num) is None
#     assert not snow[0][0].MELTING
