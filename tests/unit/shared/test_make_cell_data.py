from vic import lib as vic_lib


def test_make_cell_data():
    assert vic_lib.make_cell_data(4) is not None
