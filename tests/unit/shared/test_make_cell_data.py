from vic.vic import make_cell_data


def test_make_cell_data():
    assert make_cell_data(4) is not None
