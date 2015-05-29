import pytest
from vic.vic import initialize_filenames, initialize_fileps, filenames, filep


def test_initialize_filenames():
    assert initialize_filenames() is None
    assert filenames.init_state == 'MISSING'


def test_initialize_fileps():
    assert initialize_fileps() is None
    with pytest.raises(ValueError):
        filep.globalparam.contents
