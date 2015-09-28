from vic import lib as vic_lib
from vic import ffi


def test_initialize_filenames():
    assert vic_lib.initialize_filenames() is None
    assert ffi.string(vic_lib.filenames.init_state) == b'MISSING'


def test_initialize_fileps():
    assert vic_lib.initialize_fileps() is None
