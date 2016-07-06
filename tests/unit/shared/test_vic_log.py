import datetime
from vic import lib as vic_lib
from vic import ffi


def test_finalize_logging():
    assert vic_lib.finalize_logging() is None


def test_get_current_datetime():
    dts = ffi.new('char [2048]')
    assert len(ffi.string(dts)) == 0
    vic_lib.get_current_datetime(dts)
    assert len(ffi.string(dts)) == 14
    now = datetime.datetime.now()
    now_string = now.strftime('%Y%m%d')
    assert now_string == ffi.string(dts)[:8].decode()


def test_get_logname():
    path = ffi.new('char [2048]', b'foobar/')
    idnum = 3
    filename = ffi.new('char [2048]')
    print([type(o) for o in (path, idnum, filename)])
    assert vic_lib.get_logname(path, idnum, filename) is None
    print(ffi.string(filename))
    assert len(ffi.string(filename)) > 0


def test_initialize_log():
    assert vic_lib.initialize_log() is None
