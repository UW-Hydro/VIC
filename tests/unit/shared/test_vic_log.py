import datetime
import os
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
    assert now_string == dts[:8]


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


def test_setup_logging():
    assert vic_lib.setup_logging(1) is None


def test_setup_logging_to_file(tmpdir):
    vic_lib.filenames.log_path = (str(tmpdir.mkdir('log')) + '/').encode()
    num = 321
    assert vic_lib.setup_logging(num) is None
    log_file = os.listdir(ffi.string(vic_lib.filenames.log_path))[0]
    print(log_file)
    assert log_file.startswith(b'vic.log')
    assert log_file.endswith('{0}.txt'.format(num).encode())
    vic_lib.finalize_logging()
