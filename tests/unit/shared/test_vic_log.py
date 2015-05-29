import datetime
import os
from vic.vic import (finalize_logging, get_current_datetime, get_logname,
                     initialize_log, setup_logging, String, filenames)


def test_finalize_logging():
    assert finalize_logging() is None


def test_get_current_datetime():
    dts = String()
    assert len(dts) == 0
    get_current_datetime(dts)
    assert len(dts) == 14

    now = datetime.datetime.now()
    now_string = now.strftime('%Y%m%d')
    assert now_string == dts[:8]


# Segfaulting for an unknown reason
# def test_get_logname():
#     path = String('some/path/')
#     idnum = 3
#     filename = String()
#     print([type(o) for o in (path, idnum, filename)])
#     print(get_logname.argtypes)
#     assert get_logname(path, idnum, filename) is not None
#     print(filename)
#     assert len(filename) > 0

def test_initialize_log():
    assert initialize_log() is None


def test_setup_logging():
    assert setup_logging(1) is None


def test_setup_logging_to_file(tmpdir):
    filenames.log_path = str(tmpdir.mkdir('log')) + '/'
    num = 321
    assert setup_logging(num) is None
    log_file = os.listdir(filenames.log_path)[0]
    assert log_file.startswith('vic.log')
    assert log_file.endswith('{0}.txt'.format(num))
    finalize_logging()
