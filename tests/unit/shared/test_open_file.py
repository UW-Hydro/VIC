import pytest
import tempfile
from vic import lib as vic_lib


@pytest.fixture()
def temp_file(scope='module'):
    temp = tempfile.NamedTemporaryFile(prefix='test_file', suffix='txt',
                                       delete=False)
    return temp.name.encode()


def test_open_file_read(temp_file):
    print(temp_file)
    assert vic_lib.open_file(temp_file, b'r') is not None


def test_open_file_write(temp_file):
    assert vic_lib.open_file(temp_file, b'w') is not None
