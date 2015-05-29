import pytest
import tempfile
from vic.vic import open_file


@pytest.fixture()
def temp_file(scope='module'):
    temp = tempfile.NamedTemporaryFile(prefix='test_file', suffix='txt',
                                       delete=False)
    return temp.name


def test_open_file_read(temp_file):
    print(temp_file)
    assert open_file(temp_file, 'r') is not None


def test_open_file_write(temp_file):
    assert open_file(temp_file, 'w') is not None
