import pytest
import tempfile
from vic import lib as vic_lib


@pytest.fixture()
def param_file(scope='function'):
    p = 'GAUGE_HEIGHT 12.33\n'
    temp = tempfile.NamedTemporaryFile(prefix='test_param', suffix='txt')
    with open(temp.name, 'w') as f:
        f.write(p)
    return vic_lib.open_file(temp.name.encode(), b'r')


def test_get_parameters(param_file):
    assert vic_lib.get_parameters(param_file) is None
    assert vic_lib.param.GAUGE_HEIGHT == 12.33


def test_validate_parameters():
    assert vic_lib.validate_parameters() is None


# The System Exit Calls in vic_log.h are exiting out of the python testing
# session. This test should pass and does to some extent but we need to rethink
# the error call in vic_log.h before we can uncomment the lines below.
# def test_validate_parameters_raise():
#     vic_lib.param.LAPSE_RATE = 100
#     with pytest.raises(SystemExit):
#         vic_lib.validate_parameters()
