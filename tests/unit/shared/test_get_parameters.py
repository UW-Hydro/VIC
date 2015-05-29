import pytest
import tempfile
from vic.vic import open_file, param, get_parameters, validate_parameters


@pytest.fixture()
def param_file(scope='function'):
    p = 'GAUGE_HEIGHT 12.33\n'
    temp = tempfile.NamedTemporaryFile(prefix='test_param', suffix='txt')
    with open(temp.name, 'w') as f:
        f.write(p)
    return open_file(temp.name, 'r')


def test_get_parameters(param_file):
    assert get_parameters(param_file) is None
    assert param.GAUGE_HEIGHT == 12.33


def test_validate_parameters():
    assert validate_parameters(param) is None


# The System Exit Calls in vic_log.h are exiting out of the python testing
# session. This test should pass and does to some extent but we need to rethink
# the error call in vic_log.h before we can uncomment the lines below.
# def test_validate_parameters_raise():
#     param.LAPSE_RATE = 100
#     with pytest.raises(SystemExit):
#         validate_parameters(param)
