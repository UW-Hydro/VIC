import pytest
from vic import lib as vic_lib


@pytest.mark.xfail
def test_make_all_vars():
    all_vars = vic_lib.make_all_vars(5)
    assert all([hasattr(all_vars, member) for member in
               ('snow', 'energy', 'veg_var', 'cell')])
    assert hasattr(all_vars.cell[0][0], 'runoff')
