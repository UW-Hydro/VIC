from vic import lib as vic_lib


def test_str_to_bool():
    for s, expected in [('TRUE', True), ('FALSE', False),
                        ('true', True), ('false', False),
                        ('TrUe', True), ('FaLsE', False)]:
        assert vic_lib.str_to_bool(s) == expected


def test_str_to_agg_type():
    assert vic_lib.str_to_agg_type('') == vic_lib.AGG_TYPE_DEFAULT
    assert vic_lib.str_to_agg_type('*') == vic_lib.AGG_TYPE_DEFAULT
    for s in ['AGG_TYPE_AVG', 'AGG_TYPE_BEG', 'AGG_TYPE_END',
              'AGG_TYPE_MAX', 'AGG_TYPE_MIN', 'AGG_TYPE_SUM']:
        expected = getattr(vic_lib, s)
        assert vic_lib.str_to_agg_type(s) == expected
        assert vic_lib.str_to_agg_type(s.lower()) == expected


def test_str_to_out_type():
    assert vic_lib.str_to_out_type('') == vic_lib.OUT_TYPE_DEFAULT
    assert vic_lib.str_to_out_type('*') == vic_lib.OUT_TYPE_DEFAULT
    for s in ['OUT_TYPE_USINT', 'OUT_TYPE_SINT', 'OUT_TYPE_FLOAT',
              'OUT_TYPE_DOUBLE']:
        expected = getattr(vic_lib, s)
        assert vic_lib.str_to_out_type(s) == expected
        assert vic_lib.str_to_out_type(s.lower()) == expected


def test_str_to_out_mult():
    assert vic_lib.str_to_out_mult('') == vic_lib.OUT_MULT_DEFAULT
    assert vic_lib.str_to_out_mult('*') == vic_lib.OUT_MULT_DEFAULT
    for expected in range(0., 10000., 100.):
        assert vic_lib.str_to_out_mult(str(expected)) == float(expected)


def test_str_to_freq_flag():
    for s in ['NEVER', 'NSTEPS', 'NSECONDS', 'NMINUTES', 'NHOURS', 'NDAYS',
              'NMONTHS', 'NYEARS', 'DATE', 'END']:
        expected = getattr(vic_lib, 'FREQ_{}'.format(s))
        assert vic_lib.str_to_freq_flag(s) == expected
        assert vic_lib.str_to_freq_flag(s.lower()) == expected


def test_str_to_ascii_format():
    # TODO: figure out the best way to pass a mutable string to ffi
    pass


def test_str_to_calendar():
    for s in ['STANDARD', 'GREGORIAN', 'PROLEPTIC_GREGORIAN', 'NOLEAP',
              'NO_LEAP', '365_DAY', '360_DAY', 'JULIAN', 'ALL_LEAP', '366_DAY']:
        expected = getattr(vic_lib, 'CALENDAR_{}'.format(s))
        assert vic_lib.str_to_calendar(s) == expected
        assert vic_lib.str_to_calendar(s.lower()) == expected


def test_str_to_timeunits():
    for s in ['SECONDS', 'MINUTES', 'HOURS', 'DAYS']:
        expected = getattr(vic_lib, 'TIME_UNITS_{}'.format(s))
        assert vic_lib.str_to_timeunits(s) == expected
        assert vic_lib.str_to_timeunits(s.lower()) == expected


def test_str_from_time_units():
    # TODO: figure out the best way to pass a mutable string to ffi
    pass


def test_str_from_calendar():
    # TODO: figure out the best way to pass a mutable string to ffi
    pass


def test_cell_method_from_agg_type():
    # TODO: figure out the best way to pass a mutable string to ffi
    pass
