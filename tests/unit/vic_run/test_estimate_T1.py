from vic.vic import estimate_T1


def test_estimate_t1():
    # these values may not be reasonable (test for completeness)
    t1 = estimate_T1(2.3, 2.0, 1.8, 0.2, 0.15, 0.2, 0.2, 0.73, 4., 0.5)
    assert t1 != 0.
