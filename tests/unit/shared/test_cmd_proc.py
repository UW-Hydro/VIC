from vic.vic import print_version, print_license, print_usage, VIC_DRIVER


def test_print_version():
    assert print_version(VIC_DRIVER) is None


def test_print_license():
    assert print_license() is None


def test_print_usage():
    assert print_usage('vic_driver_string') is None
