from vic import lib as vic_lib, VIC_DRIVER

# TODO: Capture stdout (possibly point to LOGDEST?)


def test_print_version():
    assert vic_lib.print_version(VIC_DRIVER) is None


def test_print_license():
    assert vic_lib.print_license() is None


def test_print_usage():
    assert vic_lib.print_usage(b'vic_driver_string') is None
