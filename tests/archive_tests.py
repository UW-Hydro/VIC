#!/usr/bin/env python
"""VIC command line test archiving interface"""

from .run_tests import epilog

# -------------------------------------------------------------------- #
def main():
    """
    Archive VIC tests
    """


    # Parse arguments
    test_results = OrderedDict()

    parser = argparse.ArgumentParser(description=description, epilog=epilog,
                                     formatter_class=CustomFormatter)

    parser.add_argument("tests", type=str,
                        help="Test sets to run",
                        choices=['all', 'unit', 'system', 'science',
                                 'examples', 'release'],
                        default=['unit', 'system', 'science'], nargs='+')
    parser.add_argument("--output_dir", type=str,
                        help="directory to get test output",
                        default="$WORKDIR/VIC_tests_{0}".format(datetime.datetime.now().strftime("%Y%m%d")))
    parser.add_argument("--data_dir", type=str,
                        help="directory to put test data",
                        default='test_data/VIC.4.1.2')
    args = parser.parse_args()


    return
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
if __name__ == "__main__":
    main()
# -------------------------------------------------------------------- #
