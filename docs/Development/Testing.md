VIC Testing
========

The VIC Test Suite includes six main test types:

1.  **unit**:  function level tests. These tests use the Python Drivers api and mostly focus on testing the routines in `vic_run` and `drivers/shared_all`.
2.  **system**: tests that aim to address model runtime issues.  These tests are generally very quick. These tests mostly focus on testing the `classic` and `image` driver functions:
    * configuration errors - tests that address model startup and error checking.
    * restart:  tests that address model state and restart capability.
    * I/O:  tests that address model input and output functionality.
        *  forcings come out the way they come in
        *  parameter files are appropriately ingested and parameter values are correctly allocated.
3.  **science**:  tests that aim to assess the model's scientific skill.  Many of these tests are compared to observations of some kind.
4.  **examples**:  a set of examples that users may download and run.
5.  **release**:  longer, full domain simulations performed prior to release demonstrating model output for a final release.

## Test data

The **system** and **examples** tests use the [VIC sample data repository](https://github.com/UW-Hydro/VIC_sample_data). This repository includes short (e.g. 10 days) test setups for the VIC image and classic drivers.

The **release** tests are under development (as of August 2016).

## Running the VIC Test Suite

The VIC test suite uses a set of Python utilities and libraries to execute the series of tests listed above. Brief instructions for installing and running the test suite are shown below:

1. Download and run the Anaconda installer: http://continuum.io/downloads
2. Install the required dependencies:

        # Create a conda virtual environment including the required dependencies
        conda env create -n vic_test_env -f ci/requirements.yml

        # Activate the virtual environment
        source activate vic_test_env

        # Install the Python driver into the virtual environment
        python ./vic/drivers/python/setup.py install

3. Run the test suite:

        # Print the run_tests.py usage
        ./tests/run_tests.py -h

        # Run the test suite for the unit and examples cases for the classic driver
        ./tests/run_tests.py unit examples \
            --classic=vic/drivers/classic/vic_classic.exe \
            --data_dir=${SAMPLES_PATH}/data \
            --examples=./tests/examples/examples.cfg

## Travis

VIC uses the [Travis CI](http://travis-ci.org/) continuous integration system. VIC's build tests on Travis test the compilation of the main VIC drivers using a range of environments:

- *Compilers*: `gcc` and `clang`
- *Platforms*: `Linux` and `OSX`

The Travis tests also run the **unit**, **system**, and **examples** tests described above.
