VIC Testing
========

The VIC Test Suite includes six main test types:

1.  **unit**:  function level tests. These tests mostly focus on testing the routines in `vic_run` and `drivers/shared_all`.
2.  **system**: tests that aim to address model runtime issues.  These tests are generally very quick. These tests mostly focus on testing the `classic` and `image` driver functions:
    * configuration errors - tests that address model startup and error checking.
    * restart:  tests that address model state and restart capacity.
    * I/O:  tests that address model input and output functionality.
        *  forcings come out the way the come in
        *  parameter files are appropriately read in allocated.
3.  **science**:  tests that aim to assess the model's scientific skill.  Many of these tests are compared to observations of some kind.
4.  **examples**:  a set of examples that users may download and run.
5.  **release**:  longer, full domain simulations performed prior to release demonstrating model output for a final release.

## Travis

VIC uses the [Travis](http://travis-ci.org/) continuous integration system. VIC's build tests on Travis test the compilation of the main VIC drivers using a range of environments:

- *Compilers*: `gcc` and `clang`
- *Platforms*: `Linux` and `OSX`

The Travis tests also run the **unit**, **system**, and **examples** tests described above.
