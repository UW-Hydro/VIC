VIC Test Suite
=======

This is the VIC Test Suite.  There are six main test types:

1.  **unit**:  function level tests.  ***[Note: None currently implemented.]***
2.  **system**: tests that aim to address model runtime issues.  These tests are generally very quick.
    * configuration errors - tests that address model startup and error checking.
    * restart:  tests that address model state and restart capacity.
    * I/O:  tests that address model input and output functionality.
        *  forcings come out the way the come in
        *  parameter files are appropriately read in allocated.
3.  **science**:  tests that aim to assess the model's scientific skill.  Many of these tests are compared to observations of some kind.
4.  **examples**:  a set of examples that users may download and run.
5.  **release**:  longer, full domain simulations performed prior to release demonstrating model output for a final release.

-------

The VIC Test Suite was developed by Joe Hamman at the University of Washington Land Surface Hydrology Group.  For questions about the development or use of VIC, please email the VIC users list serve at vic_users@u.washington.edu.
