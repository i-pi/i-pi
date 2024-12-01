-- tests directory --
=====================================================================================================================
 * A the moment we can perform four type of tests:
   
    - unit_test: These tests check basic functionality of the code, like io from/to xyz files is working properly

    - examples: These tests search for all the examples provided in the i-pi/examples folder and check
                   the integrity of the corresponding input files by running i-PI for a few steps with a toy model.

    - regression_test: These test check complex functionality by running i-PI for a few steps and comparing the output obatained
                      against reference results. You can find all the regression tests inside i-pi/ipi_test/regression_tests/tests/ 

    - profiling: These tests are designed to evaluate the computational (overhead) cost of i-PI. 

 * To run the unitary, example and regression tests, please use i-pi/bin/i-pi-tests. For details of usage call it with "-h" option.

 * The tests can be customized to an extent by including a `test_settings.dat` file that specifies how to run the test,
   mostly by describing the driver that should be used. See some of the existing files as examples. 

 * In some situations it is desirable to perform regression/examples tests only for a selected set.
   For these cases, we also provide more flexible scripts for performing the regresion tests and testing the provided examples.
   See i-pi/ipi_tests/regression_tests/test_run.py and i-pi/ipi_tests/examples/test_examples.py, respectively. 
   For details of usage call them with "-h" option.

 * To run the profiling test, have a look a profiling/README.md


