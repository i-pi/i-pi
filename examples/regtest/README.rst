========================
Regression tests library
========================

Regression tests allow to compare the output of a program with its previous runs.
In this way, we can detect if our changes broke the code.

Running regression tests
========================
1. Create (if you have not done it before) a reference output.
`i-pi-regtest --create-reference`
2. Introduce changes to the code which you would like to test.
3. Run the tests. The script will compare the current results with the reference output.
`i-pi-regtest`

Creating a regression test case
===============================
In order to create a regression test case, modify the i-PI input: add a regtest header. Also, do not forget to change the number of steps in the i-PI input.

Input must contain a header that specifies which
files are needed to run the test case (dependencies) and what command should be used to run the
driver. You must follow the exact syntax (only whitespace does not matter).

Here is the example of the header:

<!--REGTEST
DEPENDENCIES  h5o2.dms4B.coeff.com.dat h5o2.pes4B.coeff.dat h5o2+.xyz
COMMAND(10)    i-pi-driver -u -h REGTEST_SOCKET -m zundel
COMMAND       i-pi-driver -u -h REGTEST_SOCKET -m zundel
ENDREGTEST-->

`DEPENDECIES`: files needed to run the test (except the i-PI input)
`COMMAND(n)`: command used to run the driver. The `n` is the number of
    instances of the COMMAND to be created. In the example there will be a
    total of 11 instances of the driver running at once.

This script will replace socket adresses. In case the name of the socket
matches more than once in the COMMAND instruction above, only the first one will be
changed. For example: using the command "i-pi-driver -u -m zundel -h zundel"
would only replace the first string "zundel" and cause an error.
