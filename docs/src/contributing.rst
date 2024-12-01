.. _contributing:

Contributing
============

i-PI is an open source project, and everyone is welcome to contribute
with bug fixes, documentation, examples, new features, regression tests, etc.
Any contribution from users or developers is very much welcome! These contributions
may be full new methods/algorithms that you have developed or that you just would like to
add to i-PI, or simple improvements (that would make your life easier) regarding the
code itself, the documentation, or analysis scripts.


Useful resources
----------------

For more information on i-PI in general please consult the following places:

* i-PI homepage: http://ipi-code.org
* i-PI forum: https://groups.google.com/g/ipi-users
* i-PI github repo, including issues page: https://github.com/i-pi/i-pi

What we expect from a contribution
----------------------------------

To start contributing, note the following general guidelines:

* Your contribution should be based on the main branch of the repository

* We expect you to fork the repository and make all the changes in your fork

* When done, you can open a pull request (PR) to our repository, which will be reviewed.

For a PR to be issued, we will expect that:

* You have produced code that is compliant with our formatting requirements (`black` and `flake8`, as explained in the README here: https://github.com/i-pi/i-pi)

* You have written understandable docstrings for all functions and classes you have created

* For all contributions (bug fixes, new functionalities, etc.), ensure that the continuous integration tests all pass.

* In case of a new functionality:
 
    *  You have added an example to the `examples` folder that showcases it. The inputs should be a simple but meaningful example of the functionality.

    *  You have added a regression test to `i-pi/ipi_tests/regression_tests` folder (more on this below).

Creating a regression test
--------------------------

If you code up a new functionality, it is best that it goes directly to the regression test infrastructure, so that we can ensure that it does not
break with future code restructuring or the addition of new features in the future. The information below can also be found in the README
file in the `i-pi/ipi_tests/regression_tests` folder:

To add a new regression test please provide:

   *  input.xml (and all the files required by it)

   * test_settings.dat 
     This file specifies the driver_model and flags that will be used when running
     the driver to perform the regression test.
     For examples including the usage of option flags, please see:

     tests/NVE/NVE_1/harmonic_python/test_settings.dat
     tests/NVE/NVE_1/harmonic/test_settings.dat

     If your regtest uses a new driver mode, please add it to the ``fortran_driver_models`` list in the `i-pi/ipi_tests/test_tools.py` file.

   * file_to_check.txt specifying the files serving as reference with their
     respective filenames. For an existing example, please see:

     tests/INSTANTON/100K/files_to_check.txt

   * reference files that are listed in the files_to_check.txt file. Currently,
     the available formats are .xyz and numpy-accessible which should be
     specified with the keywords 'xyz' and 'numpy'.

Important points to note:

* In the case of multiple drivers needed for a simulation, the `test_settings.dat` file has to include a line for each driver by specifying as many lines as needed of `driver_model xxxx` in the order they should be assigned to the sockets. That is followed by lines decribing the corresponding flags of each driver.  
   
* The extension `*.out` appears in the `.gitignore` file. This means that the 'ref_simulation.out' file has to be added manually to your commit. You can do that by typing:
  "git add -f ref_simulation.out"

If you are developing something new and you are unsure about how to proceed, we encourage opening an issue with your question. Depending on the development, we can open a channel of discussion with the core developers (Slack) in order to see how to best handle your case!

Adding features involving a client code
---------------------------------------

If your new development in i-PI is directly related to a specific client code or if there is the need of coding a new native client in i-PI there can be some further adjustments to do.
We very much welcome new interfaces and we will be happy to answer your questions. 

If you want to enable the communication of a new client code with i-PI, it is not difficult: Please check an example of how it was  done in ``fortran`` and ``python`` in the `drivers` folder in the repository.
It is especially simple to add a new potential energy that is evaluated in Python: it is sufficient to add a file in the `ipi/pes` folder, specifying 
`__DRIVER_NAME__` (a string that will be used to refer to the PES from the i-PI input or the command line) and `__DRIVER_CLASS__`, the name of the 
actual class, that should provide, directly or through inheritance, a `__call__(self, cell, pos)` function and return a tuple with
`(potential, forces, virial, extras)`. See any of the existing PES files to use as templates - it is particularly simple to create a class that 
piggybacs on an existing ASE-style calculator.


Getting recognition for your contribution
-----------------------------------------

We recognize that adding new features or enhancements to any project is a task that needs more recognition.
If you want your contribution advertised, we are happy to do so through our webpage. Just let us know and we
will add your name, a short description, and links to any publications here_.

.. _here: http://ipi-code.org/about/features/ 


Contributing with a bug report
------------------------------

Bug reports are a very welcome feedback! These should be posted in the i-PI Github repository_ (issues tab). Please include inputs and outputs that allow us to reproduce the bug, and a short description of the bug you are encountering. 
 
.. _repository: https://github.com/i-pi/i-pi/issues

Questions related to i-PI usage
-------------------------------

Questions related to usage or general questions are best posted in our forum (see link above to google groups). 
If we establish that there is the need to open a Github issue related to the question, we will do so.



