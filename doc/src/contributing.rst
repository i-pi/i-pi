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

- i-PI forum:
- i-PI github repo, including issues page:
- i-PI homepage:
- i-PI documentation:

What we expect from a contribution
----------------------------------

To start contributing, note the following general guidelines:

- Your contribution should be based on the master branch of the repository
- We expect you to fork the repository and make all the changes in your fork
- When done, you can open a pull request (PR) to our repository, which will be reviewed.

For a PR to be issued, we will expect that:

- You have produced code that is compliant with our formatting requirements (`black` and `flake8`, as explained in the README here: )
- You have written understandable docstrings for all functions and classes you have created
- For all contributions (bug fixes, new functionalities, etc.), ensure that the continuous integration tests all pass.
- In case of a new functionality: 
  + You have added an example to the `examples` folder that showcases it. The inputs should be a simple but meaningful example of the functionality.
  + You have added a regression test to `i-pi/ipi_tests/regression_tests` folder (more on this below).

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

   * file_to_check.txt specifying the files serving as reference with their
     respective filenames

     for an existing example, please see:
        tests/INSTANTON/100K/files_to_check.txt

   * reference files that are listed in the files_to_check.txt file. Currently,
     the available formats are .xyz and numpy-accessible which should be
     specified with the keywords 'xyz' and 'numpy'.

Important points to note:

- In the case of multiple drivers needed for a simulation, the `test_settings.dat` file has to include a line for each driver by specifying as many lines as needed of `driver_model xxxx` in the order they should be assigned to the sockets. That
is followed by lines decribing the corresponding flags of each driver.  
   
- The extension `*.out` appears in the `.gitignore` file. This means that the 'ref_simulation.out' file has to be added manually to your commit.
   You can do that by typing:

    "git add -f ref_simulation.out"


Getting recognition for your contribution
-----------------------------------------

We recognize that adding new features or enhancements to any project is a task that needs more recognition.
If you want your contribution advertised, we are happy to do so through our webpage. Just let us know and we
will add your name, a short description, and links to any publications here_:

.. _here: http://ipi-code.org/about/features/ 



Contributing with a bug report
------------------------------

Bug reports are a very welcome feedback! These should be posted in the i-PI Github repository (issues tab). Please include inputs and outputs that allow us to reproduce the bug, and a short description of the bug you are encountering. 
 

Questions related to i-PI usage
-------------------------------

Questions related to usage or general questions are best posted in our forum (see link above). 
If we establish that there is the need to open a Github issue related to the question, we will do so.



DRAFT!

List of points to be addressed
 
- Welcome contributors
- Try to do a table of contents (need to find out how)
- Link to: homepage, issues/bug reports and forum
- Where are our tests
- What is the PR protocol (add example, etc.)
- How to report bugs and ask about enhancements
- Where to post "how to use i-PI questions"
- How to get your contribution listed
- Details about client codes and how to contribute a feature involving a client




