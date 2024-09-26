Contributing
================

i-PI is an open source project, and everyone is welcome to contribute
with bug fixes, documentation, examples, new features, regression tests, etc.

Your contribution should be based on the master branch. We kindly ask you to first fork the project,
make your changes, make sure you comply with all items in our checklist below, and finally create a pull request (PR).

Checklist to create a pull request:

- The PR follows our format compliance (based on `black` and `flake8` as explained above)
- All the new classes and functions include the corresponding docstrings

(If the PR adds a new functionally, please fulfill the next two requirements as well)

- Add a working example to the `examples` folder to showcase the new functionality
- Add a regression test to the `i-pi/ipi_tests/regression_tests` folder (see the corresponding README file for further details)
- Make sure that all the automatic checks pass without any error

We are looking forward to your contribution!

Format Compliance
-----------------

i-PI code should be compliant to a minimal subset of PEP-8 recommendations.
Currently, we require the use of `black` as formatter and linter.
We also ask for the usage of `flake8` for syntactic checks, which is also
part of linting.
In most systems, both packages can be easily installed using `pip`.
BEFORE proceeding to a pull request, the minimal requirement is that you run

```
$ make lint
$ make pretty
```

This will ensure the formatting and linting requirement are applied in the whole
directory tree. Please resolve any warnings or errors that may appear. Your
commit will not pass the CI tests otherwise.

For a more flexible setup, we also provide the script `i-pi-style`, for
which instructions can be obtained by typing

```
$ i-pi-style -h
```