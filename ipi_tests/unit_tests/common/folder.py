# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


import os


def local(file=None):
    """Returns local folder of the tests directory.

    Args:
        - file: Append file to the local folder
    """
    if file is None:
        return os.path.dirname(__file__)
    else:
        return os.path.join(os.path.dirname(__file__), file)
