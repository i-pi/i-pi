try:
    # importlib is a python stdlib from python 2.7 on
    from importlib import *
except ImportError:
    from .bundledimportlib import *
