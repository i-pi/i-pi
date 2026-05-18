"""Guards the backward-compatibility shims for ipi.utils.scripting and
ipi.utils.parsing, which re-export symbols from their new home in
ipi.scripting and emit a DeprecationWarning on import.
"""

import importlib
import warnings

import pytest


def test_utils_scripting_shim_emits_deprecation_warning():
    import ipi.utils.scripting as shim

    with pytest.warns(DeprecationWarning, match="ipi.scripting"):
        importlib.reload(shim)


def test_utils_parsing_shim_emits_deprecation_warning():
    import ipi.utils.parsing as shim

    with pytest.warns(DeprecationWarning, match="ipi.scripting"):
        importlib.reload(shim)


def test_utils_scripting_shim_reexports_identical_objects():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        import ipi.utils.scripting as shim
    import ipi.scripting as new

    for name in (
        "simulation_xml",
        "motion_nvt_xml",
        "forcefield_xml",
        "InteractiveSimulation",
    ):
        assert getattr(shim, name) is getattr(new, name), name


def test_utils_parsing_shim_reexports_identical_objects():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        import ipi.utils.parsing as shim
    import ipi.scripting as new

    for name in ("read_output", "read_trajectory"):
        assert getattr(shim, name) is getattr(new, name), name
