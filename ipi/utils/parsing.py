"""Backward-compatibility shim: use :mod:`ipi.scripting` instead."""

import warnings

warnings.warn(
    "ipi.utils.parsing has moved to ipi.scripting.parsing "
    "(re-exported from ipi.scripting); "
    "update your imports to `from ipi.scripting import ...`. "
    "This compatibility shim will be removed in a future release.",
    DeprecationWarning,
    stacklevel=2,
)

from ipi.scripting.parsing import *  # noqa: F401,F403
from ipi.scripting.parsing import __all__  # noqa: F401
