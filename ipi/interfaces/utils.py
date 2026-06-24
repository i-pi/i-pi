"""Helpers shared by the driver interfaces (sockets, MPI) and FFDirect."""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.

import json


def parse_extra(mxtra):
    """Decodes the optional 'extra' string returned by a force calculation into
    a dict. Empty/whitespace payloads give an empty dict; otherwise any JSON
    fields are returned and the literal string is stashed under 'raw'."""

    if not mxtra or mxtra.isspace():
        return {}
    try:
        mxtradict = json.loads(mxtra)
    except Exception:
        # not JSON: still expose the literal string under 'raw'
        mxtradict = {}
    if "raw" in mxtradict:
        raise ValueError(
            "'raw' cannot be used as a field in a JSON-formatted extra string"
        )
    mxtradict["raw"] = mxtra
    return mxtradict
