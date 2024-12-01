"""Contains different implementations for reading/checkpointing an i-PI
simulation. For now only xml, but in future possibly also yml/json.
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


__all__ = ["io_xml"]


def read_value(s):
    """Attempt to parse a string to int or float; fallback to string."""
    s = s.strip()
    for cast in (int, float):
        try:
            return cast(s)
        except ValueError:
            continue
    if s.lower() == "false":
        return False
    if s.lower() == "true":
        return True
    return s


def read_args_kwargs(input_str):
    """
    Parses a string into positional arguments and keyword arguments.

    Args:
        input_str (str): The input string containing comma-separated values and key-value pairs.

    Returns:
        tuple: A tuple containing a list of positional arguments and a dictionary of keyword arguments.
    """
    args = []
    kwargs = {}
    tokens = input_str.split(",")
    for token in tokens:
        token = token.strip()
        if "=" in token:
            key, value = token.split("=", 1)
            kwargs[key.strip()] = read_value(value)
        elif len(token) > 0:
            args.append(read_value(token))
    return args, kwargs
