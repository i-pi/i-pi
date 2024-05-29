"""Package with functions for reading and writing files.

This module has machinery for abstract I/O handling.

The idea is that the unit conversion is done here. The default is to guess
units from the file, but it can be overridden.
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


import sys
import os
import io
import json
import numpy as np

from ipi.utils.messages import info, verbosity
from ipi.utils.units import unit_to_user
from ipi.external import importlib
from ipi.utils.decorators import cached

__all__ = [
    "io_units",
    "iter_file",
    "print_file_path",
    "print_file",
    "read_file",
    "netstring_encoded_savez",
    "netstring_encoded_loadz",
    "NumpyEncoder",
]

mode_map = {
    "bin": "binary",
}


io_map = {
    "print_path": "print_%s_path",
    "print": "print_%s",
    "read": "read_%s",
    "iter": "iter_%s",
}


class NumpyEncoder(json.JSONEncoder):
    """Special json encoder for numpy types"""

    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)


@cached
def _get_io_function(mode, io):
    """Returns io function with specified mode.

    This will be determined on the fly based on file name extension and
    available ipi/utils/io/backends/io_*.py backends.

    Args:
        mode: Which format has the file? e.g. "pdb", "xml", "bin", "extxyz", or "xyz"
        io: One of "print_path", "print", "read" or "iter"
    """

    try:
        mode = mode[mode.find(".") + 1 :]
        if mode in mode_map:
            mode = mode_map[mode]
        module = importlib.import_module("ipi.utils.io.backends.io_%s" % mode)
    except ImportError:
        print("Error: mode %s is not supported." % mode)
        sys.exit()

    try:
        func = getattr(module, io_map[io] % mode)
    except KeyError:
        print("Error: io %s is not supported with mode %s." % (io, mode))
        sys.exit()

    return func


def print_file_path_raw(
    mode, beads, cell, filedesc=sys.stdout, title="", cell_conv=1.0, atoms_conv=1.0
):
    """Prints all the bead configurations, into a `mode` formatted file.

    Prints all the replicas for each time step separately, rather than all at
    once.

    Args:
        beads: A beads object giving the bead positions.
        cell: A cell object giving the system box.
        filedesc: An open writable file object. Defaults to standard output.
    """

    return _get_io_function(mode, "print_path")(
        beads=beads,
        cell=cell,
        filedesc=filedesc,
        title=title,
        cell_conv=cell_conv,
        atoms_conv=atoms_conv,
    )


def print_file_path(
    mode,
    beads,
    cell,
    filedesc=sys.stdout,
    title="",
    key="",
    dimension="length",
    units="automatic",
    cell_units="automatic",
):
    """Prints all the bead configurations, into a `mode` formatted file.

    Prints all the replicas for each time step separately, rather than all at
    once.

    Args:
        beads: A beads object giving the bead positions.
        cell: A cell object giving the system box.
        filedesc: An open writable file object. Defaults to standard output.
    """

    if mode == "pdb":  # special case for PDB
        if dimension != "length":
            raise ValueError("PDB Standard is only designed for atomic positions")
        if units == "automatic":
            units = "angstrom"
        if cell_units == "automatic":
            cell_units = "angstrom"

    # if output mode is 'ase', "automatic" units are actually "ase". In so far, only "ase" units can be used.
    # Raises an error even if units like "angstrom" are used. This is because atomic_units to angstrom
    # conversion used in i-PI can differ from ASE conversion factors up to numerics.
    elif mode == "ase":
        if units == "automatic":
            units = "ase"
        elif units != "ase":
            raise ValueError(
                'Only "ase" or default units can be used with extended xyz format.'
            )

        if cell_units == "automatic":
            units = "ase"
        elif units != "ase":
            raise ValueError(
                'Only "ase" or default units can be used with extended xyz format.'
            )

    # in general, "automatic" units are actually "atomic_units"
    else:
        if units == "automatic":
            units = "atomic_unit"
        if cell_units == "automatic":
            cell_units = "atomic_unit"

    cell_conv = unit_to_user("length", cell_units, 1.0)
    atoms_conv = unit_to_user(dimension, units, 1.0)

    title = title + ("%s{%s}  cell{%s}" % (key, units, cell_units))

    return _get_io_function(mode, "print_path")(
        beads=beads,
        cell=cell,
        filedesc=filedesc,
        cell_conv=cell_conv,
        atoms_conv=atoms_conv,
    )


def print_file_raw(
    mode, atoms, cell, filedesc=sys.stdout, title="", cell_conv=1.0, atoms_conv=1.0
):
    """Prints atom positions, or atom-vector properties, into a `mode` formatted file,
       providing atoms and cell in the internal i-PI representation but doing no conversion.

    Args:
        atoms: An atoms object containing the positions (or properties) of the atoms
        cell: A cell object containing the system box.
        filedesc: An open writable file object. Defaults to standard output.
        title: This contains the string that will be used for the comment line.
        cell_conv: Conversion factor for the cell parameters
        atoms_conv: Conversion factors for the atomic properties
    """

    return _get_io_function(mode, "print")(
        atoms=atoms,
        cell=cell,
        filedesc=filedesc,
        title=title,
        cell_conv=cell_conv,
        atoms_conv=atoms_conv,
    )


def print_file(
    mode,
    atoms,
    cell,
    filedesc=sys.stdout,
    title="",
    key="",
    dimension="length",
    units="automatic",
    cell_units="automatic",
):
    """Prints atom positions, or atom-vector properties, into a `mode` formatted file,
       using i-PI internal representation of atoms & cell. Does conversion and prepares
       formatted title line.

    Args:
        mode: I/O file format (e.g. "xyz")
        atoms: An atoms object containing the positions (or properties) of the atoms
        cell: A cell object containing the system box.
        filedesc: An open writable file object. Defaults to standard output.
        title: This gives a string to be appended to the comment line.
        key: Description of the property that is being output
        dimension: Dimensions of the property (e.g. "length")
        units: Units for the output (e.g. "angstrom")
        cell_units: Units for the cell (dimension length, e.g. "angstrom")
    """
    if mode == "pdb":  # special case for PDB and ASE
        if dimension != "length":
            raise ValueError("PDB Standard is only designed for atomic positions")
        if units == "automatic":
            units = "angstrom"
        if cell_units == "automatic":
            cell_units = "angstrom"
        if key == "":
            key = "positions"

    # if output mode is 'ase', "automatic" units are actually "ase". In so far, only "ase" units can be used.
    # Raises an error even if units like "angstrom" are used. This is because atomic_units to angstrom
    # conversion used in i-PI can differ from ASE conversion factors up to numerics.
    elif mode == "ase":
        if units == "automatic":
            units = "ase"
        elif units != "ase":
            raise ValueError(
                'Only "ase" or default units can be used with extended xyz format.'
            )
        if cell_units == "automatic":
            cell_units = "ase"
        elif cell_units != "ase":
            raise ValueError(
                'Only "ase" or default units can be used with extended xyz format.'
            )

    # in general, "automatic" units are actually "atomic_units"
    else:
        if units == "automatic":
            units = "atomic_unit"
        if cell_units == "automatic":
            cell_units = "atomic_unit"

    cell_conv = unit_to_user("length", cell_units, 1.0)
    atoms_conv = unit_to_user(dimension, units, 1.0)

    title = title + ("%s{%s}  cell{%s}" % (key, units, cell_units))

    print_file_raw(
        mode=mode,
        atoms=atoms,
        cell=cell,
        filedesc=filedesc,
        title=title,
        cell_conv=cell_conv,
        atoms_conv=atoms_conv,
    )


def read_file_raw(mode, filedesc):
    """Reads atom positions, or atom-vector properties, from a file of mode "mode",
        returns positions and cell parameters in raw array format, without creating i-PI
        internal objects.

    Args:
        mode: I/O file format (e.g. "xyz")
        filedesc: An open readable file object.

    """
    reader = _get_io_function(mode, "read")

    comment, cell, atoms, names, masses = reader(filedesc=filedesc)

    return {
        "comment": comment,
        "data": atoms,
        "masses": masses,
        "names": names,
        "natoms": len(names),
        "cell": cell,
    }


def read_file(
    mode, filedesc, dimension="automatic", units="automatic", cell_units="automatic"
):
    """Reads one frame from an open `mode`-style file. Also performs units
        conversion as requested, or as guessed from the input comment line.

    Args:
        mode: I/O file format (e.g. "xyz")
        filedesc: An open readable file object from a `mode` formatted file.
        dimension: Dimensions of the property (e.g. "length")
        units: Units for the input (e.g. "angstrom")
        cell_units: Units for the cell (dimension length, e.g. "angstrom")

        All other args are passed directly to the responsible io function.

    Returns:
        A dictionary as returned by `process_units`.
    """

    raw_read = read_file_raw(mode=mode, filedesc=filedesc)

    # late import is needed to break an import cycle
    from .io_units import process_units

    return process_units(
        dimension=dimension, units=units, cell_units=cell_units, mode=mode, **raw_read
    )


def read_file_name(filename):
    """Read one frame from file, guessing its format from the extension.

    Args:
        filename: Name of input file.

    Returns:
        A dictionary with 'atoms', 'cell' and 'comment'.
    """

    return read_file(os.path.splitext(filename)[1], open(filename))


def iter_file_raw(mode, filedesc):
    """Takes an open `mode`-style file and yields a dictionary of positions and cell parameters in raw array format, without creating i-PI internal objects.

    Args:
        filedesc: An open readable file object from a `mode` formatted file.

    Returns:
        Generator of frames dictionaries, as returned by `process_units`.
    """

    reader = _get_io_function(mode, "read")

    try:
        while True:
            comment, cell, atoms, names, masses = reader(filedesc=filedesc)
            yield {
                "comment": comment,
                "data": atoms,
                "masses": masses,
                "names": names,
                "natoms": len(names),
                "cell": cell,
            }
    except EOFError:
        pass


def iter_file(
    mode, filedesc, dimension="automatic", units="automatic", cell_units="automatic"
):
    """Takes an open `mode`-style file and yields one Atoms object after another.

    Args:
        filedesc: An open readable file object from a `mode` formatted file.

    Returns:
        Generator of frames dictionaries, as returned by `process_units`.
    """

    # late import is needed to break an import cycle
    from .io_units import process_units

    iter_file_raw_generator = iter_file_raw(mode, filedesc)
    for raw_read in iter_file_raw_generator:
        yield process_units(
            dimension=dimension,
            units=units,
            cell_units=cell_units,
            mode=mode,
            **raw_read
        )


def iter_file_name(filename):
    """Open a trajectory file, guessing its format from the extension.

    Args:
        filename: Filename of a trajectory file.

    Returns:
        Generator of frames (dictionary with 'atoms', 'cell' and 'comment') from the trajectory in `filename`.
    """

    return iter_file(os.path.splitext(filename)[1], open(filename))


def iter_file_name_raw(filename):
    """Open a trajectory file, guessing its format from the extension.

    Args:
        filename: Filename of a trajectory file.

    Returns:
        Raw  I/O iterator
    """

    return iter_file_raw(mode=os.path.splitext(filename)[1], filedesc=open(filename))


def open_backup(filename, mode="r", buffering=-1):
    """A wrapper around `open` which saves backup files.

    If the file is opened in write mode and already exists, it is first
    backed up under a new file name, keeping all previous backups. Then,
    a new file is opened for writing.

    For reference: https://docs.python.org/2/library/functions.html#open

    Args:
        The same as for `open`.

    Returns:
        An open file as returned by `open`.
    """

    if mode.startswith("w"):
        # If writing, make sure nothing is overwritten.

        i = 0
        fn_backup = filename
        while os.path.isfile(fn_backup):
            fn_backup = "#" + filename + "#%i#" % i
            i += 1

        if fn_backup != filename:
            os.rename(filename, fn_backup)
            info(
                "Backup performed: {0:s} -> {1:s}".format(filename, fn_backup),
                verbosity.low,
            )

    else:
        # There is no need to back up.
        # `open` will sort out whether `mode` is valid.
        pass

    return open(filename, mode, buffering)


# def netstring_encoded_savez(ofile, compressed=True, *unnamed_objs, **named_objs):


def netstring_encoded_savez(ofile, compressed=True, **named_objs):
    output = io.StringIO()
    if compressed:
        # np.savez_compressed(output,*unnamed_objs,**named_objs)
        np.savez_compressed(output, **named_objs)
    else:
        # np.savez(output,*unnamed_objs,**named_objs)
        np.savez(output, **named_objs)
    content = output.getvalue()
    ofile.write(str(len(content)) + ":" + content + ",")


def netstring_encoded_loadz(ifile):
    # read string length
    c = ifile.read(1)
    if c == "0":
        raise ValueError("Reading an empty netstring")
    length = ""
    while c in ["0", "1", "2", "3", "4", "5", "6", "7", "8", "9"]:
        length += c
        c = ifile.read(1)
    if not c == ":":
        raise ValueError("Invalid netstring delimiter")
    content = ifile.read(int(length))
    if not ifile.read(1) == ",":
        raise ValueError("Invalid netstring delimiter")

    istr = io.StringIO(content)
    npz = np.load(istr)
    rdic = {}
    for a in npz.files:
        rdic[a] = npz[a]
    istr.close()

    return rdic
