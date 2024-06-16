"""Classes that deal with output of simulation data.

Holds classes to deal with the output of different properties, trajectories
and the restart files.
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.

import os

import numpy as np

from ipi.utils.messages import verbosity, info, warning
from ipi.utils.units import unit_to_user
from ipi.utils.softexit import softexit
from ipi.utils.depend import *
import ipi.utils.io as io
from ipi.utils.io.inputs.io_xml import *
from ipi.utils.io import open_backup
from ipi.engine.properties import getkey
from ipi.engine.atoms import *
from ipi.engine.cell import *

__all__ = [
    "PropertyOutput",
    "TrajectoryOutput",
    "CheckpointOutput",
    "OutputList",
    "OutputMaker",
    "BaseOutput",
]


class OutputList(list):
    """A simple decorated list to save the output prefix and bring it
    back to the initialization phase of the simulation"""

    def __init__(self, prefix, olist):
        super(OutputList, self).__init__(olist)
        self.prefix = prefix


class OutputMaker:
    """Class to create floating outputs with an appropriate prefix"""

    def __init__(self, prefix="", f_start=False):
        self.prefix = prefix
        self.f_start = f_start

    def bind(self, system):
        self.system = system

    def get_output(self, filename="out", mode=None):
        if self.prefix != "":
            filename = self.prefix + "." + filename
        rout = BaseOutput(filename)
        if mode is None:
            if self.f_start:
                mode = "w"
            else:
                mode = "a"
        rout.bind(mode)
        return rout


class BaseOutput(object):
    """Base class for outputs. Deals with flushing upon close and little more"""

    def __init__(self, filename="out", stride=1):
        """Initializes the class"""

        self.stride = stride
        self.filename = filename
        self.out = None

    def softexit(self):
        """Emergency call when i-pi must exit quickly"""

        self.close_stream()

    def close_stream(self):
        """Closes the output stream"""

        if self.out is not None:
            self.out.close()

    def open_stream(self, mode="w"):
        """Opens the output stream"""

        # Only open a new file if this is a new run, otherwise append.
        self.mode = mode
        self.out = open_backup(self.filename, self.mode)

    def bind(self, mode="w", system=None):
        """Stores a reference to system and registers for exiting"""

        self.system = system
        self.open_stream(mode)
        softexit.register_function(self.softexit)

    def force_flush(self):
        """Tries hard to flush the output stream"""

        if self.out is not None:
            self.out.flush()
            os.fsync(self.out)

    def remove(self):
        """Removes (temporary) output"""

        if self.out is not None:
            self.out.close()
            os.remove(self.filename)

    def write(self, data):
        """Writes data to file"""

        if self.out is not None:
            return self.out.write(data)

    def active(self):
        """Whether we will output at this step"""

        return (self.system.simul.step + 1) % self.stride == 0


class PropertyOutput(BaseOutput):
    """Class dealing with outputting a set of properties to file.

    Does not do any calculation, just manages opening a file, getting data
    from a Properties object and outputting with the desired stride.

    Attributes:
       filename: The name of the file to output to.
       outlist: A list of the properties to be output.
       stride: The number of steps that should be taken between outputting the
          data to file.
       flush: How often we should flush to disk.
       nout: Number of steps since data was last flushed.
       out: The output stream on which to output the properties.
       system: The system object to get the data to be output from.
    """

    def __init__(self, filename="out", stride=1, flush=1, outlist=None):
        """Initializes a property output stream opening the corresponding
        file name.

        Also writes out headers.

        Args:
           filename: A string giving the name of the file to be output to.
           stride: An integer giving how many steps should be taken between
              outputting the data to file.
           flush: Number of writes to file between flushing data.
           outlist: A list of all the properties that should be output.
        """

        super(PropertyOutput, self).__init__(filename, stride)

        if outlist is None:
            outlist = np.zeros(0, np.dtype("|U1024"))
        self.outlist = np.asarray(outlist, np.dtype("|U1024"))
        self.flush = flush
        self.nout = 0

    def bind(self, system, mode="w"):
        """Binds output proxy to System object.

        Args:
           system: A System object to be bound.
        """

        # Checks as soon as possible if some asked-for properties are
        # missing or mispelled

        for what in self.outlist:
            key = getkey(what)
            if key not in list(system.properties.property_dict.keys()):
                print(
                    "Computable properties list: ",
                    list(system.properties.property_dict.keys()),
                )
                raise KeyError(key + " is not a recognized property")

        super(PropertyOutput, self).bind(mode, system)

    def print_header(self):
        # print nice header if information is available on the properties
        icol = 1
        for what in self.outlist:
            ohead = "# "
            key = getkey(what)
            prop = self.system.properties.property_dict[key]

            if "size" not in prop or prop["size"] == 1:
                ohead += "column %3d    " % (icol)
                icol += 1
            else:
                if (type(prop["size"]) is str) or (prop["size"] <= 0):
                    raise RuntimeError("ERROR: property %s has undefined size." % key)
                elif prop["size"] > 1:
                    ohead += "cols.  %3d-%-3d" % (icol, icol + prop["size"] - 1)
                    icol += prop["size"]

            ohead += " --> %s " % (what)
            if "help" in prop:
                ohead += ": " + prop["help"]
            self.out.write(ohead + "\n")

    def write(self):
        """Outputs the required properties of the system.

        Note that properties are outputted using the same format as for the
        output to the xml checkpoint files, as specified in io_xml.

        Raises:
           KeyError: Raised if one of the properties specified in the output list
              are not contained in the property_dict member of properties.
        """

        if softexit.triggered:
            return  # don't write if we are about to exit!

        if not self.active():
            return
        self.out.write("  ")
        for what in self.outlist:
            try:
                quantity, dimension, unit = self.system.properties[what]
                if dimension != "" and unit != "":
                    quantity = unit_to_user(dimension, unit, quantity)
            except KeyError:
                raise KeyError(what + " is not a recognized property")
            if not hasattr(quantity, "__len__"):
                self.out.write(write_type(float, quantity) + "   ")
            else:
                for el in quantity:
                    self.out.write(write_type(float, el) + " ")

        self.out.write("\n")

        self.nout += 1
        if self.flush > 0 and self.nout >= self.flush:
            self.force_flush()
            self.nout = 0


class TrajectoryOutput(BaseOutput):
    """Class dealing with outputting atom-based properties as a
    trajectory file.

    Does not do any calculation, just manages opening a file, getting data
    from a Trajectories object and outputting with the desired stride.

    Attributes:
       filename: The (base) name of the file to output to.
       format: The format of the trajectory file to be created.
       what: The trajectory that needs to be output.
       stride: The number of steps that should be taken between outputting the
          data to file.
       out: The output stream on which to output the trajectories.
       flush: How often we should flush to disk.
       nout: Number of steps since data was last flushed.
       ibead: Index of the replica to print the trajectory of.
       cell_units: The units that the cell parameters are given in.
       system: The System object to get the data to be output from.
    """

    def __init__(
        self,
        filename="out",
        stride=1,
        flush=1,
        what="",
        format="xyz",
        cell_units="atomic_unit",
        ibead=-1,
        extra_type="raw",
    ):
        """Initializes a trajectory output stream opening the corresponding
        file name.

        Also writes out headers.

        Args:
           filename: A string giving the name of the file to be output to.
           stride: An integer giving how many steps should be taken between
              outputting the data to file.
           flush: How often we should flush to disk
           what: A string specifying what trajectory should be output.
           format: A string specifying the type of trajectory file to be created.
           cell_units: A string specifying the units that the cell parameters are
              given in.
           ibead: If positive, prints out only the selected bead. If negative, prints out one file per bead.
           extra_type: Specifies the type of extras string that is printed in the file
        """

        super(TrajectoryOutput, self).__init__(filename, stride)
        self.what = what
        self.flush = flush
        self.ibead = ibead
        self.format = format
        self.cell_units = cell_units
        self.out = None
        self.nout = 0
        self.extra_type = extra_type

    def bind(self, system, mode="w"):
        """Binds output proxy to System object.

        Args:
           system: A System object to be bound.
        """

        # Checks as soon as possible if some asked-for trajs are missing or misspelled
        key = getkey(self.what)
        if key not in list(system.trajs.traj_dict.keys()):
            print(
                "Computable trajectories list: ",
                list(system.trajs.traj_dict.keys()),
            )
            raise KeyError(key + " is not a recognized output trajectory")

        super(TrajectoryOutput, self).bind(mode, system)

    def print_header(self):
        """No headers for trajectory files"""
        pass

    def open_stream(self, mode):
        """Opens the output stream(s)."""

        # prepare format string for zero-padded number of beads,
        # including underscore
        fmt_bead = (
            "{0:0"
            + str(int(1 + np.floor(np.log(self.system.beads.nbeads) / np.log(10))))
            + "d}"
        )

        if getkey(self.what) in [
            "positions",
            "velocities",
            "forces",
            "forces_spring",
            "Eforces",
            "extras",
            # "extras_component_raw", write out a single file as we don't know how to do contraction here
            "extras_bias",
            "forces_sc",
            "momenta",
            "becx",
            "becy",
            "becz",
        ]:
            # must write out trajectories for each bead, so must create b streams

            # prepare format string for file name
            if getkey(self.what)[:6] == "extras":
                fmt_fn = self.filename + "_" + fmt_bead
            elif self.format == "ase":
                fmt_fn = self.filename + "_" + fmt_bead + ".extxyz"
            else:
                fmt_fn = self.filename + "_" + fmt_bead + "." + self.format

            # open all files
            self.out = []
            for b in range(self.system.beads.nbeads):
                if (self.ibead < 0 and (b % (-self.ibead) == 0)) or (self.ibead == b):
                    self.out.append(open_backup(fmt_fn.format(b), mode))
                else:
                    # Create null outputs if a single bead output is chosen.
                    self.out.append(None)
        else:
            # open one file
            filename = self.filename
            # prepare format string for file name
            if getkey(self.what)[:6] != "extras":
                if self.format == "ase":
                    filename += ".extxyz"
                else:
                    filename += "." + self.format
            self.out = open_backup(filename, mode)

    def close_stream(self):
        """Closes the output stream."""

        try:
            if hasattr(self.out, "__getitem__"):
                for o in self.out:
                    if o is not None:
                        o.close()
            else:
                self.out.close()
        except AttributeError:
            # This gets called on softexit. We want to carry on to shut down as cleanly as possible
            warning(
                "Exception while closing output stream " + str(self.out), verbosity.low
            )

    def write(self):
        """Writes out the required trajectories."""

        if softexit.triggered:
            return  # don't write if we are about to exit!
        if not self.active():
            return

        doflush = False
        self.nout += 1
        if self.flush > 0 and self.nout >= self.flush:
            doflush = True
            self.nout = 0

        data, dimension, units = self.system.trajs[
            self.what
        ]  # gets the trajectory data that must be printed
        # quick-and-dirty way to check if a trajectory is "global" or per-bead
        # Checks to see if there is a list of files or just a single file.
        if hasattr(self.out, "__getitem__"):
            if self.ibead < 0:
                for b in range(len(self.out)):
                    if self.out[b] is not None:
                        self.write_traj(
                            data,
                            self.what,
                            self.out[b],
                            b,
                            format=self.format,
                            dimension=dimension,
                            units=units,
                            cell_units=self.cell_units,
                            flush=doflush,
                        )
            elif self.ibead < len(self.out):
                self.write_traj(
                    data,
                    self.what,
                    self.out[self.ibead],
                    self.ibead,
                    format=self.format,
                    dimension=dimension,
                    units=units,
                    cell_units=self.cell_units,
                    flush=doflush,
                )
            else:
                raise ValueError(
                    "Selected bead index "
                    + str(self.ibead)
                    + " does not exist for trajectory "
                    + self.what
                )
        else:
            self.write_traj(
                data,
                getkey(self.what),
                self.out,
                b=0,
                format=self.format,
                dimension=dimension,
                units=units,
                cell_units=self.cell_units,
                flush=doflush,
            )

    def write_traj(
        self,
        data,
        what,
        stream,
        b=0,
        format="xyz",
        dimension="",
        units="automatic",
        cell_units="automatic",
        flush=True,
    ):
        """Prints out a frame of a trajectory for the specified quantity and bead.

        Args:
           what: A string specifying what to print.
           b: The bead index. Defaults to 0.
           stream: A reference to the stream on which data will be printed.
           format: The output file format.
           cell_units: The units used to specify the cell parameters.
           flush: A boolean which specifies whether to flush the output buffer
              after each write to file or not.
        """

        key = getkey(what)
        if key in ["extras", "extras_component_raw", "extras_bias"]:
            if key == "extras_component_raw":
                stream.write(
                    " #%s(%s)# Step:  %10d  Bead:  %5d  \n"
                    % (key.upper(), self.extra_type, self.system.simul.step + 1, b)
                )
            else:
                stream.write(
                    " #%s(%s)# Step:  %10d \n"
                    % (key.upper(), self.extra_type, self.system.simul.step + 1)
                )
            if self.extra_type in data:
                try:
                    if key == "extras_component_raw":
                        # don't partition into beads as there might be a different number when contracting
                        floatarray = np.asarray(
                            data[self.extra_type], dtype=float
                        ).squeeze()
                    else:
                        # picks up the desired bead
                        floatarray = np.asarray(data[self.extra_type][b], dtype=float)
                    if floatarray.ndim == 2:
                        stream.write(
                            "\n".join(
                                [
                                    "      ".join(
                                        ["{:15.8f}".format(item) for item in row]
                                    )
                                    for row in floatarray
                                ]
                            )
                        )
                    elif floatarray.ndim == 1:
                        stream.write("      ".join("%15.8f" % el for el in floatarray))
                    else:
                        raise ValueError(
                            "No specialized writer for arrays of dimension > 2"
                        )
                except:
                    stream.write("%s" % data[self.extra_type][b])
                stream.write("\n")
            elif self.extra_type == "raw":
                stream.write(str(data))
                stream.write("\n")
            else:
                raise KeyError(
                    "Extra type '"
                    + self.extra_type
                    + "' is not among the quantities returned by any of the forcefields."
                )
            if flush:
                stream.flush()
                os.fsync(stream)
            return
        elif getkey(what) in [
            "positions",
            "velocities",
            "forces",
            "forces_spring",
            "Eforces",
            "forces_sc",
            "momenta",
            "becx",
            "becy",
            "becz",
        ]:
            fatom = Atoms(self.system.beads.natoms)
            fatom.names[:] = self.system.beads.names
            fatom.q[:] = data[b]
        else:
            fatom = Atoms(self.system.beads.natoms)
            fatom.names[:] = self.system.beads.names
            fatom.q[:] = data

        fcell = Cell()
        fcell.h = self.system.cell.h

        if units == "":
            units = "automatic"
        if cell_units == "":
            cell_units = "automatic"
        io.print_file(
            format,
            fatom,
            fcell,
            stream,
            title=("Step:  %10d  Bead:   %5d " % (self.system.simul.step + 1, b)),
            key=key,
            dimension=dimension,
            units=units,
            cell_units=cell_units,
        )
        if flush:
            stream.flush()
            os.fsync(stream)


class CheckpointOutput:
    """Class dealing with outputting checkpoints.

    Saves the complete status of the simulation at regular intervals.

    Attributes:
       filename: The (base) name of the file to output to.
       step: the number of times a checkpoint has been written out.
       stride: The number of steps that should be taken between outputting the
          data to file.
       overwrite: If True, the checkpoint file is overwritten at each output.
          If False, will output to 'filename_step'. Note that no check is done
          on whether 'filename_step' exists already.
       simul: The simulation object to get the data to be output from.
       status: An input simulation object used to write out the checkpoint file.
    """

    def __init__(self, filename="restart", stride=1000, overwrite=True, step=0):
        """Initializes a checkpoint output proxy.

        Args:
           filename: A string giving the name of the file to be output to.
           stride: An integer giving how many steps should be taken between
              outputting the data to file.
           overwrite: If True, the checkpoint file is overwritten at each output.
              If False, will output to 'filename_step'. Note that no check is done
              on whether 'filename_step' exists already.
           step: The number of checkpoint files that have been created so far.
        """

        self.filename = filename
        self.stride = stride
        self._step = depend_value(name="step", value=step)
        self.overwrite = overwrite
        self._storing = False
        self._continued = False

    def bind(self, simul):
        """Binds output proxy to simulation object.

        Args:
           simul: A simulation object to be bound.
        """

        self.simul = simul
        import ipi.inputs.simulation as isimulation

        self.status = isimulation.InputSimulation()
        self.status.store(simul)

    def active(self):
        """Whether we will output at this step"""

        return (self.simul.step + 1) % self.stride == 0

    def store(self):
        """Stores the current simulation status.

        Used so that, if halfway through a step a kill signal is received,
        we can output a checkpoint file corresponding to the beginning of the
        current step, which is the last time that both the velocities and
        positions would have been consistent.
        """

        self._storing = True
        self.status.store(self.simul)
        self._storing = False

    def write(self, store=True):
        """Writes out the required trajectories.

        Used for both the checkpoint files and the soft-exit restart file.
        We have slightly different behaviour for these two different types of
        checkpoint file, as the soft-exit files have their store() function
        called automatically, and we do not want this to be updated as the
        status of the simulation after a soft-exit call is unlikely to be in
        a consistent state. On the other hand, the standard checkpoint files
        are not automatically updated in this way, and we must manually store the
        current state of the system before writing them.

        Args:
           store: A boolean saying whether the state of the system should be
              stored before writing the checkpoint file.
        """

        if self._storing:
            info(
                "@ CHECKPOINT: Write called while storing. Force re-storing",
                verbosity.low,
            )
            self.store()

        if not self.active():
            return

        # function to use to open files
        open_function = open_backup

        if self.overwrite:
            filename = self.filename
            if self._continued:
                open_function = open
        else:
            filename = self.filename + "_" + str(self.step)

        # Advance the step counter before saving, so next time the correct index will be loaded.
        if store:
            self.step += 1
            self.store()
            self.status.step.store(self.simul.step + 1)

        with open_function(filename, "w") as check_file:
            check_file.write(self.status.write(name="simulation"))

        # Do not use backed up file open on subsequent writes.
        self._continued = True


dproperties(CheckpointOutput, ["step"])
