Troubleshooting
===============

Input errors
------------

-  *not well-formed (invalid token)*: Seen if the input file does not
   have the correct xml syntax. Should be accompanied by a line number
   giving the point in the file where the syntax is incorrect.

-  *mismatched tag*: One of the closing tags does not have the same name
   as the corresponding opening tag. Could be caused either by a
   misspelling of one of the tags, or by having the closing tag in the
   wrong place. This last one is a standard part of the xml syntax, if
   the opening tag of one item is after the opening tag of a second,
   then its closing tag should be before the closing tag of the second.
   Should be accompanied by a line number giving the position of the
   closing tag.

-  *Uninitialized value of type \_\_* or *Attribute/Field name \_\_ is
   mandatory and was not found in the input for property \_\_*: The xml
   file is missing a mandatory tag, i.e. one without which the
   simulation cannot be initialized. Find which tag name is missing and
   add it.

-  *Attribute/tag name \_\_ is not a recognized property of \_\_
   objects*: The first tag should not be found within the second set of
   tags. Check that the first tag is spelt correctly, and that it has
   been put in the right place.

-  *\_\_ is not a valid option (___)*: This attribute/tag only allows a
   certain range of inputs. Pick one of the items from the list given.

-  *\_\_ is an undefined unit for kind \_\_* or *\_\_ is not a valid
   unit prefix* or *Unit \_\_ is not structured with a prefix+base
   syntax*: The unit input by the user is not correct. Make sure it
   corresponds to the correct dimensionality, and is spelt correctly.

-  *Invalid literal for int() with base 10: \_\_* or *Invalid literal
   for float(): \_\_* or *\_\_ does not represent a bool value*: The
   data input by the user does not have the correct data type. See
   section `3.1.1 <#ifilestructure>`__ for what constitutes a valid
   integer/float/boolean value.

-  *Error in list syntax: could not locate delimiters*: The array input
   data did not have the required braces. For a normal array use [], for
   a dictionary use {}, and for a tuple use ().

-  *The number of atom records does not match the header of xyz file*:
   Self-explanatory.

-  *list index out of range*: This will normally occur if the
   configuration is initialized from an invalid input file. This will
   either cause the code to try to read part of the input file that does
   not exist, or to set the number of beads to zero which causes this
   error in a different place. Check that the input file has the correct
   syntax.

Initialization errors
---------------------

-  *Negative \_\_ parameter specified.*: Self-explanatory.

-  *If you are initializing cell from cell side lengths you must pass
   the ’cell’ tag an array of 3 floats*: If you are attempting to
   initialize a cell using the “abc” mode, the code expects three floats
   corresponding to the three side lengths.

-  *If you are initializing cell from cell side lengths and angles you
   must pass the ’cell’ tag an array of 6 floats*: If you are attempting
   to initialize a cell using the “abcABC” mode, the code expects six
   floats corresponding to the three side lengths, followed by the three
   angles in degrees.

-  *Cell objects must contain a 3x3 matrix describing the cell
   vectors.*: If you are attempting to initialize a cell using the
   “manual” mode, the code expects nine floats corresponding to the cell
   vector matrix side lengths. Note that the values of the
   lower-diagonal elements will be set to zero.

-  *Array shape mismatch in q/p/m/names in beads input*: The size of the
   array in question does not have the correct number of elements given
   the number of atoms and the number of beads used in the rest of the
   input. If the number of beads is nbeads and the number of atoms
   natoms, then q and p should have shape (nbeads, 3\*natoms) and m and
   names should have shape (natoms,).

-  *No thermostat/barostat tag provided for NVT/NPT simulation*: Some
   ensembles can only be sampled if a thermostat and/or barostat have
   been defined, and so for simulations at constant temperature and/or
   pressure these tags are mandatory. If you wish to not use a
   thermostat/barostat, but still want to keep the ensemble the same,
   then use “dummy” mode thermostat/barostat, which simply does nothing.

-  *Pressure/Temperature should be supplied for constant
   pressure/temperature simulation*: Since in this case the ensemble is
   defined by these parameters, these must be input by the user. Add the
   appropriate tags to the input file.

-  *Manual path mode requires (nbeads-1) frequencies, one for each
   internal mode of the path.*: If the “mode” tag of “normal_modes” is
   set to “manual”, it will expect an array of frequencies, one for each
   of the internal normal modes of the ring polymers.

-  *PA-CMD mode requires the target frequency of all the internal
   modes.*: If the “mode” tag of “normal_modes” is set to “pa-cmd”, it
   will expect an array of one frequency, to which all the internal
   modes of the ring polymers will be set.

-  *WMAX-CMD mode requires [wmax, wtarget]. The normal modes will be
   scaled such that the first internal mode is at frequency wtarget and
   all the normal modes coincide at frequency wmax.*: If the “mode” tag
   of “normal_modes” is set to “wmax-cmd”, it will expect an array of
   two frequencies, one two set the lowest frequency normal mode, and
   one for the other normal mode frequencies.

-  *Number of beads \_\_ doesn’t match GLE parameter nb= \_\_*: The
   matrices used to define the generalized Langevin equations of motion
   do not have the correct first dimension. If matrices have been
   downloaded from `<http://http://imx-cosmo.github.io/gle4md/>`_
   make sure that you have input the correct number of beads.

-  *Initialization tries to match up structures with different atom
   numbers*: If in the initialization any of the matrices has an array
   shape which do not correspond to the same number of atoms, then they
   cannot correspond to the same system. Check the size of the arrays
   specified if they have been input manually.

-  *Cannot initialize single atom/bead as atom/bead index \_\_ is
   larger than the number of atoms/beads*: Self-explanatory. However,
   note that indices are counted from 0, so the first replica/atom is
   defined by an index 0, the second by an index 1, and so on.

-  *Cannot initialize the momenta/masses/labels/single atoms before the
   size of the system is known.*: In the code, a beads object is created
   to hold all the information related to the configuration of the
   system. However, until a position vector has been defined, this
   object is not created. Therefore, whichever arrays are being
   initialized individually, the position vector must always be
   initialized first.

-  *Trying to resample velocities before having masses.*: A
   Maxwell-Boltzmann distribution is partly defined by the atomic
   masses, and so the masses must be defined before the velocities can
   be resampled from this distribution.

-  *Cannot thermalize a single bead*: It does not make sense to
   initialize the momenta of only one of the beads, and so i-PI does not
   give this functionality.

-  *Initializer could not initialize \_\_*: A property of the system
   that is mandatory to properly run the simulation has not been
   initialized in either the “initialize” section or the appropriate
   section in beads.

-  *Ensemble does not have a thermostat to initialize* or *There is
   nothing to initialize in non-GLE thermostats* or *Checkpoint file
   does not contain usable thermostat data*: These are raised if the
   user has tried to initialize the matrices for the GLE thermostats
   with a checkpoint file that either does not have a GLE thermostat or
   does not have a thermostat at all.

-  *Size mismatch in thermostat initialization data*: Called if the
   shape of the GLE matrices defined in the checkpoint file is different
   from those defined in the new simulation.

-  *Replay can only read from PDB or XYZ files – or a single frame from
   a CHK file*: If the user specifies a replay ensemble, the state of
   the system must be defined by either a configuration file or a
   checkpoint file, and cannot be specified manually.

Output errors
-------------

-  *The stride length for the \_\_ file output must be positive.*:
   Self-explanatory

-  *\_\_ is not a recognized property/output trajectory*: The string as
   defined in the “properties”/”trajectory” tag does not correspond to
   one of the available trajectories. Make sure that both the syntax is
   correct, and that the property has been spelt correctly.

-  *Could not open file \_\_ for output*: Raised if there is a problem
   opening the file defined by the “filename” attribute.

-  *Selected bead index \_\_ does not exist for trajectory \_\_*: You
   have asked for the trajectory of a bead index greater than the number
   of the replicas of the system. Note that indices are counted from 0,
   so the first replica is defined by an index 0, the second by an index
   1, and so on.

-  *Incorrect format in unit specification \_\_*: Usually raised if one
   of the curly braces has been neglected.

-  *Incorrect format in argument list \_\_*: This will be raised either
   if one of the brackets has been neglected, or if the delimiters
   between arguments, in this case “;”, are not correct. This is usually
   raised if, instead of separating the arguments using “;”, they are
   instead separated by “,”, since this causes the property array to be
   parsed incorrectly.

-  *\_\_ got an unexpected keyword argument \_\_*: This will occur if
   one of the argument lists of one of the properties specified by the
   user has a keyword argument that does not match any of those in the
   function to calculate it. Check the properties.py module to see which
   property this function is calculating, and what the correct keyword
   arguments are. Then check the “properties” tag, and find which of the
   arguments has been misspelt.

-  *Must specify the index of atom_vec property*: Any property which
   prints out a vector corresponding to one atom needs the index of that
   atom, as no default is specified.

-  *Cannot output \_\_ as atom/bead index \_\_ is larger than the
   number of atoms/beads*: Self-explanatory. However, note that indices
   are counted from 0, so the first replica/atom is defined by an index
   0, the second by an index 1, and so on.

-  *Couldn’t find an atom that matched the argument of \_\_*: For
   certain properties, you can specify an atom index or label, so that
   the property is averaged only over the atoms that match it. If
   however no atom labels match the argument given, then the average
   will be undefined. Note that for properties which are cumulatively
   counted rather than averaged, this error is not raised, and if no
   atom matches the label given 0 will be returned.

Socket errors
-------------

-  *Address already in use*: This is called if the server socket is
   already being used by the host network. There are several possible
   reasons for getting this error. Firstly, it might simply be that two
   simulations are running concurrently using the same host and port
   number. In this case simply change the port number of one of the
   simulations. Secondly, you can get this error if you try to rerun a
   simulation that previously threw an exception, since it takes a
   minute or so before the host will disconnect the server socket if it
   is not shut down cleanly. In this case, simply wait for it to
   disconnect, and try again. Finally, you will get this error if you
   try to use a restricted port number (i.e. below 1024) while not root.
   You should always use a non-restricted port number for i-PI
   simulations.

-  *Error opening unix socket. Check if a file /tmp/ipi__\_ exists, and
   remove it if unused.*: Similar to the above error, but given if you
   are using a unix socket rather than an internet socket. Since this
   binds locally the socket can be removed by the user, which means that
   it is not necessary to wait for the computer to automatically
   disconnect an unused server socket.

-  *Port number \_\_ out of acceptable range*: The port number must be
   between 1 and 65535, and should be greater than 1024. Change the port
   number accordingly.

-  *Slot number \_\_ out of acceptable range*: The slot number must be
   between 1 and 5. Change the slot number accordingly.

-  *’NoneType’ object has no attribute ’Up’*: This is called if an
   exception is raised during writing the data to output, and so the
   thread that deals with the socket is terminated uncleanly. Check the
   stack trace for the original exception, since this will be the actual
   source of the problem. Also note that, since the socket thread was
   not cleaned up correctly, the server socket may not have been
   disconnected properly and you may have to wait for a minute before
   you can restart a simulation using the same host and port number.

Mathematical errors
-------------------

-  *math domain error*: If the cell parameters are defined using the
   side lengths and angles, with either a pdb file or using the “abcABC”
   initialization mode, then for some value of the angles it is
   impossible to construct a valid cell vector matrix. This will cause
   the code to attempt to take the square root of a negative number,
   which gives this exception.

-  *overflow encountered in exp*: Sometimes occurs in *NPT* runs when
   the simulation box “explodes”. Make sure you have properly
   equilibrated the system before starting and that the timestep is
   short enough to not introduce very large integration errors.
