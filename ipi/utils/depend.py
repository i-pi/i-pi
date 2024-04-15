"""Tools that deal with the dependency detection, value caching and automatic
updating of variables.

The classes defined in this module overload the standard __get__ and __set__
routines of the numpy ndarray class and standard library object class so that
they automatically keep track of whether anything they depend on has been
altered, and so only recalculate their value when necessary.

Basic quantities that depend on nothing else can be manually altered in the
usual way, all other quantities are updated automatically and cannot be changed
directly.

The exceptions to this are synchronized properties, which are in effect
multiple basic quantities all related to each other, for example the bead and
normal mode representations of the positions and momenta. In this case any of
the representations can be set manually, and all the other representations
must keep in step.

For a more detailed discussion, see the reference manual.
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2017 i-PI developers
# See the "licenses" directory for full license information.


import weakref
import threading

from copy import deepcopy

import numpy as np

from ipi.utils.messages import verbosity, warning


__all__ = [
    "depend_value",
    "depend_array",
    "synchronizer",
    "dpipe",
    "dcopy",
    "dstrip",
    "depraise",
    "dproperties",
    "ddot",
    "noddot",
]


class synchronizer(object):
    """Class to implement synched objects.

    Holds the objects used to keep two or more objects in step with each other.
    This is shared between all the synched objects.

    Attributes:
        synched: A dictionary containing all the synched objects, of the form
            {"name": depend object}.
        manual: A string containing the name of the object being manually
            changed.
    """

    def __init__(self, deps=None):
        """Initialises synchronizer.

        Args:
           deps: Optional dictionary giving the synched objects of the form
                 {"name": depend object}.
        """

        if deps is None:
            self.synced = dict()
        else:
            self.synced = deps

        self.manual = None


# TODO put some error checks in the init to make sure that the object is initialized from consistent synchro and func states
class depend_base(object):
    """Base class for dependency handling.

    Builds the majority of the machinery required for the different depend
    objects. Contains functions to add and remove dependencies, the tainting
    mechanism by which information about which objects have been updated is
    passed around the dependency network, and the manual and automatic update
    functions to check that depend objects with functions are not manually
    updated and that synchronized objects are kept in step with the one
    manually changed.

    Attributes:
        _tainted: An array containing one boolean, which is True if one of the
            dependencies has been changed since the last time the value was
            cached.
        _func: A function name giving the method of calculating the value,
            if required. None otherwise.
        _name: The name of the depend base object.
        _synchro: A synchronizer object to deal with synched objects, if
            required. None otherwise.
        _dependants: A list containing all objects dependent on the self.
    """

    def __init__(
        self,
        name,
        synchro=None,
        func=None,
        dependants=None,
        dependencies=None,
        tainted=None,
    ):
        """Initialises depend_base.

        An unusual initialisation routine, as it has to be able to deal with
        the depend array mechanism for returning slices as new depend arrays.

        This is the reason for the penultimate if statement; it automatically
        taints objects created from scratch but does nothing to slices which
        are not tainted.

        Also, the last if statement makes sure that if a synchronized property
        is sliced, this initialization routine does not automatically set it to
        the manually updated property.

        Args:
            name: A string giving the name of self.
            tainted: An optional array containing one boolean which is True if
                one of the dependencies has been changed.
            func: An optional argument that can be specified either by a
                function name, or for synchronized values a dictionary of the
                form {"name": function name}; where "name" is one of the other
                synched objects and function name is the name of a function to
                get the object "name" from self.
            synchro: An optional synchronizer object.
            dependants: An optional list containing objects that depend on self.
            dependencies: An optional list containing objects that self
                depends upon.
        """

        if tainted is None:
            tainted = np.array([True], bool)
        if dependants is None:
            dependants = []
        if dependencies is None:
            dependencies = []

        self._tainted = tainted
        self._func = func
        self._name = name
        self._threadlock = threading.RLock()
        self._dependants = []
        self._synchro = None

        self.add_synchro(synchro)

        for item in dependencies:
            item.add_dependant(self, tainted)

        # Convert dependants to weakreferences consistently
        for item in dependants:
            if not isinstance(item, weakref.ref):
                dependants.remove(item)
                dependants.append(weakref.ref(item))
        self._dependants = dependants

        # Don't taint self if the object is a primitive one.
        # However, do propagate tainting to dependants if required.
        if tainted[0]:
            if self._func is None:
                self.taint(taintme=False)
            else:
                self.taint(taintme=tainted)

    def __deepcopy__(self, memo):
        """Overrides deepcopy behavior, to avoid copying the (uncopiable) RLOCK"""

        newone = type(self)(None)

        for member in newone.__dict__:
            if member == "_threadlock":
                continue
            setattr(newone, member, deepcopy(getattr(self, member), memo))

        return newone

    def add_synchro(self, synchro=None):
        """Links depend object to a synchronizer."""

        assert (
            self._synchro is None
        ), "This object must not have a previous synchronizer!"

        self._synchro = synchro
        if self._synchro is not None and self._name not in self._synchro.synced:
            self._synchro.synced[self._name] = self
            self._synchro.manual = self._name

    def add_dependant(self, newdep, tainted=True):
        """Adds a dependant property.

        Args:
            newdep: The depend object to be added to the dependency list.
            tainted: A boolean that decides whether newdep should be tainted.
               True by default.
        """

        newdep.add_dependency(self, tainted=tainted)

    def add_dependency(self, newdep, tainted=True):
        """Adds a dependency.

        Args:
            newdep: The depend object self now depends upon.
            tainted: A boolean that decides whether self should
                be tainted. True by default.
        """

        newdep._dependants.append(weakref.ref(self))
        if tainted:
            self.taint(taintme=True)

    def taint(self, taintme=True):
        """Recursively sets tainted flag on dependent objects.

        The main function dealing with the dependencies. Taints all objects
        further down the dependency tree until either all objects have been
        tainted, or it reaches only objects that have already been tainted.
        Note that in the case of a dependency loop the initial setting of
        _tainted to True prevents an infinite loop occuring.

        Also, in the case of a synchro object, the manually set quantity is not
        tainted, as it is assumed that synchro objects only depend on each
        other.

        Args:
           taintme: A boolean giving whether self should be tainted at the end.
              True by default.
        """

        # propagates dependency
        for item in self._dependants:
            item = item()
            if not item._tainted[0]:
                item.taint()

        if self._synchro is None:
            self._tainted[0] = taintme
        else:
            self._tainted[0] = False
            for v in self._synchro.synced.values():
                if not v._tainted[0]:
                    v._tainted[0] = v._name != self._synchro.manual
                    # do the propagation on the dependants here
                    # so we don't iterate multiple times over synchro
                    for item in v._dependants:
                        item = item()
                        if not item._tainted[0]:
                            item.taint()
            self._tainted[0] = taintme and (not self._name == self._synchro.manual)

    def tainted(self):
        """Returns tainted flag."""

        return self._tainted[0]

    def update_auto(self):
        """Automatic update routine.

        Updates the value when get has been called and self has been tainted.
        """

        if self._synchro is not None:
            if not self._name == self._synchro.manual:
                self.set(self._func[self._synchro.manual](), manual=False)
            else:
                warning(
                    self._name + " probably shouldn't be tainted (synchro)",
                    verbosity.low,
                )
        elif self._func is not None:
            self.set(self._func(), manual=False)
        else:
            warning(
                self._name + " probably shouldn't be tainted (value)", verbosity.low
            )

    def update_man(self):
        """Manual update routine.

        Updates the value when the value has been manually set. Also raises an
        exception if a calculated quantity has been manually set. Also starts
        the tainting routine.

        Raises:
            NameError: If a calculated quantity has been manually set.
        """

        if self._synchro is not None:
            self._synchro.manual = self._name
        elif self._func is not None:
            raise NameError(
                "Cannot set manually the value of the automatically-computed property <"
                + self._name
                + ">"
            )
        self.taint(taintme=False)

    def set(self, value, manual=False):
        """Dummy setting routine."""

        raise ValueError("Undefined set function for base depend class")

    def get(self):
        """Dummy getting routine."""

        raise ValueError("Undefined get function for base depend class")


class depend_value(depend_base):
    """Depend class for scalar values.

    Attributes:
        _value: The value associated with self.
    """

    def __init__(
        self,
        name,
        value=None,
        synchro=None,
        func=None,
        dependants=None,
        dependencies=None,
        tainted=None,
    ):
        """Initialises depend_value.

        Args:
            name: A string giving the name of self.
            value: The value of the object. Optional.
            tainted: An optional array giving the tainted flag. Default is [True].
            func: An optional argument that can be specified either by a function
                name, or for synchronized values a dictionary of the form
                {"name": function name}; where "name" is one of the other
                synched objects and function name is the name of a function to
                get the object "name" from self.
            synchro: An optional synchronizer object.
            dependants: An optional list containing objects that depend on self.
            dependencies: An optional list containing objects that self
                depends upon.
        """

        self._value = value
        super(depend_value, self).__init__(
            name, synchro, func, dependants, dependencies, tainted
        )

    def get(self):
        """Returns value, after recalculating if necessary."""

        return self.__get__(self, self.__class__)

    def __get__(self, instance=None, owner=None):
        """Overwrites standard get function

        Returns a cached value, recomputing only if the value is tainted.
        """

        if self._tainted[0]:
            with self._threadlock:
                self.update_auto()
                self.taint(taintme=False)

        return self._value

    def set(self, value, manual=True):
        """Alters value and taints dependencies.

        Overwrites the standard method of setting value, so that dependent
        quantities are tainted, and so we check that computed quantities are not
        manually updated.
        """

        with self._threadlock:
            self._value = value
            if manual:
                self.update_man()

    def __set__(self, instance, value):
        """Overwrites standard set function."""

        self.set(value)


class depend_array(np.ndarray, depend_base):
    """Depend class for arrays.

    Differs from depend_value as arrays handle getting items in a different
    way to scalar quantities, and as there needs to be support for slicing an
    array. Initialisation is also done in a different way for ndarrays.

    Attributes:
        _bval: The base deparray storage space. Equal to dstrip(self) unless
            self is a slice.
    """

    def __new__(
        cls,
        value,
        name,
        synchro=None,
        func=None,
        dependants=None,
        dependencies=None,
        tainted=None,
        base=None,
    ):
        """Creates a new array from a template.

        Called whenever a new instance of depend_array is created. Casts the
        array base into an appropriate form before passing it to
        __array_finalize__().

        Args:
           See __init__().
        """

        obj = np.asarray(value).view(cls)
        return obj

    def __init__(
        self,
        value,
        name,
        synchro=None,
        func=None,
        dependants=None,
        dependencies=None,
        tainted=None,
        base=None,
    ):
        """Initialises depend_array.

        Note that this is only called when a new array is created by an
        explicit constructor.

        Args:
            name: A string giving the name of self.
            value: The (numpy) array to serve as the memory base.
            tainted: An optional array giving the tainted flag. Default is [True].
            func: An optional argument that can be specified either by a function
                name, or for synchronized values a dictionary of the form
                {"name": function name}; where "name" is one of the other
                synched objects and function name is the name of a function to
                get the object "name" from self.
            synchro: An optional synchronizer object.
            dependants: An optional list containing objects that depend on self.
            dependencies: An optional list containing objects that self
                depends upon.
        """

        super(depend_array, self).__init__(
            name, synchro, func, dependants, dependencies, tainted
        )

        if base is None:
            self._bval = value
        else:
            self._bval = base

    def copy(self, order="C", maskna=None):
        """Wrapper for numpy copy mechanism."""

        # Sets a flag and hands control to the numpy copy
        self._fcopy = True
        return super(depend_array, self).copy(order)

    def __array_finalize__(self, obj):
        """Deals with properly creating some arrays.

        In the case where a function acting on a depend array returns a ndarray,
        this casts it into the correct form and gives it the
        depend machinery for other methods to be able to act upon it. New
        depend_arrays will next be passed to __init__ ()to be properly
        initialized, but some ways of creating arrays do not call __new__() or
        __init__(), so need to be initialized.
        """

        depend_base.__init__(self, name="")

        if type(obj) is depend_array:
            # We are in a view cast or in new from template. Unfortunately
            # there is no sure way to tell (or so it seems). Hence we need to
            # handle special cases, and hope we are in a view cast otherwise.
            if hasattr(obj, "_fcopy"):
                del obj._fcopy  # removes the "copy flag"
                self._bval = dstrip(self)
            else:
                # Assumes we are in view cast, so copy over the attributes from the
                # parent object. Typical case: when transpose is performed as a
                # view.
                super(depend_array, self).__init__(
                    obj._name,
                    obj._synchro,
                    obj._func,
                    obj._dependants,
                    None,
                    obj._tainted,
                )
                self._bval = obj._bval
        else:
            # Most likely we came here on the way to init.
            # Just sets a defaults for safety
            self._bval = dstrip(self)

    def __array_prepare__(self, arr, context=None):
        """Prepare output array for ufunc.

        Depending on the context we try to understand if we are doing an
        in-place operation (in which case we want to keep the return value a
        deparray) or we are generating a new array as a result of the ufunc.
        In this case there is no way to know if dependencies should be copied,
        so we strip and return a ndarray.
        """

        if context is None or len(context) < 2 or not type(context[0]) is np.ufunc:
            # It is not clear what we should do. If in doubt, strip dependencies.
            return np.ndarray.__array_prepare__(
                self.view(np.ndarray), arr.view(np.ndarray), context
            )
        elif len(context[1]) > context[0].nin and context[0].nout > 0:
            # We are being called by a ufunc with a output argument, which is being
            # actually used. Most likely, something like an increment,
            # so we pass on a deparray
            return super(depend_array, self).__array_prepare__(arr, context)
        else:
            # Apparently we are generating a new array.
            # We have no way of knowing its
            # dependencies, so we'd better return a ndarray view!
            return np.ndarray.__array_prepare__(
                self.view(np.ndarray), arr.view(np.ndarray), context
            )

    def __array_wrap__(self, arr, context=None):
        """Wraps up output array from ufunc.

        See docstring of __array_prepare__().
        """

        if context is None or len(context) < 2 or not type(context[0]) is np.ufunc:
            return np.ndarray.__array_wrap__(
                self.view(np.ndarray), arr.view(np.ndarray), context
            )
        elif len(context[1]) > context[0].nin and context[0].nout > 0:
            return super(depend_array, self).__array_wrap__(arr, context)
        else:
            return np.ndarray.__array_wrap__(
                self.view(np.ndarray), arr.view(np.ndarray), context
            )

    # whenever possible in compound operations just return a regular ndarray
    __array_priority__ = -1.0

    def reshape(self, newshape):
        """Changes the shape of the base array.

        Args:
           newshape: A tuple giving the desired shape of the new array.

        Returns:
           A depend_array with the dimensions given by newshape.
        """

        return depend_array(
            dstrip(self).reshape(newshape),
            name=self._name,
            synchro=self._synchro,
            func=self._func,
            dependants=self._dependants,
            tainted=self._tainted,
            base=self._bval,
        )

    def flatten(self):
        """Makes the base array one dimensional.

        Returns:
            A flattened array.
        """

        return self.reshape(self.size)

    @staticmethod
    def __scalarindex(index, depth=1):
        """Checks if an index points at a scalar value.

        Used so that looking up one item in an array returns a scalar, whereas
        looking up a slice of the array returns a new array with the same
        dependencies as the original, so that changing the slice also taints
        the global array.

        Arguments:
            index: the index to be checked.
            depth: the rank of the array which is being accessed. Default value
                is 1.

        Returns:
            A logical stating whether a __get__ instruction based
            on index would return a scalar.
        """

        if hasattr(index, "__len__"):
            if len(index) == depth and isinstance(index, tuple):
                # if the index is a tuple check it does not contain slices
                for i in index:
                    if isinstance(i, slice):
                        return False
                return True
        elif depth <= 1:
            return True

        return False

    def __getitem__(self, index):
        """Returns value[index], after recalculating if necessary.

        Overwrites the standard method of getting value, so that value
        is recalculated if tainted. Scalar slices are returned as an ndarray,
        so without depend machinery. If you need a "scalar depend" which
        behaves as a slice, just create a 1x1 matrix, e.g b=a(7,1:2)

        Args:
           index: A slice variable giving the appropriate slice to be read.
        """

        if self._tainted[0]:
            with self._threadlock:
                self.update_auto()
                self.taint(taintme=False)

        if self.__scalarindex(index, self.ndim):
            return self.view(np.ndarray)[index]
        else:
            return super(depend_array, self).__getitem__(index)

    def __getslice__(self, i, j):
        """Overwrites standard get function."""

        return self.__getitem__(slice(i, j, None))

    def get(self):
        """Alternative to standard get function."""

        return self.__getitem__(slice(None, None, None))

    def __get__(self, instance=None, owner=None):
        """Overwrites standard get function."""

        # It is worth duplicating this code that is also used in __getitem__ as this
        # is called most of the time, and we avoid creating a load of copies pointing to the same depend_array

        if self._tainted[0]:
            with self._threadlock:
                self.update_auto()
                self.taint(taintme=False)

        return self

    def __setitem__(self, index, value, manual=True):
        """Alters value[index] and taints dependencies.

        Overwrites the standard method of setting value, so that dependent
        quantities are tainted, and so we check that computed quantities are not
        manually updated.

        Args:
           index: A slice variable giving the appropriate slice to be read.
           value: The new value of the slice.
           manual: Optional boolean giving whether the value has been changed
              manually. True by default.
        """

        with self._threadlock:
            if manual:
                self.view(np.ndarray)[index] = value
                self.update_man()
            elif index == slice(None, None, None):
                self._bval[index] = value
                self.taint(taintme=False)
            else:
                raise IndexError(
                    "Automatically computed arrays should span the whole parent"
                )

    def __setslice__(self, i, j, value):
        """Overwrites standard set function."""

        return self.__setitem__(slice(i, j), value)

    def set(self, value, manual=True):
        """Alterative to standard set function.

        Args:
            See __setitem__().
        """

        self.__setitem__(slice(None, None), value=value, manual=manual)

    def __set__(self, instance, value):
        """Overwrites standard set function."""

        self.__setitem__(slice(None, None), value=value)


# BEGINS NUMPY FUNCTIONS OVERRIDE
# np.dot and other numpy.linalg functions have the nasty habit to
# view cast to generate the output. Since we don't want to pass on
# dependencies to the result of these functions, and we can't use
# the ufunc mechanism to demote the class type to ndarray, we must
# overwrite np.dot and other similar functions.

# ** np.dot
noddot = np.dot


def ddot(da, db):
    a = dstrip(da)
    b = dstrip(db)

    return noddot(a, b)


# activate this to have safe dstripping
# np.dot = ddot
# activate this if you want to assume dot will almost always be applied on dstripped vectors
np.dot = noddot

# ENDS NUMPY FUNCTIONS OVERRIDE


def dstrip(da):
    """Removes dependencies from a depend_array.

    Takes a depend_array and returns its value as a ndarray, effectively
    stripping the dependencies from the ndarray. This speeds up a lot of
    calculations involving these arrays. Must only be used if the value of the
    array is not going to be changed.

    Args:
        deparray: A depend_array.

    Returns:
        A ndarray with the same value as deparray.
    """

    try:
        return da.view(np.ndarray)
    except:  # TODO: remove all remaining stray dstrip so we don't need to check here
        warning("dstrip should only be called on `depend_array`s", verbosity.low)
        return da


def dpipe(dfrom, dto, item=-1):
    """Synchonizes two depend objects.

    Takes two depend objects, and makes one of them depend on the other in such
    a way that both keep the same value. Used for attributes such as
    temperature that are used in many different modules, and so need different
    depend objects in each, but which should all have the same value.

    Args:
        dfrom: The object that is depend on.
        dto: The object that depends on the former one.
    """

    if item < 0:
        dto._func = lambda: dfrom.__get__(dfrom, dfrom.__class__)
    else:
        dto._func = lambda i=item: dfrom.__getitem__(i)
    dto.add_dependency(dfrom)


def dcopy(dfrom, dto):
    """Copies the dependencies of one depend object to another.

    Args:
        see dpipe.
    """
    dto._dependants = dfrom._dependants
    dto._synchro = dfrom._synchro
    dto.add_synchro(dfrom._synchro)
    dto._tainted = dfrom._tainted
    dto._func = dfrom._func
    if hasattr(dfrom, "_bval"):
        dto._bval = dfrom._bval


def depraise(exception):
    raise exception


def _inject_depend_property(cls, attr_name):
    private_name = f"_{attr_name}"

    def getter(self):
        return getattr(self, private_name).__get__(self, cls)

    def setter(self, value):
        return getattr(self, private_name).__set__(self, value)

    setattr(cls, attr_name, property(getter, setter))


def dproperties(cls, attr_names):
    """
    Adds to a class property wrappers to access its depend members.
    Assumes that depend objects are named with an underscore prefix, and
    creates getters and setters with the given names.

    If a class `A` contains a `_val` depend object, one should use
    `dproperties(A, ["val"])`.
    """
    if not isinstance(attr_names, list):
        attr_names = [attr_names]

    for attr_name in attr_names:
        _inject_depend_property(cls, attr_name)
