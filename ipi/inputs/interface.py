"""Creates objects that create the interface to drivers."""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


from ipi.utils.messages import verbosity, warning
from ipi.utils.inputvalue import *
from ipi.interfaces.sockets import *


__all__ = ["InputInterfaceSocket"]


class InputInterfaceSocket(Input):
    """Interface input class.

    Handles generating the apporopriate interface class from the xml
    input file, and generating the xml checkpoint tags and data from an
    instance of the object.

    Attributes:
       mode: A string giving the type of socket used.
       pbc: A boolean giving whether the atom positions will be folded back
          into the unit cell before being sent through the socket or not.
          Default is false.

    Fields:
       address: A string giving the host name.
       port: An integer giving the port used by the socket.
       slots: An integer giving the maximum allowed backlog of queued clients.
       latency: A float giving the number of seconds that the interface waits
          before updating the client list.
       timeout: A float giving a number of seconds after which a calculation core
          is considered dead. Defaults to one hour.
    """

    fields = {
        "address": (
            InputValue,
            {
                "dtype": str,
                "default": "localhost",
                "help": "This gives the server address that the socket will run on.",
            },
        ),
        "port": (
            InputValue,
            {
                "dtype": int,
                "default": 65535,
                "help": "This gives the port number that defines the socket.",
            },
        ),
        "slots": (
            InputValue,
            {
                "dtype": int,
                "default": 4,
                "help": "This gives the number of client codes that can queue at any one time.",
            },
        ),
        "latency": (
            InputValue,
            {
                "dtype": float,
                "default": 1e-3,
                "help": "This gives the number of seconds between each check for new clients.",
            },
        ),
        "max_workers": (
            InputValue,
            {
                "dtype": int,
                "default": 128,
                "help": "This gives the maximum number of job threads that are active simultaneously.",
            },
        ),
        "timeout": (
            InputValue,
            {
                "dtype": float,
                "default": 3600.0,
                "help": "This gives the number of seconds before assuming a calculation has died. If 0 there is no timeout.",
            },
        ),
    }
    attribs = {
        "mode": (
            InputAttribute,
            {
                "dtype": str,
                "options": ["unix", "inet"],
                "default": "inet",
                "help": "Specifies whether the driver interface will listen onto a internet socket [inet] or onto a unix socket [unix].",
            },
        ),
        "pbc": (
            InputAttribute,
            {
                "dtype": bool,
                "default": False,
                "help": "Applies periodic boundary conditions to the atoms coordinates before passing them on to the driver code.",
            },
        ),
    }

    default_help = "Specifies the parameters for the socket interface."
    default_label = "INTERFACE"

    def store(self, iface):
        """Takes an Interface instance and stores a minimal representation of it.

        Args:
           iface: An interface object.
        """

        super(InputInterfaceSocket, self).store(iface)
        self.latency.store(iface.latency)
        self.mode.store(iface.mode)
        self.address.store(iface.address)
        self.port.store(iface.port)
        self.slots.store(iface.slots)
        self.timeout.store(iface.timeout)
        self.pbc.store(iface.dopbc)
        self.max_workers.store(iface.max_workers)

    def fetch(self):
        """Creates an InterfaceSocket object.

        Returns:
           An interface object with the appropriate socket given the attributes
           of the InputInterfaceSocket object.
        """

        super(InputInterfaceSocket, self).fetch()
        return InterfaceSocket(
            address=self.address.fetch(),
            port=self.port.fetch(),
            slots=self.slots.fetch(),
            mode=self.mode.fetch(),
            latency=self.latency.fetch(),
            timeout=self.timeout.fetch(),
            dopbc=self.pbc.fetch(),
            max_workers=self.max_workers.fetch(),
        )

    def check(self):
        """Function that deals with optional arguments."""

        super(InputInterfaceSocket, self).check()
        if self.port.fetch() < 1 or self.port.fetch() > 65535:
            raise ValueError(
                "Port number " + str(self.port.fetch()) + " out of acceptable range."
            )
        elif self.port.fetch() < 1025:
            warning(
                "Low port number being used, this may interrupt important system processes.",
                verbosity.low,
            )

        if self.slots.fetch() < 1 or self.slots.fetch() > 5:
            raise ValueError(
                "Slot number " + str(self.slots.fetch()) + " out of acceptable range."
            )
        if self.latency.fetch() < 0:
            raise ValueError("Negative latency parameter specified.")
        if self.timeout.fetch() < 0.0:
            raise ValueError("Negative timeout parameter specified.")
