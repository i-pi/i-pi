"""Holds the algorithms to perform replica exchange.

Algorithms implemented by Robert Meissner and Riccardo Petraglia, 2016
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2016 i-PI developers
# See the "licenses" directory for full license information.


from ipi.engine.motion import Motion


__all__ = ["TemperatureRamp", "PressureRamp"]


class TemperatureRamp(Motion):
    """Temperature ramp (quench/heat).

    Attributes:

    """

    def __init__(
        self,
        fixcom=False,
        fixatoms_dof=None,
        t_start=1.0,
        t_end=1.0,
        total_steps=0,
        current_step=0,
        logscale=True,
    ):
        """Initialises a temperature ramp motion

        Args:
        """

        super(TemperatureRamp, self).__init__()
        self.t_start = t_start
        self.t_end = t_end
        self.total_steps = total_steps
        self.current_step = current_step
        self.logscale = logscale

    def bind(self, ens, beads, nm, cell, bforce, prng, omaker):
        super(TemperatureRamp, self).bind(ens, beads, nm, cell, bforce, prng, omaker)

    def step(self, step=None):
        """Updates ensemble temperature. Yes, that's it everything else should follow through"""

        self.current_step += 1
        # just sets the temperature
        if self.current_step >= self.total_steps:
            self.ensemble.temp = self.t_end
        else:
            if self.logscale:
                self.ensemble.temp = self.t_start * (self.t_end / self.t_start) ** (
                    self.current_step * 1.0 / self.total_steps
                )
            else:
                self.ensemble.temp = (
                    self.t_start
                    + self.current_step * (self.t_end - self.t_start) / self.total_steps
                )


class PressureRamp(Motion):
    """Pressure ramp (quench/heat).

    Attributes:

    """

    def __init__(
        self,
        fixcom=False,
        fixatoms_dof=None,
        p_start=1.0,
        p_end=1.0,
        total_steps=0,
        current_step=0,
        logscale=True,
    ):
        """Initialises a temperature ramp motion

        Args:
        """

        super(PressureRamp, self).__init__()
        self.p_start = p_start
        self.p_end = p_end
        self.total_steps = total_steps
        self.current_step = current_step
        self.logscale = logscale

    def bind(self, ens, beads, nm, cell, bforce, prng, omaker):
        super(PressureRamp, self).bind(ens, beads, nm, cell, bforce, prng, omaker)

    def step(self, step=None):
        """Updates ensemble temperature. Yes, that's it everything else should follow through"""

        self.current_step += 1
        # just sets the temperature
        if self.current_step >= self.total_steps:
            self.ensemble.pext = self.p_end
        else:
            if self.logscale:
                self.ensemble.pext = self.p_start * (self.p_end / self.p_start) ** (
                    self.current_step * 1.0 / self.total_steps
                )
            else:
                self.ensemble.pext = (
                    self.p_start
                    + self.current_step * (self.p_end - self.p_start) / self.total_steps
                )
