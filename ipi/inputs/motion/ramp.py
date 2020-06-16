"""Deals with creating the ensembles class.

Copyright (C) 2013, Joshua More and Michele Ceriotti

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http.//www.gnu.org/licenses/>.


Classes:
   InputEnsemble: Deals with creating the Ensemble object from a file, and
      writing the checkpoints.
"""

from ipi.engine.motion import *
from ipi.utils.inputvalue import *
from ipi.utils.units import *

__all__ = ["InputTemperatureRamp", "InputPressureRamp"]


class InputTemperatureRamp(InputDictionary):
    """Temperature ramp options.

    Contains options controlling a temperature ramp (quench/heating)

    """

    fields = {
        "t_start": (
            InputValue,
            {
                "dtype": float,
                "dimension": "energy",
                "default": 1.0,
                "help": "Initial temperature",
            },
        ),
        "t_end": (
            InputValue,
            {
                "dtype": float,
                "dimension": "energy",
                "default": 1.0,
                "help": "Final temperature",
            },
        ),
        "logscale": (
            InputValue,
            {
                "dtype": bool,
                "default": False,
                "help": "Change temperature on a logarihthmic scale.",
            },
        ),
        "total_steps": (
            InputValue,
            {"dtype": int, "default": 0, "help": "Total number of steps for the ramp"},
        ),
        "current_step": (
            InputValue,
            {"dtype": int, "default": 0, "help": "Current step along the ramp"},
        ),
    }

    default_help = """TemperatureRamp Motion class. It just updates the ensemble
                    temperature in steps, between the indicated temperatures, and
                    then holds to the highest value. It should typically be combined
                    with a dynamics class and thermostats, using a MultiMotion."""
    default_label = "TRAMP"

    def store(self, ramp):
        if ramp == {}:
            return
        self.t_start.store(ramp.t_start)
        self.t_end.store(ramp.t_end)
        self.logscale.store(ramp.logscale)
        self.total_steps.store(ramp.total_steps)
        self.current_step.store(ramp.current_step)

    def fetch(self):
        rv = super(InputTemperatureRamp, self).fetch()
        return rv


class InputPressureRamp(InputDictionary):
    """Pressure ramp options.

    Contains options controlling a pressure ramp

    """

    fields = {
        "p_start": (
            InputValue,
            {
                "dtype": float,
                "dimension": "pressure",
                "default": 1.0,
                "help": "Initial pressure",
            },
        ),
        "p_end": (
            InputValue,
            {
                "dtype": float,
                "dimension": "pressure",
                "default": 1.0,
                "help": "Final pressure",
            },
        ),
        "logscale": (
            InputValue,
            {
                "dtype": bool,
                "default": False,
                "help": "Change pressure on a logarihthmic scale.",
            },
        ),
        "total_steps": (
            InputValue,
            {"dtype": int, "default": 0, "help": "Total number of steps for the ramp"},
        ),
        "current_step": (
            InputValue,
            {"dtype": int, "default": 0, "help": "Current step along the ramp"},
        ),
    }

    default_help = """PressureRamp Motion class. It just updates the ensemble
                    pressure in steps, between the indicated values, and
                    then holds to the highest value. It should typically be combined
                    with a dynamics class and barostats, using a MultiMotion."""

    default_label = "PRAMP"

    def store(self, ramp):
        if ramp == {}:
            return
        self.p_start.store(ramp.p_start)
        self.p_end.store(ramp.p_end)
        self.logscale.store(ramp.logscale)
        self.total_steps.store(ramp.total_steps)
        self.current_step.store(ramp.current_step)

    def fetch(self):
        rv = super(InputPressureRamp, self).fetch()
        return rv
