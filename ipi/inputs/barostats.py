"""Creates objects that deal with constant pressure and stress simulations."""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


import numpy as np

import ipi.engine.thermostats
from ipi.engine.barostats import *
from ipi.utils.inputvalue import *
from ipi.utils.array_backend import xp
from ipi.inputs.thermostats import *
from ipi.engine.cell import Cell
from ipi.inputs.cell import *

__all__ = ["InputBaro"]


class InputBaro(Input):
    """Barostat input class.

    Handles generating the appropriate barostat class from the xml input file,
    and generating the xml checkpoint tags and data from an
    instance of the object.

    Attributes:
       mode: An optional string giving the type of barostat used. Defaults to
          'dummy'.

    Fields:
       thermostat: A thermostat object giving the cell thermostat.
       tau: The time constant associated with the dynamics of the cell.
       bulk_modulus: An estimate of the bulk modulus, which sets the mass of
          the cell, or the size of the moves of stochastic cell rescaling.
       v0: The reference volume for the mass of the cell.
       ebaro: The accumulated stochastic cell-rescaling effective-energy
          correction.
       p: The momentum of the cell.
       h0: The reference cell of the anisotropic barostat.
       hfix: The elements of the cell that are kept fixed.
       vol_constraint: Whether the volume of the cell is kept fixed.
    """

    attribs = {
        "mode": (
            InputAttribute,
            {
                "dtype": str,
                "default": "dummy",
                "help": """The type of barostat. 'isotropic' implements the Bussi-Zykova-Parrinello barostat [doi:10.1063/1.3073889], which scales the cell isotropically while sampling the isothermal-isobaric ensemble, and is suitable for liquids; the implementation details are given in [doi:10.1016/j.cpc.2013.10.027]. 'stochastic-rescaling' implements the Bernetti-Bussi stochastic cell-rescaling barostat [doi:10.1063/5.0020514] for classical dynamics, with the reversible Trotter scheme and the instantaneous kinetic pressure; it requires an estimate of the bulk modulus and no cell thermostat. 'sc-isotropic' implements the Bussi-Zykova-Parrinello barostat for Suzuki-Chin path integral molecular dynamics [10.1021/acs.jctc.8b01297], for the 'scnpt' dynamics. 'flexible' implements the path integral version of the Martyna-Tuckerman-Tobias-Klein barostat [doi:10.1063/1.478193], which changes the lengths and angles of the cell at constant pressure, for anisotropic systems such as molecular solids. 'anisotropic' implements the Raiteri-Gale-Bussi barostat [10.1088/0953-8984/23/33/334213], which changes the lengths and angles of the cell at a constant (possibly non-diagonal) stress, with the 'nst' dynamics; it requires a reference cell h0, and is valid within the elastic limit of small strains. For diagonal stresses (or external pressures) and small strains, 'flexible' and 'anisotropic' give very similar results. 'dummy' does nothing.""",
                "options": [
                    "dummy",
                    "isotropic",
                    "stochastic-rescaling",
                    "flexible",
                    "anisotropic",
                    "sc-isotropic",
                ],
            },
        )
    }
    fields = {
        "thermostat": (
            InputThermo,
            {
                "default": input_default(factory=ipi.engine.thermostats.Thermostat),
                "help": "The thermostat for the cell. Keeps the cell velocity distribution at the correct temperature. Note that the 'pile_l', 'pile_g', 'nm_gle' and 'nm_gle_g' options will not work for this thermostat, and that 'stochastic-rescaling' uses none.",
            },
        ),
        "tau": (
            InputValue,
            {
                "default": 1.0,
                "dtype": float,
                "dimension": "time",
                "help": "The time scale of the barostat. For 'stochastic-rescaling' it is the relaxation time of the volume. For the other modes it sets the mass of the cell (see bulk_modulus): if bulk_modulus is given, the volume oscillates with a period of about tau, for any material and temperature, and otherwise with a period of 2 pi tau / sqrt(3 B / (rho kB T)) for a material of bulk modulus B and atomic density rho (about 1.6 tau for water and 0.2 tau for copper at 300 K). With 'flexible' and 'anisotropic', the shear deformations of the cell oscillate with frequencies set by the shear moduli rather than by B.",
            },
        ),
        "bulk_modulus": (
            InputValue,
            {
                "default": -1.0,
                "dtype": float,
                "dimension": "pressure",
                "help": "An estimate of the isothermal bulk modulus B of the system, which only affects the dynamics, not the ensemble that is sampled. It must be given for 'stochastic-rescaling', where it sets the size of the volume moves. For the other modes it sets the mass of the cell to 9 V B (tau / 2 pi)^2 for 'isotropic' and 'sc-isotropic', and to 3 V B (tau / 2 pi)^2 for each element of the cell momentum for 'flexible' and 'anisotropic', with V the reference volume v0, or the volume of h0 for 'anisotropic', so that tau is the period of the volume oscillations. Without it the masses are 3 N kB T tau^2 and N kB T tau^2, with N the number of atoms. Both masses are multiplied by the number of beads.",
            },
        ),
        "v0": (
            InputValue,
            {
                "default": -1.0,
                "dtype": float,
                "dimension": "volume",
                "help": "The reference volume V in the mass of the cell (see bulk_modulus), for the 'isotropic', 'sc-isotropic' and 'flexible' modes with a bulk_modulus. Defaults to the volume of the cell at the start, and is written in restarts so that they continue with the same mass. A value close to the average volume makes the dynamics independent of the initial cell.",
            },
        ),
        "ebaro": (
            InputValue,
            {
                "default": 0.0,
                "dtype": float,
                "dimension": "energy",
                "help": "Accumulated stochastic cell-rescaling effective-energy correction. This is written to checkpoints and should normally not be set manually.",
            },
        ),
        "p": (
            InputArray,
            {
                "dtype": float,
                "default": input_default(factory=np.zeros, args=(0,)),
                "help": "The momentum of the cell: one value for 'isotropic' and 'sc-isotropic', and a 3x3 upper-triangular matrix for 'flexible' and 'anisotropic'.",
                "dimension": "momentum",
            },
        ),
        "h0": (
            InputCell,
            {
                "dtype": float,
                "default": input_default(factory=Cell),
                "help": "The reference cell of the 'anisotropic' barostat, which requires it. The external stress acts on the strain of the cell relative to h0, and equals the stress on the cell only when the cell is close to h0, which should therefore be close to the average cell.",
                "dimension": "length",
            },
        ),
        "hfix": (
            InputArray,
            {
                "default": input_default(factory=np.zeros, args=(0, str)),
                "dtype": str,
                "help": "A list of the cell entries that should be held fixed (xx, yy, zz, xy, xz, yz), whose velocity is set to zero. 'offdiagonal' is an alias for xy, xz, yz. xz can only be fixed together with xy or yz, as otherwise it changes when both of them do.",
            },
        ),
        "vol_constraint": (
            InputValue,
            {
                "default": False,
                "dtype": bool,
                "help": "If True, the (N, V, sigma_a = 0, T) ensemble [doi:10.1021/acs.jctc.5b00748] is sampled, in which the shape of the cell changes at constant volume. Only for the 'flexible' barostat, and not with fixed diagonal elements of the cell.",
            },
        ),
    }

    default_help = "Simulates an external pressure bath."
    default_label = "BAROSTAT"

    def store(self, baro):
        """Takes a barostat instance and stores a minimal representation of it.

        Args:
           baro: A barostat object.
        """

        super(InputBaro, self).store(baro)
        self.thermostat.store(baro.thermostat)
        self.tau.store(baro.tau)
        if type(baro) is BaroSCR:
            self.mode.store("stochastic-rescaling")
            self.bulk_modulus.store(baro.bulk_modulus)
            self.ebaro.store(baro.ebaro)
        elif type(baro) in (BaroBZP, BaroSCBZP):
            self.mode.store("isotropic" if type(baro) is BaroBZP else "sc-isotropic")
            self.p.store(baro.p)
            if baro.bulk_modulus > 0.0:
                self.bulk_modulus.store(baro.bulk_modulus)
                self.v0.store(baro.v0)
        elif type(baro) is BaroMTK:
            self.mode.store("flexible")
            self.p.store(baro.p)
            self.hfix.store(baro.hfix)
            self.vol_constraint.store(baro.vol_constraint)
            if baro.bulk_modulus > 0.0:
                self.bulk_modulus.store(baro.bulk_modulus)
                self.v0.store(baro.v0)
        elif type(baro) is BaroRGB:
            self.mode.store("anisotropic")
            self.p.store(baro.p)
            self.h0.store(baro.h0)
            self.hfix.store(baro.hfix)
            if baro.bulk_modulus > 0.0:
                self.bulk_modulus.store(baro.bulk_modulus)
        elif type(baro) is Barostat:
            self.mode.store("dummy")
        else:
            raise TypeError(
                "The type " + type(baro).__name__ + " is not a valid barostat type"
            )

    def fetch(self):
        """Creates a barostat object.

        Returns:
           A barostat object of the appropriate type and with the appropriate
           thermostat given the attributes of the InputBaro object.
        """

        super(InputBaro, self).fetch()
        if self.vol_constraint.fetch() and self.mode.fetch() != "flexible":
            raise ValueError(
                "Constraining the volume is not supported for the "
                + self.mode.fetch()
                + " barostat, use flexible instead."
            )
        if self.vol_constraint.fetch() and any(
            _ in self.hfix.fetch() for _ in ("xx", "yy", "zz")
        ):
            raise ValueError(
                "Diagonal cell entries cannot be fixed while constraining the volume."
            )
        if self.bulk_modulus._explicit and self.mode.fetch() == "dummy":
            raise ValueError("The bulk modulus is not used by the dummy barostat.")
        if self.v0._explicit and (
            self.mode.fetch() not in ("isotropic", "sc-isotropic", "flexible")
            or self.bulk_modulus.fetch() <= 0.0
        ):
            raise ValueError(
                "The reference volume v0 is only used by the 'isotropic', "
                "'sc-isotropic' and 'flexible' barostats with a bulk modulus."
            )
        if self.h0._explicit and self.mode.fetch() != "anisotropic":
            raise ValueError(
                "The reference cell h0 is only used by the 'anisotropic' barostat."
            )
        if self.mode.fetch() == "stochastic-rescaling":
            if len(self.hfix.fetch()) > 0:
                raise ValueError(
                    "Cannot fix individual cell components with a "
                    "'stochastic-rescaling' barostat"
                )
            if self.thermostat.mode.fetch() != "":
                raise ValueError(
                    "The 'stochastic-rescaling' barostat does not use a separate "
                    "cell thermostat. Remove the thermostat nested inside the "
                    "barostat tag."
                )
            if not self.bulk_modulus._explicit:
                raise ValueError(
                    "The bulk modulus must be specified explicitly for the "
                    "'stochastic-rescaling' barostat."
                )
            bulk_modulus = self.bulk_modulus.fetch()
            if bulk_modulus <= 0.0:
                raise ValueError(
                    "The stochastic cell-rescaling bulk modulus must be positive."
                )
            if self.tau.fetch() <= 0.0:
                raise ValueError(
                    "The stochastic cell-rescaling relaxation time must be positive."
                )
            baro = BaroSCR(
                thermostat=self.thermostat.fetch(),
                tau=self.tau.fetch(),
                bulk_modulus=bulk_modulus,
                ebaro=self.ebaro.fetch(),
            )
        elif self.mode.fetch() == "isotropic":
            if len(self.hfix.fetch()) > 0:
                raise ValueError(
                    "Cannot fix individual cell components with an 'isotropic' barostat"
                )
            baro = BaroBZP(
                thermostat=self.thermostat.fetch(),
                tau=self.tau.fetch(),
                bulk_modulus=self.bulk_modulus.fetch(),
                v0=self.v0.fetch(),
            )
            if self.p._explicit:
                baro.p = xp.asarray(self.p.fetch())
        elif self.mode.fetch() == "sc-isotropic":
            if len(self.hfix.fetch()) > 0:
                raise ValueError(
                    "Cannot fix individual cell components with an 'sc-isotropic' barostat"
                )
            baro = BaroSCBZP(
                thermostat=self.thermostat.fetch(),
                tau=self.tau.fetch(),
                bulk_modulus=self.bulk_modulus.fetch(),
                v0=self.v0.fetch(),
            )
            if self.p._explicit:
                baro.p = xp.asarray(self.p.fetch())
        elif self.mode.fetch() == "flexible":
            baro = BaroMTK(
                thermostat=self.thermostat.fetch(),
                tau=self.tau.fetch(),
                hfix=self.hfix.fetch(),
                vol_constraint=self.vol_constraint.fetch(),
                bulk_modulus=self.bulk_modulus.fetch(),
                v0=self.v0.fetch(),
            )
            if self.p._explicit:
                baro.p = xp.asarray(self.p.fetch())
        elif self.mode.fetch() == "anisotropic":
            baro = BaroRGB(
                thermostat=self.thermostat.fetch(),
                tau=self.tau.fetch(),
                hfix=self.hfix.fetch(),
                bulk_modulus=self.bulk_modulus.fetch(),
            )
            if self.p._explicit:
                baro.p = xp.asarray(self.p.fetch())
            if self.h0._explicit:
                baro.h0 = self.h0.fetch()
            else:
                raise ValueError(
                    "Reference cell MUST be specified for an anisotropic barostat"
                )
        elif self.mode.fetch() == "dummy":
            baro = Barostat(thermostat=self.thermostat.fetch(), tau=self.tau.fetch())
        else:
            raise ValueError(self.mode.fetch() + " is not a valid mode of barostat")

        return baro
