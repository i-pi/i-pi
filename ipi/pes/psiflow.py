import warnings
import tempfile
import os

import numpy as np
from dataclasses import dataclass
from typing import ClassVar, Union, Optional
from pathlib import Path
import json
from typing import get_type_hints
import typeguard
from ase import Atoms
from ase.data import chemical_symbols, atomic_masses
from ase.io import read
from ase.units import fs, kJ, mol, nm

from .dummy import Dummy_driver
from ipi.utils.units import unit_to_internal, unit_to_user

__DRIVER_NAME__ = "psiflow"
__DRIVER_CLASS__ = "Psiflow_driver"

class Psiflow_driver(Dummy_driver):

    def __init__(self, template, hamiltonian, *args, **kwargs):
        self.template = template
        self.hamiltonian = hamiltonian
        super().__init__(*args, **kwargs)

    def check_parameters(self):
        self.template_geometry = Geometry.from_atoms(read(self.template))
        self.function = function_from_json(self.hamiltonian)

    def __call__(self, cell, pos):

        pos = unit_to_user("length", "angstrom", pos)
        cell = unit_to_user("length", "angstrom", cell.T)

        self.template_geometry.per_atom.positions[:] = pos
        if self.template_geometry.periodic:
            self.template_geometry.cell[:] = cell

        outputs = self.function(self.template_geometry)
        energy = outputs["energy"]
        forces = outputs["forces"]
        stress = outputs["stress"]

        # converts to internal quantities
        pot_ipi = np.asarray(
            unit_to_internal("energy", "electronvolt", energy), np.float64
        )
        force_ipi = np.asarray(unit_to_internal("force", "ev/ang", forces), np.float64)
        vir_calc = -stress * self.template_geometry.volume
        vir_ipi = np.array(
            unit_to_internal("energy", "electronvolt", vir_calc.T), dtype=np.float64
        )
        extras = ""

        return pot_ipi, force_ipi, vir_ipi, extras


@typeguard.typechecked
@dataclass
class Function:
    outputs: ClassVar[tuple]

    def __call__(
        self,
        geometry
    ) -> dict[str, float | np.ndarray]:
        raise NotImplementedError


@dataclass
class EnergyFunction(Function):
    outputs: ClassVar[tuple[str, ...]] = ("energy", "forces", "stress")

@dataclass
class EinsteinCrystalFunction(EnergyFunction):
    force_constant: float
    centers: np.ndarray
    volume: float = 0.0

    def __call__(
        self,
        geometry,
    ) -> dict[str, float | np.ndarray]:
        delta = geometry.per_atom.positions - self.centers
        energy = self.force_constant * np.sum(delta**2) / 2
        grad_pos = (-1.0) * self.force_constant * delta
        if geometry.periodic and self.volume > 0.0:
            delta = np.linalg.det(geometry.cell) - self.volume
            _stress = self.force_constant * np.eye(3) * delta
        else:
            _stress = np.zeros((3, 3))
        grad_cell = _stress
        return {"energy": energy, "forces": grad_pos, "stress": grad_cell}


@typeguard.typechecked
@dataclass
class PlumedFunction(EnergyFunction):
    plumed_input: str
    external: Optional[Union[str, Path]] = None

    def __post_init__(self):
        self.plumed_instances = {}

    def __call__(
        self,
        geometry,
    ) -> dict[str, float | np.ndarray]:
        plumed_input = self.plumed_input
        if self.external is not None:
            assert self.external in plumed_input

        key = self._geometry_to_key(geometry)
        if key not in self.plumed_instances:
            from plumed import Plumed

            tmp = tempfile.NamedTemporaryFile(
                prefix="plumed_", mode="w+", delete=False
            )
            # plumed creates a back up if this file would already exist
            os.remove(tmp.name)
            plumed_ = Plumed()
            ps = 1000 * fs
            plumed_.cmd("setRealPrecision", 8)
            plumed_.cmd("setMDEnergyUnits", mol / kJ)
            plumed_.cmd("setMDLengthUnits", 1 / nm)
            plumed_.cmd("setMDTimeUnits", 1 / ps)
            plumed_.cmd("setMDChargeUnits", 1.0)
            plumed_.cmd("setMDMassUnits", 1.0)
            plumed_.cmd("setLogFile", tmp.name)
            plumed_.cmd("setRestart", True)
            plumed_.cmd("setNatoms", len(geometry))
            plumed_.cmd("init")
            for line in plumed_input.split("\n"):
                plumed_.cmd("readInputLine", line)
            os.remove(tmp.name)  # remove whatever plumed has created
            self.plumed_instances[key] = plumed_

        # input system
        plumed_ = self.plumed_instances[key]
        plumed_.cmd("setStep", 0)
        masses = np.array([atomic_masses[n] for n in geometry.per_atom.numbers])
        plumed_.cmd("setMasses", masses)
        copied_positions = geometry.per_atom.positions.astype(np.float64, copy=True)
        plumed_.cmd("setPositions", copied_positions)
        if geometry.periodic:
            cell = geometry.cell.astype(np.float64, copy=True)
            plumed_.cmd("setBox", cell)

        # perform calculation
        energy = np.zeros((1,))
        forces = np.zeros((len(geometry), 3))
        virial = np.zeros((3, 3))
        plumed_.cmd("setForces", forces)
        plumed_.cmd("setVirial", virial)
        plumed_.cmd("prepareCalc")
        plumed_.cmd("performCalcNoUpdate")
        plumed_.cmd("getBias", energy)
        if geometry.periodic:
            stress = virial / np.linalg.det(geometry.cell)
        else:
            stress = np.zeros((3, 3))
        return {"energy": float(energy.item()), "forces": forces, "stress": stress}

    @staticmethod
    def _geometry_to_key(geometry) -> tuple:
        return tuple([geometry.periodic]) + tuple(geometry.per_atom.numbers)


@typeguard.typechecked
@dataclass
class ZeroFunction(EnergyFunction):
    def __call__(
        self,
        geometry
    ) -> dict[str, float | np.ndarray]:
        return {"energy": 0., "forces": np.zeros(shape=(len(geometry), 3)), "stress": np.zeros(shape=(3, 3))}


@typeguard.typechecked
@dataclass
class HarmonicFunction(EnergyFunction):
    positions: np.ndarray
    hessian: np.ndarray
    energy: Optional[float] = None

    def __call__(
        self,
        geometry
    ) -> dict[str, float | np.ndarray]:
        delta = geometry.per_atom.positions.reshape(-1) - self.positions.reshape(-1)
        grad = np.dot(self.hessian, delta)
        energy = 0.5 * np.dot(delta, grad)
        if self.energy is not None:
            energy += self.energy
        return {"energy": energy, "forces": (-1.0) * grad.reshape(-1, 3), "stress": np.zeros(shape=(3, 3))}

@dataclass
class MACEFunction(EnergyFunction):
    model_path: str
    ncores: int
    device: str
    dtype: str
    atomic_energies: dict[str, float]
    env_vars: Optional[dict[str, str]] = None

    def __post_init__(self):
        import logging

        # import environment variables before any nontrivial imports
        if self.env_vars is not None:
            for key, value in self.env_vars.items():
                os.environ[key] = value

        import torch
        from mace.tools import torch_tools, utils

        torch_tools.set_default_dtype(self.dtype)
        if self.device == "gpu":  # when it's not a specific GPU, use any
            self.device = "cuda"
        self.device = torch_tools.init_device(self.device)

        torch.set_num_threads(self.ncores)
        model = torch.load(f=self.model_path, map_location="cpu")
        if self.dtype == "float64":
            model = model.double()
        else:
            model = model.float()
        model = model.to(self.device)
        model.eval()
        self.model = model
        self.r_max = float(self.model.r_max)
        self.z_table = utils.AtomicNumberTable(
            [int(z) for z in self.model.atomic_numbers]
        )

        # remove unwanted streamhandler added by MACE / torch!
        logging.getLogger("").removeHandler(logging.getLogger("").handlers[0])

    def get_atomic_energy(self, geometry):
        total = 0
        numbers, counts = np.unique(geometry.per_atom.numbers, return_counts=True)
        for idx, number in enumerate(numbers):
            symbol = chemical_symbols[number]
            try:
                total += counts[idx] * self.atomic_energies[symbol]
            except KeyError:
                warnings.warn(f'(MACEFunction) No atomic energy entry for symbol "{symbol}". Are you sure?')
        return total

    def __call__(
        self,
        geometry
    ) -> dict[str, float | np.ndarray]:
        from mace import data
        from mace.tools.torch_geometric.batch import Batch

        # TODO: is this call necessary?
        # torch_tools.set_default_dtype(self.dtype)

        energy, forces, stress = 0.0, np.zeros(shape=(len(geometry), 3)), np.zeros(shape=(3, 3))

        # compute offset if possible
        if self.atomic_energies:
            energy += self.get_atomic_energy(geometry)

        cell = np.copy(geometry.cell) if geometry.periodic else None
        atoms = Atoms(
            positions=geometry.per_atom.positions,
            numbers=geometry.per_atom.numbers,
            cell=cell,
            pbc=geometry.periodic,
        )
        config = data.config_from_atoms(atoms)
        data = data.AtomicData.from_config(config, z_table=self.z_table, cutoff=self.r_max)
        batch = Batch.from_data_list([data]).to(device=self.device)
        out = self.model(batch.to_dict(), compute_stress=cell is not None)
        energy += out["energy"].detach().cpu().item()
        forces += out["forces"].detach().cpu().numpy()
        if cell is not None:
            stress += out["stress"].detach().cpu().numpy().squeeze()
        return {"energy": energy, "forces": forces, "stress": stress}

def function_from_json(path: Union[str, Path], **kwargs) -> Function:
    functions = [
        EinsteinCrystalFunction,
        HarmonicFunction,
        MACEFunction,
        PlumedFunction,
        None,
    ]
    with open(path, "r") as f:
        data = json.loads(f.read())
    assert "function_name" in data
    for function_cls in functions:
        if data["function_name"] == function_cls.__name__:
            break
    data.pop("function_name")
    for name, type_hint in get_type_hints(function_cls).items():
        if type_hint is np.ndarray:
            data[name] = np.array(data[name])
    for key, value in kwargs.items():
        if key in data:
            data[key] = value
    function = function_cls(**data)
    return function

per_atom_dtype = np.dtype(
    [
        ("numbers", np.uint8),
        ("positions", np.float64, (3,)),
        ("forces", np.float64, (3,)),
    ]
)

QUANTITIES = [
    "positions",
    "cell",
    "numbers",
    "energy",
    "per_atom_energy",
    "forces",
    "stress",
    "delta",
    "logprob",
    "phase",
    "identifier",
]

class Geometry:
    """
    Represents an atomic structure with associated properties.

    This class encapsulates the atomic structure, including atom positions, cell parameters,
    and various physical properties such as energy and forces.

    Attributes:
        per_atom (np.recarray): Record array containing per-atom properties.
        cell (np.ndarray): 3x3 array representing the unit cell vectors.
        order (dict): Dictionary to store custom ordering information.
        energy (Optional[float]): Total energy of the system.
        stress (Optional[np.ndarray]): Stress tensor of the system.
        delta (Optional[float]): Delta value, if applicable.
        phase (Optional[str]): Phase information, if applicable.
        logprob (Optional[np.ndarray]): Log probability values, if applicable.
        stdout (Optional[str]): Standard output information, if applicable.
        identifier (Optional[int]): Unique identifier for the geometry.
    """

    per_atom: np.recarray
    cell: np.ndarray
    order: dict
    energy: Optional[float]
    stress: Optional[np.ndarray]
    delta: Optional[float]
    phase: Optional[str]
    logprob: Optional[np.ndarray]
    stdout: Optional[str]
    identifier: Optional[int]

    def __init__(
            self,
            per_atom: np.recarray,
            cell: np.ndarray,
            order: Optional[dict] = None,
            energy: Optional[float] = None,
            stress: Optional[np.ndarray] = None,
            delta: Optional[float] = None,
            phase: Optional[str] = None,
            logprob: Optional[np.ndarray] = None,
            stdout: Optional[str] = None,
            identifier: Optional[int] = None,
        ):
            """
            Initialize a Geometry instance, though the preferred way of instantiating
            proceeds via the `from_data` or `from_atoms` class methods

            Args:
                per_atom (np.recarray): Record array containing per-atom properties.
                cell (np.ndarray): 3x3 array representing the unit cell vectors.
                order (Optional[dict], optional): Custom ordering information. Defaults to None.
                energy (Optional[float], optional): Total energy of the system. Defaults to None.
                stress (Optional[np.ndarray], optional): Stress tensor of the system. Defaults to None.
                delta (Optional[float], optional): Delta value. Defaults to None.
                phase (Optional[str], optional): Phase information. Defaults to None.
                logprob (Optional[np.ndarray], optional): Log probability values. Defaults to None.
                stdout (Optional[str], optional): Standard output information. Defaults to None.
                identifier (Optional[int], optional): Unique identifier for the geometry. Defaults to None.
            """
            self.per_atom = per_atom.astype(per_atom_dtype)  # copies data
            self.cell = cell.astype(np.float64)
            assert self.cell.shape == (3, 3)
            if order is None:
                order = {}
            self.order = order
            self.energy = energy
            self.stress = stress
            self.delta = delta
            self.phase = phase
            self.logprob = logprob
            self.stdout = stdout
            self.identifier = identifier

    def __len__(self):
        """
        Get the number of atoms in the geometry.

        Returns:
            int: The number of atoms.
        """
        return len(self.per_atom)

    def copy(self):
        """
        Create a deep copy of the Geometry instance.

        Returns:
            Geometry: A new Geometry instance with the same data.
        """
        return Geometry.from_string(self.to_string())

    @classmethod
    def load(cls, path_xyz: Union[Path, str]):
        """
        Load a Geometry instance from an XYZ file.

        Args:
            path_xyz (Union[Path, str]): Path to the XYZ file.

        Returns:
            Geometry: A new Geometry instance loaded from the file.
        """
        assert path_xyz.exists()
        with open(path_xyz, "r") as f:
            content = f.read()
        return cls.from_string(content)

    @property
    def periodic(self):
        """
        Check if the geometry is periodic.

        Returns:
            bool: True if the geometry is periodic, False otherwise.
        """
        return np.any(self.cell)
    
    @property
    def volume(self):
        """
        Calculate the volume of the unit cell.

        Returns:
            float: Volume of the unit cell for periodic systems, np.nan for non-periodic systems.
        """
        if not self.periodic:
            return np.nan
        else:
            return np.linalg.det(self.cell)

    @classmethod
    def from_atoms(cls, atoms: Atoms):
        """
        Create a Geometry instance from an ASE Atoms object.

        Args:
            atoms (Atoms): ASE Atoms object.

        Returns:
            Geometry: A new Geometry instance.
        """
        per_atom = np.recarray(len(atoms), dtype=per_atom_dtype)
        per_atom.numbers[:] = atoms.numbers.astype(np.uint8)
        per_atom.positions[:] = atoms.get_positions()
        per_atom.forces[:] = atoms.arrays.get("forces", np.nan)
        if np.any(atoms.pbc):
            cell = np.array(atoms.cell)
        else:
            cell = np.zeros((3, 3))
        geometry = cls(per_atom, cell)
        geometry.energy = atoms.info.get("energy", None)
        geometry.stress = atoms.info.get("stress", None)
        geometry.delta = atoms.info.get("delta", None)
        geometry.phase = atoms.info.get("phase", None)
        geometry.logprob = atoms.info.get("logprob", None)
        geometry.stdout = atoms.info.get("stdout", None)
        geometry.identifier = atoms.info.get("identifier", None)
        return geometry

def create_outputs(quantities: list[str], data: list[Geometry]) -> list[np.ndarray]:
    """
    Create output arrays for specified quantities from a list of Geometry instances.

    Args:
        quantities (list[str]): List of quantity names to extract.
        data (list[Geometry]): List of Geometry instances.

    Returns:
        list[np.ndarray]: List of arrays containing the requested quantities.
    """
    order_names = list(set([k for g in data for k in g.order]))
    assert all([q in QUANTITIES + order_names for q in quantities])
    natoms = np.array([len(geometry) for geometry in data], dtype=int)
    max_natoms = np.max(natoms)
    nframes = len(data)
    nprob = 0
    max_phase = 0
    for state in data:
        if state.logprob is not None:
            nprob = max(len(state.logprob), nprob)
        if state.phase is not None:
            max_phase = max(len(state.phase), max_phase)

    arrays = []
    for quantity in quantities:
        if quantity in ["positions", "forces"]:
            array = np.empty((nframes, max_natoms, 3), dtype=np.float64)
            array[:] = np.nan
        elif quantity in ["cell", "stress"]:
            array = np.empty((nframes, 3, 3), dtype=np.float64)
            array[:] = np.nan
        elif quantity in ["numbers"]:
            array = np.empty((nframes, max_natoms), dtype=np.uint8)
            array[:] = 0
        elif quantity in ["energy", "delta", "per_atom_energy"]:
            array = np.empty((nframes,), dtype=np.float64)
            array[:] = np.nan
        elif quantity in ["phase"]:
            array = np.empty((nframes,), dtype=(np.unicode_, max_phase))
            array[:] = ""
        elif quantity in ["logprob"]:
            array = np.empty((nframes, nprob), dtype=np.float64)
            array[:] = np.nan
        elif quantity in ["identifier"]:
            array = np.empty((nframes,), dtype=np.int64)
            array[:] = -1
        elif quantity in order_names:
            array = np.empty((nframes,), dtype=np.float64)
            array[:] = np.nan
        else:
            raise AssertionError("missing quantity in if/else")
        arrays.append(array)
    return arrays

