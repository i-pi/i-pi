import json
import numpy as np
from typing import Any, Dict, Tuple, List, Union
from ipi.utils.units import unit_to_internal, unit_to_user


# --------------------------------------- #
class JSONLogger:

    def __init__(self, file: str = None):
        self.file = file
        self.enabled = file is not None

    def _serialize(self, obj: Any) -> Any:
        """Recursively convert objects into JSON-serializable types."""

        # NumPy arrays → lists
        if isinstance(obj, np.ndarray):
            return obj.tolist()

        # NumPy scalar types → Python scalars
        if isinstance(obj, (np.floating, np.integer)):
            return obj.item()

        # Standard JSON scalars
        if isinstance(obj, (int, float, str, bool)) or obj is None:
            return obj

        # Dicts → recursively serialize values
        if isinstance(obj, dict):
            return {k: self._serialize(v) for k, v in obj.items()}

        # Lists / tuples / sets → recursively process each element
        if isinstance(obj, (list, tuple, set)):
            return [self._serialize(x) for x in obj]

        # Fallback: store as string if unrecognized
        return str(obj)

    def save(self, results: Dict[str, Any], file: str = None):
        if not self.enabled:
            return

        safe_data = self._serialize(results)

        if file is None:
            file = self.file
        # ToDo: change it
        with open(file, "a") as f:
            f.write(json.dumps(safe_data, indent=4))
            f.write("\n")


# --------------------------------------- #
Parent = Dict[str, Union[float, np.ndarray]]


class StructureResults(dict):
    """
    Stores the results produced by one model for one structure.

    Parameters
    ----------
    natoms : int
        Number of atoms in this structure. Used to expand dynamic
        shape entries equal to the string "natoms".
    shapes : Dict[str, Tuple]
        Desired output shapes for each property. Tuples may contain the
        literal "natoms", which is replaced by the actual integer.
    """

    def __init__(self, natoms: int, shapes: Dict[str, Tuple]):
        self.natoms = natoms
        # Expand "natoms" once at initialization
        self.shapes = {k: self._expand(shape) for k, shape in shapes.items()}

    def _expand(self, shape: Tuple) -> Tuple:
        """Replace 'natoms' in shape by the actual integer."""
        if not isinstance(shape, tuple):
            raise TypeError(f"Expected tuple, got {type(shape)}")
        return tuple(self.natoms if d == "natoms" else d for d in shape)

    def store(self, key: str, value: Union[float, np.ndarray]):
        """Store a property value, reshaping it to the expected shape."""
        if key not in self.shapes:
            raise KeyError(f"Unknown property '{key}'")

        shape = self.shapes[key]
        if shape == ():  # scalar
            try:
                self[key] = float(value)
            except Exception:
                raise ValueError(f"Value for '{key}' must be convertible to float")
        else:
            try:
                self[key] = np.asarray(value).reshape(shape)
            except Exception:
                raise ValueError(
                    f"Value for '{key}' has wrong shape: expected {shape}, got {np.asarray(value).shape}"
                )

    def __add__(self, other: "StructureResults") -> "StructureResults":
        out = StructureResults(self.natoms, self.shapes)
        out.update({k: self[k] + other[k] for k in self})
        return out

    def __truediv__(self, divisor: float) -> "StructureResults":
        out = StructureResults(self.natoms, self.shapes)
        out.update({k: self[k] / divisor for k in self})
        return out

    def as_dict(self, copy_arrays: bool = True) -> Parent:
        """Return stored properties as a plain dictionary, optionally copying arrays."""
        return {
            k: (v.copy() if copy_arrays and isinstance(v, np.ndarray) else v)
            for k, v in self.items()
        }


class ModelResults:
    """
    Handle results returned by one model for multiple structures.
    """

    def __init__(self, shapes: Dict[str, Tuple]):
        self._shapes = shapes
        self._results: List[StructureResults] = []

    def store(self, natoms: List[int], results: Dict[str, Any]):
        """
        Store results of one model over multiple structures.
        """
        unknown = {
            key: value for key, value in results.items() if key not in self._shapes
        }
        if unknown:
            raise ValueError(self._unknown_properties_message(unknown, natoms))

        ptr = np.cumsum([0] + natoms)
        new_structs = [StructureResults(n, self._shapes) for n in natoms]

        for key, value in results.items():
            if "natoms" in self._shapes[key]:
                value = np.split(value, ptr[1:], axis=0)[:-1]

            for i, s in enumerate(new_structs):
                s.store(key, value[i])

        self._results.extend(new_structs)

    @staticmethod
    def _shape_candidates(value: Any, natoms: List[int]) -> Tuple[Tuple, List[List]]:
        """Return the raw tensor shape and plausible settings-file shapes."""

        raw_shape = tuple(np.asarray(value).shape)
        shape_candidates = []

        if raw_shape and raw_shape[0] == sum(natoms):
            shape_candidates.append(["natoms", *raw_shape[1:]])
        if raw_shape and raw_shape[0] == len(natoms):
            per_structure = list(raw_shape[1:])
            if per_structure not in shape_candidates:
                shape_candidates.append(per_structure)

        return raw_shape, shape_candidates

    @classmethod
    def _unknown_properties_message(
        cls, unknown: Dict[str, Any], natoms: List[int]
    ) -> str:
        """Report every unregistered output and consolidated JSON remedies."""

        details = []
        registrations = {}
        ambiguous = []

        for key in sorted(unknown):
            raw_shape, candidates = cls._shape_candidates(unknown[key], natoms)
            if len(candidates) == 1:
                inferred = candidates[0]
                registrations[key] = inferred
                kind = (
                    "per-atom"
                    if inferred and inferred[0] == "natoms"
                    else "per-structure"
                )
                details.append(
                    f"- '{key}': raw shape {raw_shape}; inferred {kind} shape "
                    f"{json.dumps(inferred)}"
                )
            elif len(candidates) > 1:
                ambiguous.append(key)
                choices = " or ".join(json.dumps(shape) for shape in candidates)
                details.append(
                    f"- '{key}': raw shape {raw_shape}; ambiguous shape, choose "
                    f"{choices} based on the property's semantics"
                )
            else:
                inferred = list(raw_shape)
                registrations[key] = inferred
                details.append(
                    f"- '{key}': raw shape {raw_shape}; no atom/batch leading "
                    f"dimension detected, suggested shape {json.dumps(inferred)} "
                    "(verify manually)"
                )

        registration = json.dumps({"ase_like_properties": registrations}, indent=2)
        ignore = json.dumps({"instructions": {"ignore": sorted(unknown)}}, indent=2)
        manual_note = ""
        if ambiguous:
            manual_note = (
                "\nAmbiguous properties are omitted from the consolidated "
                "registration block; add one of the shapes shown above manually."
            )

        return (
            f"Unknown model properties for a batch with natoms={natoms}:\n"
            + "\n".join(details)
            + "\nTo retain the unambiguous outputs, merge this block into the "
            "MACE settings JSON and verify the inferred shapes:\n"
            + registration
            + manual_note
            + "\nIf none of these outputs is needed, merge this ignore block "
            "instead:\n" + ignore
        )

    def __len__(self) -> int:
        return len(self._results)

    def natoms(self) -> List[int]:
        return [s.natoms for s in self._results]

    def shapes(self) -> List[Dict[str, Tuple]]:
        return [s.shapes for s in self._results]

    @staticmethod
    def mean(models: List["ModelResults"]) -> "ModelResults":
        if not models:
            raise ValueError("Cannot compute mean of empty list")

        ref_natoms = models[0].natoms()
        ref_shapes = models[0].shapes()
        shapes = models[0]._shapes
        n_structures = len(ref_natoms)
        n_models = len(models)

        # Validate consistency
        for m in models[1:]:
            if m.natoms() != ref_natoms or m.shapes() != ref_shapes:
                raise ValueError(
                    "All ModelResults must have the same natoms and shapes per structure"
                )

        out = ModelResults(shapes)
        for s in range(n_structures):
            # Initialize with a zeroed StructureResults
            summed = StructureResults(ref_natoms[s], shapes)
            for key in models[0]._results[s].keys():
                summed[key] = np.zeros(summed.shapes[key])

            # Sum over models
            for m in models:
                summed += m._results[s]

            out._results.append(summed / n_models)

        return out

    def __getitem__(self, i: int) -> Parent:
        return self._results[i].as_dict()


class Instructions:

    dimensions = {}
    units = {"length": "atomic_unit", "energy": "atomic_unit"}

    def __init__(self, instructions: dict, *argc, **argv):

        super().__init__(*argc, **argv)

        # Read instructions from file or use provided dictionary
        if isinstance(instructions, str):
            with open(instructions, "r") as f:
                instructions = json.load(f)
        elif isinstance(instructions, dict):
            pass
        else:
            raise ValueError("`instructions` can be `str` or `dict` only.")

        # convert parameters to the required units
        to_delete = list()
        for k, dimension in self.dimensions.items():
            variable = f"{k}_unit"
            if variable in instructions:
                to_delete.append(variable)
                if instructions[variable] is not None:
                    factor = convert(
                        1,
                        dimension,
                        _from=instructions[variable],
                        _to=self.units[dimension],
                    )
                    instructions[k] = process_input(instructions[k]) * factor
        for k in to_delete:
            del instructions[k]

        self.instructions = instructions


# ---------------------- #
def convert(
    what: float,
    family: str = None,
    _from: str = "atomic_unit",
    _to: str = "atomic_unit",
) -> float:
    """
    Converts a physical quantity between units of the same type (length, energy, etc.)
    Example:
    value = convert(7.6,'length','angstrom','atomic_unit')
    arr = convert([1,3,4],'energy','atomic_unit','millielectronvolt')
    """
    # from ipi.utils.units import unit_to_internal, unit_to_user
    if family is not None:
        factor = unit_to_internal(family, _from, 1)
        factor *= unit_to_user(family, _to, 1)
        return what * factor
    else:
        return what


# ---------------------- #
def process_input(value):
    """
    Standardizes user input into numerical format (float or np.array).
    """
    if isinstance(value, float):
        return value
    elif isinstance(value, int):
        return value
    elif isinstance(value, list):
        return np.array(value)
    else:
        raise TypeError("Input must be a float or a list.")
