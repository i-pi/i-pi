import json
import numpy as np
from typing import Any, Dict, Tuple, List, Union

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
