# `FFDielectric` examples

This directory contains MACE-POLAR water examples:

| Folder | Field | Where the coupling is applied |
|---|---|---|
| `mace-polar+E-server` | static | i-PI, using response tensors from MACE |
| `mace-polar+E-client` | static | MACE, after receiving `electric_field` from i-PI |
| `mace-polar+E-resonant-server` | resonant plane wave | i-PI, using response tensors from MACE |
| `mace-polar+E-resonant-client` | resonant plane wave | MACE, after receiving `electric_field` from i-PI |
| `mace-polar+E-ramp-server` | periodic triangular ramp | i-PI, using response tensors from MACE |
| `mace-polar+E-client-socket` | static | MACE socket client, after receiving `electric_field` through `EXTRADATA` |
| `mace-polar+E-resonant-client-socket` | resonant plane wave | MACE socket client, after receiving `electric_field` through `EXTRADATA` |

## Standalone extxyz evaluation

`ipi.pes.extmace` can also calculate field-dependent energies, forces, and
full stress tensors without starting an i-PI simulation. Store one Cartesian
field vector in the `Atoms.info` dictionary of every extxyz frame, for example:

```text
applied_E="_JSON [0.0, 0.0, 0.1]"
```

Then select that key on the command line:

```bash
python -m ipi.pes.extmace \
  -m MACE-POLAR-1-M.model \
  -mk mace_kwargs.json \
  -i input.extxyz \
  -o output.extxyz \
  --electric-field-key applied_E \
  --field-units V/ang
```

The output contains total `MACE_energy` (eV), `MACE_forces` (eV/angstrom), and
`MACE_stress` (eV/angstrom^3); the input metadata is retained. Each frame may
specify a different field.

For a fixed Cartesian electric displacement, use
`--electric-displacement-key applied_D` instead. The displacement is expressed
in the units selected by `--field-units` (also `V/ang` by default in the
Gaussian convention used here). Fixed-D evaluation additionally needs a
symmetric, nonsingular clamped-ion dielectric tensor, either returned by the
model as `epsilon_infinity` or configured in `mace_kwargs.json`:

```json
{
  "instructions": {
    "epsilon_infinity": [
      [2.0, 0.0, 0.0],
      [0.0, 2.0, 0.0],
      [0.0, 0.0, 2.0]
    ]
  }
}
```

The fixed-D functional is

```text
U_D = U_0 + Omega/(8 pi) (D - 4 pi P)^T epsilon_infinity^-1 (D - 4 pi P),
P = mu/Omega.
```

This is fixed Cartesian D. The volume used by the correction is held constant
when differentiating strain, so the result does not include the Maxwell-stress
term associated with fixed reduced D.

The two key-selection options are mutually exclusive. If neither is present,
the script retains its original zero-field batched behavior.

The in-process examples each contain an `input.xml`, water starting
coordinates, MACE options, a model-download script, and a README. The socket
examples reuse the corresponding client example's model assets and launch
`i-pi-py_driver`; its `requires_extra=true` parameter makes it advertise the
opt-in `NEEDEXTRA` capability.

## Client-side field acknowledgement

When `where='client'`, the driver must include every field it applied in its
returned extras JSON. For an electric field, return

```json
{"applied_fields": ["electric_field"]}
```

i-PI stops with an error if a field it sent is not acknowledged, preventing a
driver that silently ignores the field from producing an incorrect trajectory.

## Response tensors

With one or more `<electric_field>` entries and `where='server'`, the driver
must return the dipole, Born effective charges (BECs), and proper piezoelectric
tensor in its `extras` dictionary. The Cartesian order is always `x, y, z`.
Fields can use a built-in function from `ipi.pes.electric_field`, or a custom
Python function selected with the optional `file` attribute.

The built-in `ramp` field accepts `amplitude`, `period`, and `period_units`
(default: `atomic_unit`). In every period it follows a linear
`0 -> +amplitude -> 0 -> -amplitude -> 0` profile (up, down, down, up).
`phase` is accepted but ignored.

## Expected data

The default keys, shapes, axes, and units are:

| Key | Shape | Axis order | Default units |
|---|---:|---|---|
| `dipole` | `(3,)` | dipole component `i` | `eang` |
| `BEC` | `(natoms, 3, 3)` | atom `a`, dipole `i`, displacement `j` | `e` |
| `piezoelectric` | `(3, 3, 3)` | dipole `i`, strain `j`, strain `k` | `e/ang2` |

Equivalently,

```text
BEC[a, i, j]              = d mu_i / d R_(a,j)
piezoelectric[i, j, k]     = e_proper[i, j, k]
```

i-PI converts these quantities to atomic units and adds

```text
Delta U                  = -sum_i mu_i E_i
Delta F[a, j]            =  sum_i BEC[a, i, j] E_i
Delta stress[j, k]       = -sum_i e_proper[i, j, k] E_i
Delta virial[j, k]       =  Omega sum_i e_proper[i, j, k] E_i.
```

The last equality follows from i-PI's convention
`virial = -Omega * stress`. The proper piezoelectric tensor must satisfy

```text
e_proper[i, j, k] = e_proper[i, k, j].
```

It must be provided as the full Cartesian `(3, 3, 3)` tensor, not in Voigt
notation.

The default extractors may be written explicitly as

```xml
<dipole units="eang"   key="dipole" />
<bec    units="e"      key="BEC" />
<piezo  units="e/ang2" key="piezoelectric" />
```

These lines can be omitted when the driver uses the default keys and units.

## Comparison with driven dynamics

The [driven-dynamics electric-field example](../driven_dynamics/efield-water/README.md)
uses the same physical BEC but a different array layout:

```text
BEC_driven[3*a + j, i] = (1/q_e) d mu_i / d R_(a,j).
```

Driven dynamics therefore uses coordinate rows and dipole columns, whereas
`FFDielectric` uses dipole rows and coordinate columns within each atom. The
conversion is

```python
# (3*natoms, 3) driven-dynamics layout -> (natoms, 3, 3) FFDielectric layout
bec_ff = bec_driven.reshape(natoms, 3, 3).transpose(0, 2, 1)

# FFDielectric -> driven dynamics
bec_driven = bec_ff.transpose(0, 2, 1).reshape(3 * natoms, 3)
```

Thus the physical convention is the same, but the stored per-atom matrices
are transposed. When driven-dynamics values are dimensionless and
`FFDielectric` uses `units="e"`, their numerical values otherwise agree.
