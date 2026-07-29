# `FFDielectric` examples

This directory contains four MACE-POLAR water examples:

| Folder | Field | Where the coupling is applied |
|---|---|---|
| `mace-polar+E-server` | static | i-PI, using response tensors from MACE |
| `mace-polar+E-client` | static | MACE, after receiving `Efield` from i-PI |
| `mace-polar+E-resonant-server` | resonant plane wave | i-PI, using response tensors from MACE |
| `mace-polar+E-resonant-client` | resonant plane wave | MACE, after receiving `Efield` from i-PI |

Each folder contains an `input.xml`, water starting coordinates, MACE options,
a model-download script, and a README with the run instructions. The examples
use the in-process `extmace` driver, so they do not require a socket client.

## Response tensors

With one or more `<electric_field>` entries and `where='server'`, the driver
must return the dipole, Born effective charges (BECs), and proper piezoelectric
tensor in its `extras` dictionary. The Cartesian order is always `x, y, z`.
Fields can use a built-in function from `ipi.pes.electric_field`, or a custom
Python function selected with the optional `file` attribute.

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
