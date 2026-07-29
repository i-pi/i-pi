# BaTiO3, MACE-POLAR, and a homogeneous electric field

This example exercises i-PI's in-process MACE driver, MACE-POLAR dielectric
response, and the `FFDielectric` fixed-electric-field wrapper. No socket client
is needed. The calculation uses a local `MACE-POLAR-1-M.model` checkpoint.

MACE-POLAR-1 was trained on molecular data from OMol25. The periodic BaTiO3
cell used here is therefore a wiring and smoke test, not a validated BaTiO3
potential. In addition, a total dipole is branch-dependent under periodic
boundary conditions. Do not interpret this short trajectory as a physical
prediction for BaTiO3.

## Install

Use one Python environment for MACE and i-PI. MACE-POLAR support requires a
recent MACE version and its long-range graph dependency:

```bash
python -m pip install -U 'mace-torch>=0.3.16'
python -m pip install 'git+https://github.com/WillBaldwin0/graph_electrostatics.git'
python -m pip install -e ../../../..
```

The MACE-POLAR checkpoints use the Academic Software License (ASL). Review and
accept its terms before downloading a checkpoint.

## Download and run

From this directory:

```bash
./download_mace_polar.sh
i-pi input.xml
```

The download script writes `MACE-POLAR-1-M.model`. The relevant force-field
part of `input.xml` is:

```xml
<ffdielectric name='mu+E' mode='E' where='server'>
  <field mode='pw'>
    <pw>
      <amplitude units='V/ang'>[0, 0, 1]</amplitude>
      <freq units='THz'>0</freq>
    </pw>
  </field>
  <ffdirect name='mace' pbc='True'>
    <pes>extmace</pes>
    <parameters>{
      template:start.extxyz,
      model:MACE-POLAR-1-M.model,
      device:cpu,
      mace_kwargs:mace_kwargs.json
    }</parameters>
  </ffdirect>
</ffdielectric>
```

`where='server'` means that MACE returns the zero-field energy, forces,
stress, dipole, Born effective charges, and piezoelectric response. i-PI then
adds the field-dependent contribution. With `where='client'`, the underlying
driver would instead be responsible for applying the field.

`mace_kwargs.json` sets `model_type` to `PolarMACE`. It also enables response
derivatives inside the `instructions` dictionary:

```json
{
    "model_type": "PolarMACE",
    "default_dtype": "float32",
    "instructions": {
        "compute_BEC": true,
        "ignore": []
    }
}
```

The placement of `compute_BEC` is important: it is an instruction to the i-PI
MACE wrapper, not a top-level `MACECalculator` argument.

Use `device:cuda` for a CUDA GPU. Float32 is appropriate for this MD smoke
test. For response validation or finite-difference comparisons, use
`default_dtype: float64`.

`start.extxyz` supplies the global electronic state expected by MACE-POLAR:
neutral charge, singlet spin multiplicity, and zero external field. Change
these metadata for a different charge or spin state.

## Notation and tensor layout

The following notation is used below:

- `a` labels an atom and `b` labels a structure in a batch.
- `i` labels a dipole or electric-field Cartesian component.
- `j,k` label Cartesian position or deformation components.
- `mu_i` is a cell dipole and `Omega` is the cell volume.
- `eta_jk` is the full infinitesimal deformation gradient.
- `epsilon_jk = (eta_jk + eta_kj)/2` is symmetric strain.

Before removal of the batch dimension, the calculated tensors have layouts

```text
BEC[i, a, j]              = d mu_i / d R_(a,j)
dmu_deta[i, b, j, k]      = d corrected_mu_i / d eta_(b,j,k)
piezoelectric[b, i, j, k] = dmu_deta[i, b, j, k] / Omega_b
```

The BEC sent to `FFDielectric` is rearranged to `(natoms, 3, 3)`, with axes
`(atom, dipole_component, position_component)`. The piezoelectric tensor is
stored per structure with shape `(3, 3, 3)`.

## How MACE introduces strain

MACE creates a zero-valued, unconstrained autograd leaf called
`displacement`, here denoted by `eta`. It then constructs

```python
epsilon = 0.5 * (eta + eta.transpose(-1, -2))
```

and deforms positions and the cell as

```text
R' = R + R epsilon
h' = h + h epsilon.
```

Only `epsilon` enters the model geometry. The raw `eta` leaf is returned in
`data["displacement"]` so that PyTorch can differentiate with respect to it.
Consequently, the raw model dipole satisfies

```text
mu(eta) = mu(sym(eta)),
```

and is insensitive to an antisymmetric perturbation of `eta`. Its derivative
therefore obeys

```text
D^epsilon_(i,j,k) = d mu_i / d eta_jk
                  = D^epsilon_(i,k,j).
```

This distinction matters: `data["displacement"]` is not itself the
symmetrized tensor, even though the model dipole depends on it only through
its symmetric part. A later expression that uses `data["displacement"]`
directly must preserve that convention explicitly.

## Improper and proper piezoelectric response

Vanderbilt defines the improper response using the full deformation gradient:

```text
c_(i,j,k) = d P_i / d eta_jk,
```

where `P_i = mu_i/Omega` is the polarization. The proper response is

```text
e^proper_(i,j,k)
    = c_(i,j,k) + delta_jk P_i - delta_ij P_k.
```

See D. Vanderbilt, *Berry-phase theory of proper piezoelectric response*,
J. Phys. Chem. Solids 61, 147-151 (2000),
<https://arxiv.org/abs/cond-mat/9903137>.

Since

```text
d Omega / d eta_jk = Omega delta_jk,
```

the extensive form of Vanderbilt's expression is

```text
Omega e^proper_(i,j,k)
    = D^eta_(i,j,k) - delta_ij mu_k,
```

where `D^eta` is the dipole derivative with respect to the full deformation,
including rotations.

The full derivative is rotationally covariant. Its antisymmetric part is

```text
D^eta_(i,j,k) - D^eta_(i,k,j)
    = delta_ij mu_k - delta_ik mu_j.
```

This rotational contribution exactly cancels the antisymmetric part of the
geometrical correction. The physical proper tensor is therefore symmetric in
its mechanical indices:

```text
e^proper_(i,j,k) = e^proper_(i,k,j).
```

MACE directly provides only the symmetric-strain derivative `D^epsilon`, not
the rotational part of `D^eta`. Rotational covariance gives

```text
D^eta_(i,j,k)
    = D^epsilon_(i,j,k)
      + 1/2 (delta_ij mu_k - delta_ik mu_j).
```

Substitution into Vanderbilt's expression gives the formula used here:

```text
e^proper_(i,j,k) = (1 / Omega) * [
    D^epsilon_(i,j,k)
    - 1/2 (delta_ij mu_k + delta_ik mu_j)
].
```

It is manifestly symmetric under `j <-> k`.

## Why `proper_dipole` uses symmetric strain

The code obtains the preceding correction by defining, to first order around
zero strain,

```text
corrected_mu_i = mu_i - sum_k epsilon_ik mu_k,
```

with `epsilon = sym(eta)`. Differentiating with respect to the raw autograd
leaf gives

```text
d corrected_mu_i / d eta_jk
    = D^epsilon_(i,j,k)
      - 1/2 (delta_ij mu_k + delta_ik mu_j).
```

Dividing this extensive derivative by `Omega` produces the proper
piezoelectric tensor in units of charge per area. The explicit
symmetrization inside `proper_dipole` is required because MACE returns the raw
`eta` leaf, even though its model outputs depend only on `sym(eta)`.

At zero strain, `corrected_mu == mu`. Therefore the correction does not alter
the dipole value or the BEC. It changes only the strain derivative needed for
the proper piezoelectric response.

## Volume normalization and units

`dmu_deta` is a dipole derivative with respect to a dimensionless strain, so
it is extensive and has dipole units. The piezoelectric tensor is intensive:

```text
piezoelectric = dmu_deta / abs(det(cell)).
```

For MACE outputs in electron-Angstrom and cells in Angstrom, this converts
`e Ang` to `e/Ang^2`, which is the unit expected by the default
`FFDielectric` piezoelectric extractor.

The volume normalization must happen once. Do not divide the tensor by the
volume again in `FFDielectric` or in output post-processing.

## Field coupling in `FFDielectric`

For fixed electric field, i-PI uses the model response to add terms of the
form

```text
U_E = -mu_i E_i,
F_(a,j),E = Z*_(a,i,j) E_i.
```

The piezoelectric response supplies the cell/virial derivative of the field
coupling. With the proper piezoelectric tensor, the field-induced stress and
the corresponding i-PI virial are

```text
sigma_E[j, k] = -sum_i piezoelectric[i, j, k] * Efield[i]
virial_E[j, k] = -Omega * sigma_E[j, k]
                = Omega * sum_i piezoelectric[i, j, k] * Efield[i].
```

The plus sign in the virial follows from i-PI's `virial = -Omega * stress`
convention. The factor `Omega` restores an extensive, energy-like virial from
the intensive piezoelectric tensor. It is applied exactly once in
`FFDielectric`; `_mace.py` has already divided the dipole-strain derivative by
`Omega` when constructing `piezoelectric`.

## Selecting MACE outputs

MACE-POLAR can return many optional tensors in addition to energy, forces,
stress, dipole, BEC, and piezoelectric response. The exact set varies between
MACE versions and checkpoints. Transferring every tensor from a GPU is
unnecessary, so `mace_kwargs.json` ignores optional outputs not used here.

At the first force evaluation, i-PI prints all model outputs. Registered
outputs are grouped into per-atom and per-structure properties; unregistered
outputs are printed with their raw shapes. If unsupported outputs remain,
i-PI reports them together and suggests JSON registrations and ignore entries.

To ignore outputs, merge their names into the existing list:

```json
{
    "instructions": {
        "ignore": [
            "atomic_embedding",
            "graph_tensor",
            "spin_density"
        ]
    }
}
```

To retain an output, register its shape at the top level and do not ignore it:

```json
{
    "model_type": "PolarMACE",
    "default_dtype": "float32",
    "ase_like_properties": {
        "atomic_embedding": ["natoms", 16],
        "graph_tensor": [3, 3],
        "spin_density": ["natoms", 2, 4]
    },
    "instructions": {
        "compute_BEC": true,
        "ignore": ["other_unused_property"]
    }
}
```

`"natoms"` marks a per-atom leading dimension. A shape without `"natoms"`
is per-structure; i-PI removes the leading batch dimension. Shape suggestions
are inferred from dimensions alone, so verify their physical meaning.

## Diagnostics and common configuration errors

### `NoneType` has no attribute `actual_time`

`FFDielectric` evaluates a time-dependent field using the system motion's
`actual_time`. Force evaluation uses a contracted `Beads` object even when no
ring-polymer contraction is requested. That object must retain the parent
motion reference.

### `The key 'BEC' is not in the provided dictionary`

Check that `"compute_BEC": true` is inside `"instructions"` in
`mace_kwargs.json`. Also confirm that the output summary lists `BEC` among the
copied per-atom properties.

### BEC has the wrong rank

`FFDielectric` expects `(natoms, 3, 3)`, not a flattened `(natoms, 9)` array.
The extended MACE property registration preserves the full tensor shape.

### Dipole-strain derivative is not symmetric

The code reports the Frobenius norm of

```text
dmu_deta - transpose_jk(dmu_deta).
```

A large value usually means that an unsymmetrized `eta` was used directly in
the dipole correction, or that a model path depends on the antisymmetric part
of the deformation. For float32, tiny residuals can also arise from numerical
roundoff; repeat the response calculation in float64 before changing the
tolerance.

## Recommended validation

Before using this machinery for production calculations:

1. Run the response evaluation in float64.
2. Verify `BEC[i,a,j]` against finite differences of `mu_i` with respect to
   atomic displacement `R_(a,j)`.
3. Verify `e^proper_(i,j,k)` against finite symmetric cell strains.
4. Check symmetry under `j <-> k` and the expected crystal point-group
   relations.
5. Check energy, force, and virial field derivatives against finite
   differences using the same unit and sign conventions.
6. Establish that the chosen polarization branch is meaningful for the
   periodic system under study.
