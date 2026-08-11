# Resonantly driven water with the field applied by MACE

This example drives a MACE-POLAR-optimized water monomer with

```text
E_z(t) = 0.1 cos(2*pi*114485.814 GHz*t) V/angstrom,
```

where `t` is in atomic units. The cyclic frequency corresponds to the
MACE-POLAR symmetric O-H stretch at 3818.84 cm^-1 (period 8.735 fs). The
initial molecular dipole points along -z, so it is antiparallel to the +z
field at `t=0`.

## Files in this folder

- `input.xml`: resonant plane-wave field with `where='client'`.
- `start.extxyz`: water coordinates read by i-PI and MACE; `start.xyz` contains
  the same geometry in plain XYZ form.
- `mace_kwargs.json`: PolarMACE options for applying the field in MACE.
- `download_mace_polar.sh`: downloads `MACE-POLAR-1-M.model` locally.
- `README.md`: frequency, coupling location, and run instructions.

Here `where='client'`: i-PI sends `Efield` to `extmace`, and MACE returns the
field-dependent energy, forces, and stress. Download the model and run with

```bash
./download_mace_polar.sh
i-pi input.xml
```

The last three columns of `i-pi.properties.out` contain the applied field in
V/angstrom. This 10 fs trajectory is a compact functionality example rather
than a converged vibrational spectroscopy calculation.
