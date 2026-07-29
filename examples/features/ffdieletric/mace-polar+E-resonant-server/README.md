# Resonantly driven water with the field applied by i-PI

This example drives a MACE-POLAR-optimized water monomer with

```text
E_z(t) = 0.1 cos(0.719335584407 t) V/angstrom,
```

where `t` is in femtoseconds. The angular frequency corresponds to the
MACE-POLAR symmetric O-H stretch at 3818.84 cm^-1 (period 8.735 fs). The
initial molecular dipole points along -z, so it is antiparallel to the +z
field at `t=0`.

## Files in this folder

- `input.xml`: resonant plane-wave field with `where='server'`.
- `start.extxyz`: water coordinates read by i-PI and MACE; `start.xyz` contains
  the same geometry in plain XYZ form.
- `mace_kwargs.json`: PolarMACE calculator and response-tensor options.
- `download_mace_polar.sh`: downloads `MACE-POLAR-1-M.model` locally.
- `README.md`: frequency, coupling location, and run instructions.

Here `where='server'`: MACE returns the zero-field response tensors and i-PI
adds the field interaction. Download the model and run with

```bash
./download_mace_polar.sh
i-pi input.xml
```

The last three columns of `i-pi.properties.out` contain the applied field in
V/angstrom. This 10 fs trajectory is a compact functionality example rather
than a converged vibrational spectroscopy calculation.
