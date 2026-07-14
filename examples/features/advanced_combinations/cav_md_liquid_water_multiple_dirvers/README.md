Cavity molecular dynamics of liquid water under vibrational strong coupling
=======================

This example extends the [cav_md_liquid_water](../cav_md_liquid_water/README.md) setup for liquid water under vibrational strong coupling. As in the base example, an infrared cavity mode at 3550 cm<sup>-1</sup> is strongly coupled to the O-H stretch mode of liquid water, but here the CavMD simulation is generalized to propagate under advanced dipole surface drivers beyond fixed point charges.

## More advanced dipole surfaces

A more advanced way to propagate CavMD is to use more accurate dipole surfaces. The following input file and bash script use a dipole-induced-dipole water dipole surface (`water_dip_pol`, which is supported in i-PI):
- **input_multiple_drivers.xml**
- **run_multiple_drivers.sh**

Here, we need to make **two modifications** of conventional **ffcavphsocket** i-pi simulations.
- 1. **evaluate_photon**: determines whether to evaluate cavity forces;
- 2. **dipole_surface**: determines whether to use dipoles and dipole derivates from extenal driver, e.g. `water_dip_pol`.


In this case, i-PI is connected to two separate **ffcavphsocket** force evaluators: the **qtip4pf** driver provides the nuclear forces, while the **dipole** driver provides the dipole and dipole derivatives used for CavMD simulations. In this setup, nuclear forces can come from any conventional MD driver supported by i-PI; and users need to provide or implement a dipole driver for accurate dipole evaluations.

In practice, the input file reads as follows:
```xml
  <ffcavphsocket name='qtip4pf' mode='unix' pbc='False'>
    <address>h2o-cl-cavmd</address>
    <latency> 1e-3 </latency>
    <apply_photon> True </apply_photon>
    <evaluate_photon> False </evaluate_photon>
    <dipole_surface> False </dipole_surface>
    <E0> 4e-4 </E0>
    <omega_c units='inversecm'> 3550.0 </omega_c>
  </ffcavphsocket>
  <ffcavphsocket name='dipole' mode='unix' pbc="false">
    <address> h2o-dipole-cl-cavmd </address>
    <latency> 1e-3 </latency>
    <apply_photon> True </apply_photon>
    <evaluate_photon> True </evaluate_photon>
    <dipole_surface> True </dipole_surface>
    <E0> 4e-4 </E0>
    <omega_c units='inversecm'> 3550.0 </omega_c>
  </ffcavphsocket>
  <system>
    <forces>
      <force forcefield='qtip4pf'> </force>
      <force forcefield='dipole'> </force>
    </forces>
  </system>
```

The accompanying **run_multiple_drivers.sh** script starts i-PI first and then connects the two drivers to their corresponding sockets:

```bash
i-pi-driver -m qtip4pf -u -a h2o-cl-cavmd > /dev/null &
i-pi-driver -u -h h2o-dipole-cl-cavmd -m water_dip_pol -o 0 &> out
```
