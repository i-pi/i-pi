Cavity molecular dynamics of liquid water under vibrational strong coupling
=======================

This runs a cavity molecular dynamics (CavMD) simulation of liquid water under vibrational strong coupling. 

In this simulation, an infrared cavity mode at 3550 cm<sup>-1</sup> is strongly coupled to the O-H stretch mode of liquid water along both the x and y directions. Due to this strong coupling, the O-H stretch band of liquid water is splitted to two eigenstates, known as polaritons.

See the paper [Proc. Natl. Acad. Sci., 2020, 117 (31), 18324–18331](https://doi.org/10.1073/pnas.2009272117) or [https://github.com/TaoELi/cavity-md-ipi](https://github.com/TaoELi/cavity-md-ipi) for more details of CavMD.

This example is the same as the [Rabi splitting tutorial](https://github.com/TaoELi/cavity-md-ipi/tree/master/tutorials/Rabi_splitting) in [https://github.com/TaoELi/cavity-md-ipi](https://github.com/TaoELi/cavity-md-ipi). The only difference is that in the original [Rabi splitting tutorial](https://github.com/TaoELi/cavity-md-ipi/tree/master/tutorials/Rabi_splitting), the parameters for the photons are controlled by a separate input file **photon_params.json**. Here, instead, the photon parameters are controlled in the **ffcavphsocket** section of the i-pi input file **input.xml**. 

## Fixed point charges as the dipole surface

The simplest CavMD simulations assume that the molecular dipole moments are calculated using a linear combination of fixed point charges times the nuclear coordinates. For this case, the following input file and bash script are used for simulations:
- **input.xml**
- **run.sh**

Here, we need to make **two modifications** of conventional i-pi simulations.

1. We use the following force evaluator **ffcavphsocket** instead of the original **ffsocket**. In this new force evaluator, **E0** denotes the effective light-matter coupling strength in atomic units (corresponding to the \tilde{varepsilon} parameter in the [original paper](https://doi.org/10.1073/pnas.2009272117)). **omega_c_cminv** denotes the cavity photon mode frequency in wave number. **charge_array** defines the partial charges of all atoms in the same order as the input xyz file. Note that the total size of **charge_array** should match the total number of nuclei, otherwise a **ValueError** will be raised.

```
  <ffcavphsocket name='lammps' mode='unix' pbc='True'>
    <address>h2o-cl-cavmd</address>
    <latency> 1e-3 </latency>
    <apply_photon> True </apply_photon>
    <E0> 4e-4 </E0>
    <omega_c units='inversecm'> 3550.0 </omega_c>
    <charge_array> [-0.8192, 0.4096, 0.4096, ...] </charge_array>
  </ffcavphsocket>
```

2. We also need to put photons in the input xyz file.

```
650

       O  2.38972e+00  2.62791e-01  5.65120e+00
       H  1.41821e+00  1.87797e-01  5.65303e+00
       H  2.50853e+00  2.13428e-01  4.67340e+00
       ...
       L -8.65101e+00  1.11541e+00  4.56823e-01
       L  5.35376e-01  1.20389e+01 -1.19497e-01
```
Here, we have 216 water molecules (648 nuclei) and two photons (labeled as **L**). Please always put the photons at the very end of the xyz file. The mass of each photon is 1.0 atomic units by default. We assume a z-direction oriented optical cavity coupled to the molecules, so the k-vector of the cavity photons is along the z-direction, and the cavity photons can polarize along either the x- or y-direction. The first photon is coupled to the x direction of the total dipole moment of water, and the second photon is coupled to the y direction. Hence, for the first photon "atom", only the x coordinate is influenced by the molecules; for the second photon "atom", only the y coordinate is influenced by the molecules. In the current implementation, only two cavity photons are supported (more advanced cavity setups will be included in an updated version of the code).

## More advanced dipole surfaces

A more advanced way to propagate CavMD is to use more accurate dipole surfaces. The following input file and bash script use a dipole-induced-dipole water dipole surface (`water_dip_pol`, which is supported in i-PI):
- **input_multiple_drivers.xml**
- **run_multiple_drivers.sh**

In this case, i-PI is connected to two separate **ffcavphsocket** force evaluators: the **qtip4pf** driver provides the nuclear forces, while the **dipole** driver provides the dipole and dipole derivatives used for CavMD simulations. In this setup, nuclear forces can come from any conventional MD driver supported by i-PI; and users need to provide or implement a dipole driver for accurate dipole evaluations.

In practice, the input file reads as follows:
```xml
  <ffcavphsocket name='qtip4pf' mode='unix' pbc='False'>
    <address>h2o-cl-cavmd</address>
    <latency> 1e-3 </latency>
    <apply_photon> True </apply_photon>
    <E0> 4e-4 </E0>
    <omega_c units='inversecm'> 3550.0 </omega_c>
    <charge_array> [0.0, 0.0, 0.0, ...] </charge_array>
  </ffcavphsocket>
  <ffcavphsocket name='dipole' mode='unix' pbc="false">
    <address> h2o-dipole-cl-cavmd </address>
    <latency> 1e-3 </latency>
    <apply_photon> True </apply_photon>
    <E0> 4e-4 </E0>
    <omega_c units='inversecm'> 3550.0 </omega_c>
    <charge_array> [] </charge_array>
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
