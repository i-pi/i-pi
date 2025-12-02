PIMD_friction Calculations:
=================================
NVE-f: 
Contains the PIMD simulation with two types of ipi driver 1D_double_well and 1D_double_well_friction. The bafab and bafab in the name of folders demonstrate the inclusion of the friction term in the propagator and the type of splitting.

    DW-NVE: reference calculation without friction in NVE ensemble and 1D double well potential  
    DW-NVE-bafab: calculation with friction using bafab splitting for the propagation step in NVE ensemble and 1D double well potential  
    DW-NVE-fbabf: calculation with friction using fbabf splitting for the propagation step in NVE ensemble and 1D double well potential  
    DWFriction-NVE: reference calculation without friction in NVE ensemble and 1D double well with friction potential 
    DWFriction-NVE-bafab: calculation with friction using bafab splitting for the propagation step in NVE ensemble and 1D double well with friction potential  
    DWFriction-NVE-fbabf: calculation with friction using fbabf splitting for the propagation step in NVE ensemble and 1D double well with friction potential  

    reference-gle-rpmd: Reference result of RPMD calculation with friction implemented in the GLE. https://github.com/GeorgeTrenins/rpmd-gle/tree/main/rpmdgle

    NVE-f.ipynb: Analysis notebook for plotting and comparing the calculation results.



Technical background:
----------------------------

To enable the friction feature in i-PI, the vibrational frequencies and the corresponding spectral density must be provided. These quantities are required to compute the alpha parameters and later the mean-force potentials using the equations described in the reference paper "Non-Markovian Effects in Quantum Rate Calculations of Hydrogen Diffusion with Electronic Friction, PRL,134, 226201 (2025)
" (see Eq. 8a–b in https://doi.org/10.1103/PhysRevLett.134.226201). For a more detailed explanation, please refer to the reference.

$$
\tilde{V}^{\text{ren}}_{N}(\tilde{Q})
= \tilde{V}^{\text{sys}}_{N}(\tilde{Q})
+ \frac{1}{2} \sum_{n} \alpha^{(n)} [\tilde{F}^{(n)}]^2 ,
$$

$$
\alpha^{(n)} =
\frac{2}{\pi} \int_{0}^{\infty}
\frac{\Lambda(\omega)\, \tilde{\omega}_n^{\,2}}
{\omega^{2} + \tilde{\omega}_n^{\,2}}
\, d\omega .
$$

or one can provide alpha as direct input with corresponding normal mode frequencies. 


The ohmic spectral density used as input can be calculated with frictiontool.py (located in ipi/utils/frictiontools.py) by providing the parameters eta and omega_cut (the cutoff frequency).
