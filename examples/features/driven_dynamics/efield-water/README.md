# Light-driven water molecule
This example provide a way to excite a vibrational mode of an **isolated water molecule** by using a **THz electric field**.

## The physics behind this example
This example perform an ```NVE``` simulation of an isolated water molecule.
The provided initial configuration (```start.xyz```) is the ground state of the system, and no initial velocities are provided.
Then, by integrating the Hamilton's equations of motion, we expect the molecule to stay in the initial configuration (forever). 

This example we added an external (time-dependent, oscillating) electric field, which couples to the **dipole** $\mathbf{d}$ of the system (according to the **Electric Dipole Approximation**) and it can then excite the **vibrational modes** of the systems (if the field is in resonance with one of the modes, and the mode is actually IR active).

By running this simulation you will then see the water molecule having one of its vibrational mode getting excited as long as the applied electric field intensity increases.

As you might have already understood, this approach relies on the approximation that the coupling with the field happens <ins>"at the molecular dynamics level"</ins> only: 
- the only contribution from the electrons come from the value of the dipole $\mathbf{d}$
- the electric field is NOT inserted into the DFT equations (or any other infrastructure to compute the energy and the forces) 
- the electrons are assumed to stay in their unperturbed (Born-Oppenheimer) groud-state along the whole trajectory.

<ins>To summarize</ins>: any physical quantity (energy, forces, stresses, dipole, etc ...) only depend on the nuclear coordinates $\mathbf{R}$, but the equations of motions (for $\mathbf{R}$) are modified by the presence of the electric field.


## Overview 
#### Electric field
The electric field-driving will be performed on top of a ```NVE``` simulation.
In this case the adopted electric field is a **plane wave** with a **gaussian envelope function**.
The **frequency** $\omega$ and **direction** $\hat{\mathbf{E}}_{amp}$ of the plane wave (```115THz``` along the ```y```-axis) are chose to be resonant with one of the three vibrational mode of the system.
The **amplitude** $|\mathbf{E}_{amp}|$ of the field is chose to be not too intense.
The gaussian envelope function is determined by its **mean** $\mu$ (when the electric field will reach its maximum value) and **std** $\sigma$ (the duration of the pulse).
This is the analytic expression of the applied electric field is:
$$
\mathbf{E}\left(t\right) = \mathbf{E}_{amp}\cos \left( \omega t + \phi \right) \frac{1}{{\sigma \sqrt{2\pi}}} e^{-\frac{{(x - \mu)^2}}{{2\sigma^2}}}
$$

The parameters that define the shape of the electric field are defined in the ```ensemble``` field in the ```input.xml``` file. 
Some comments are provided in that file as well.

The parameters of the electric field are provided in the `input.xml` as follows (their names are self-explaining):
```xml
<motion mode='dynamics'> 
    <dynamics mode='eda-nve'>                  <!-- Electric Dipole Approximation (EDA) on top of a NVE simulation -->
        <timestep units='femtosecond'> 1 </timestep>
        <efield>
            <amp> [ 0 , 5e-2 , 0 ] </amp>             <!-- amplitude of the electric field pulse in a.u. -->
            <freq units="thz"> 115 </freq>            <!-- frequency of the electric field pulse in THz -->
            <peak units="picosecond"> 0.2 </peak>     <!-- peak of the electric field pulse -->
            <sigma units="picosecond"> 0.05 </sigma>  <!-- width (standard deviation of the gaussian) of the electric field pulse -->
        </efield>
        <bec mode="file"> bec.txt  </bec>
    </dynamics>
</motion>
```

#### Dipole and Born Effective Charge tensors
The dipole $\mathbf{d}$ of the system is not actually computed in this example.
In the same way we do not *actually* need the energy to integrate Hamilton's equations of motion but only the forces, we just need the Born Effective Charge tensors $Z^*_{ij}$ to run a light-driven simulation, i.e. the derivative of the dipole $\mathbf{d}$ w.r.t. nuclear coordinates $\mathbf{R}$:
$$
Z^*_{ij} = \frac{1}{q_e} \frac{\partial d_i}{\partial R_j}
$$

In this example we will provide the Born Effective Charge tensors $Z^*_{ij}$ using the input file ```bec.txt```, assuming (for semplicity) that their value is constant for each nuclear configuration $\mathbf{R}$.

There are two different ways to specify fixed Born Effective Charge tensors:
- by storing their values to file (as done in this example), and specify the filepath in `input.xml`
```xml
<motion mode='dynamics'> 
    <dynamics mode='eda-nve'>                  <!-- Electric Dipole Approximation (EDA) on top of a NVE simulation -->
        <timestep units='femtosecond'> 1 </timestep>
        <efield>
            <amp> [ 0 , 5e-2 , 0 ] </amp>             <!-- amplitude of the electric field pulse in a.u. -->
            <freq units="thz"> 115 </freq>            <!-- frequency of the electric field pulse in THz -->
            <peak units="picosecond"> 0.2 </peak>     <!-- peak of the electric field pulse -->
            <sigma units="picosecond"> 0.05 </sigma>  <!-- width (standard deviation of the gaussian) of the electric field pulse -->
        </efield>
        <bec mode="file"> bec.txt  </bec>
    </dynamics>
</motion>
```
`bec.txt` has to contain a matrix of shape $3N_a\times 3$, where the columns correspond to the direction of the dipole, and the rows to the coordinates w.r.t. which the derivative is computed. We are giving an explicit example in the following (neglicting the factor $1/q_e$ just for the clarity of exposition):
$$
(\text{as saved to file})
\quad
q_e Z_{ij} = 
\begin{bmatrix}
    \frac{\partial P_x}{\partial R^1_x} &  
    \frac{\partial P_y}{\partial R^1_x} &
    \frac{\partial P_z}{\partial R^1_x} \\ 
    \frac{\partial P_x}{\partial R^1_y} &  
    \frac{\partial P_y}{\partial R^1_y} &
    \frac{\partial P_z}{\partial R^1_y} \\
    \frac{\partial P_x}{\partial R^1_z} &  
    \frac{\partial P_y}{\partial R^1_z} &
    \frac{\partial P_z}{\partial R^1_z} \\
    \frac{\partial P_x}{\partial R^2_x} &  
    \frac{\partial P_y}{\partial R^2_x} &
    \frac{\partial P_z}{\partial R^2_x} \\ 
    \vdots & \vdots & \vdots \\
    \frac{\partial P_x}{\partial R^N_z} &  
    \frac{\partial P_y}{\partial R^N_z} &
    \frac{\partial P_z}{\partial R^N_z} \\    
\end{bmatrix}
$$
- by setting their values "manually" in `input.xml`:
```xml
<motion mode='dynamics'> 
    <dynamics mode='eda-nve'>                  <!-- Electric Dipole Approximation (EDA) on top of a NVE simulation -->
        <timestep units='femtosecond'> 1 </timestep>
        <efield>
            <amp> [ 0 , 5e-2 , 0 ] </amp>             <!-- amplitude of the electric field pulse in a.u. -->
            <freq units="thz"> 115 </freq>            <!-- frequency of the electric field pulse in THz -->
            <peak units="picosecond"> 0.2 </peak>     <!-- peak of the electric field pulse -->
            <sigma units="picosecond"> 0.05 </sigma>  <!-- width (standard deviation of the gaussian) of the electric field pulse -->
        </efield>
        <bec mode="file"> 
        [    0.7221015439,	-0.0000028201,	-0.0000000012,
            -0.0000028201,	-0.4070007918,	 0.0001345078,
             0.0000000000,	-0.0000003644,	-0.6096787984,
             0.3610505794,	 0.0000014101,	 0.0000003858,
             0.0000014101,	 0.2035008826,	-0.0431014831,
             0.0000006580,	-0.0735216582,	 0.3047408429,
             0.3610509645,	 0.0000014101,	-0.0000003846,
             0.0000014101,	 0.2034999092,	 0.0429669754,
            -0.0000006580,	 0.0735220226,	 0.3049379554 ]  
        </bec>
    </dynamics>
</motion>
```
where the array that has to be provided in `input.xml` is just the "flattened" version of the matrix contained in `bec.txt`.

#### Energy and forces
In this example the energy and forces of the isolated water molecule are provided by the ```pswater``` model, but you are free to use your favorite DFT code.
 
## How to run
Run ```i-PI``` as usual, and then call the driver (using a UNIX socket) for the energy and forces:
```bash
i-pi-driver -u -a driver -m pswater # > pswater.out &
```
In case you want to use an INET socket, modify ```input.xml``` as described by the comments therein, and run the driver with:
```bash
i-pi-driver -a ipaddress -p 27683 -m pswater # > pswater.out &
```
where ```ipaddress``` is your ```ip``` address and the port number (```-p```) should respect the criteria provided by the ```i-PI``` manual.

## Observations
Having a look at the output file of the properties ```i-pi.properties.out``` you will see that the systems stay in the initial configuration for the first of the dynamics.
But when the electric field starts to have non-negligible values, the kientic energy of the system increases.
If you have the pssibility to visualize the trajectory of the system with a software (`XCrySDen`, `ovito`, etc...) you will clearly see one (and only one) of the vibrational modes getting more and more excited as long as the electric field intensity increases.
(In the last step of the dynamics the electric field is nearly zero)
