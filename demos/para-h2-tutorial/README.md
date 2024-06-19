An introductory tutorial based on the p-H2 driver
=================================================

Authors: `Joshua More`, `Michele Ceriotti <michele.ceriotti@gmail.com>`


This gives an example of simulations para-H2 with the isotropic Silvera-Goldman pair
potential, as discussed in the tutorial section of the user manual. The examples include
classical NVT and NPT simulations, as well as imaginary-time path integral MD. 
Please refer to [this section](https://ipi-code.org/i-pi/tutorials.html)
for more detailed instructions on how to run this example.

To run tutorial 1:

```bash
source <i-pi-root>/env.sh
i-pi tutorial-1.xml
i-pi-driver -m sg -a localhost -o 15 -p 31415
```
