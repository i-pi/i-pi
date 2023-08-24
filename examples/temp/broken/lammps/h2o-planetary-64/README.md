Planetary model example
======================

Background
-----------

The planetary model was developed by Smith et al.[1,2] and developed further by Willatt
et al.[3] The model is based on a locally-harmonic approximation of the fluctuation of an
imaginary-time point (or planet) about the imaginary-time path centroid. The effective
frequency matrix that governs this harmonic motion depends on the centroid configuration
and is therefore time dependent. An important feature of the model is that the coupling
between the planets and centroids is one directional; the centroid motion is independent
of the planet motion, and the time-evolution is modified to ensure energy conservation.
This means that planet evolution can be treated as a post-processing step.

Originally the model was framed in the context of the Feynman-Kleinert approximation for
the centroid potential of mean force[1,2]. This approach is useful if the Gaussian
integrals associated with the effective frequency matrix can be calculated analytically.
In Ref. 3, we modified the model to use a different effective frequency matrix, given by a
centroid-constrained ring-polymer average of the Hessian, and TRPMD for the centroid
dynamics instead of Feynman-Kleinert CMD. The advantage of this approach is that the
centroid statistics are then exact. The planetary model is implemented in i-PI in this
modified form. The code was written by M. J. Willatt, R. Benson, and M. Ceriotti.


Usage
---------
In the xml input, one should include both <motion mode="dynamics"> and <motion
mode="planetary">. The former describes the TRPMD simulation, and the latter provides
details of the centroid-constrained frequency matrix sampling, e.g.

   <motion mode="multi">
    <motion mode="dynamics">
      <dynamics mode='nvt' splitting="baoab">
	<thermostat mode="pile_g">
          <tau units="femtosecond">100</tau>
          <pile_lambda> 0.5 </pile_lambda>
	</thermostat>
	<timestep units="femtosecond"> 0.50 </timestep>
      </dynamics>
    </motion>
    <motion mode="planetary">
      <planetary mode='md'>
        <thermostat mode="pile_g">
          <tau units="femtosecond">1e10</tau>
          <pile_lambda> 0.5 </pile_lambda>
        </thermostat>
        <timestep units="femtosecond"> 0.50 </timestep>
        <stride> 2 </stride>
        <nsamples> 128 </nsamples>
      </planetary>
    </motion>
  </motion>

In the above example, 256 frequency matrix samples are averaged for each centroid
configuration. The sampling is performed with constrained-centroid dynamics with
thermostatting of the non-centroid ring polymer modes. The frequency matrix is sampled
after every stride centroid position updates (once every femtosecond in the above
example). The output of the simulation is binary file called [prefix].omega2, which contains the
effective frequency matrices from the simulation in a compressed sparse column format. One
should also record the centroid positions and momenta, e.g.

<trajectory filename="xc" stride="2" format="xyz"> x_centroid </trajectory>
<trajectory filename="pc" stride="2" format="xyz"> p_centroid </trajectory>

With the [prefix].omega2 file and associated xyz files, one can calculate planetary model
time-correlation functions using the Python program i-pi-planetary. As described in
i-pi-planetary (tools/py/planetary.py), one must provide functions for the observables A
and B.

Here we provide an example to compute the dipole moment autocorrelation function of qTIP4P
water molecules. You should run 

i-pi-planetary PLANETARY 0.25 300 4 8 16 12345 ./estmod_dipole

The real part of the planetary model dipole moment autocorrelation function is output in
PLANETARY_pl_qTIP4P-cmumu-re.dat.

References
----------

1. http://www.doi.org/10.1063/1.4922887
2. http://www.doi.org/10.1063/1.4922888
3. http://www.doi.org/10.1063/1.5004808
