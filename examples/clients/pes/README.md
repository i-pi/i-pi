Examples for the PES modules
============================

These examples show how to run a simulation using one of the PES defined in
`ipi.pes`, evaluated in-process through a `<ffdirect>` clause (no socket).

harmonic-direct/ shows how to run using a `<ffdirect>` clause in the input
harmonic-batched/ shows how to run using `<ffdirect>` in a way that sends a batch of structures at the same time

custom/ and custom-direct/ show how to run using a "custom" PES file that is not part of those distributed with i-PI

For running the same PES drivers over a socket or MPI (unix, inet, shm, ffmpi),
see `../communication/`.
