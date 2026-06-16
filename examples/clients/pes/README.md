Examples for the PES modules
============================

These examples show how to run a simulation using one of the PES defined in
`ipi.pes`.

harmonic/ shows how to run using the Python driver `i-pi-py_driver`
harmonic-direct/ shows how to run using a `<ffdirect>` clause in the input
harmonic-batched/ shows how to run using `<ffdirect>` in a way that sends a batch of structures at the same time
harmonic-batched-socket/ shows how to send a batch of structures per request over a socket (`i-pi-py_driver`), by setting `<batch_size>` on `<ffsocket>`
gas-batched-socket/ the same batched-socket path with the vectorized "gas" (ideal-gas, no-op) PES, a minimal-cost target to gauge the communication overhead

In the batched-socket examples the `i-pi-py_driver` command is exactly the same as in the non-batched case: i-PI announces the batch size to the driver through the INIT string, and a batch-aware driver (e.g. `harmonic`, `gas`) switches to batched evaluation automatically.

custom/ and custom-direct/ show how to run using a "custom" PES file that is not part of those distributed with i-PI
