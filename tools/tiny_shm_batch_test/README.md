Tiny batched-SHM MPI reproduction based on:
`/scratch/c_ccmd/mb/test/ipi-develop/ipi-stable/socket_profiling3/pimd-32_nvt_ffbatch_shm`

Purpose:
- keep the same `ffsocket mode='shm' batch='true' matching='lock'` path
- keep the MPI driver launch path (`i-pi-py_mpidriver`)
- minimize system size and run length so communication/debug tests are cheap

Files:
- `input.xml`: tiny 2-bead NVT run using the dummy Python PES
- `init.xyz`: 2-atom initial structure
- `run.sh`: launches i-PI and the MPI driver with the same basic environment assumptions as the original profiling run

Usage:
```bash
cd tools/tiny_shm_batch_test
bash run.sh
```

Optional SHM debug logging:
```bash
cd tools/tiny_shm_batch_test
IPI_SHM_DEBUG=1 IPI_SHM_DEBUG_FILE=/tmp/ipi_shm_debug.log bash run.sh
cat /tmp/ipi_shm_debug.log
```
