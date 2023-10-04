NEB example
-------------

Once the env.sh file in the i-pi root directory has been sourced with

```bash
source <i-pi-root>/env.sh
```

Run i-pi with:

```bash
i-pi input.xml
```

And the driver:

```bash
i-pi-driver -u -h driver -p 20614 -m zundel
```

You can analyze the results using provided gnuplot scripts `plot*.plt`:
((Tested with gnuplot 5.4 patchlevel 1)
```bash
gnuplot plot-converg.plt
gnuplot plot-bells.plt
```

You can extract certain path as a single .xyz file using `tools/bash/get-certain-path.sh`.

Important Note:
   The code randomly crashes here with some versions of Numpy based on OpenBLAS.
   If this happens, use Numpy based on MKL, e.g. from Anaconda.
   YL (23Aug23) Tested numpy versions 1.21.5, 1.24.2

 
