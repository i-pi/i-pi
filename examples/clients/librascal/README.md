librascal - i-PI example
========================

Run a simulation of the Zundel cation using a ML potential based on 
SOAP features and Gaussian processes/KRR. The potential is obtained
using [librascal](https://github.com/cosmo-epfl/librascal) for fitting, and  
a dedicated interface to compute energy and forces, still using librascal.
The potential can be trained using an example in the librascal source
tree, or fetched as

```bash
wget https://zenodo.org/record/7927211/files/zundel_model.json
```

After having obtained the potential file, the example can be run as usual

```bash
i-pi input.xml &
i-pi-py_driver -m rascal -a zundel -u -o zundel_model.json,h5o2+.extxyz
```
