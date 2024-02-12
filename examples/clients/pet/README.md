PET - i-PI example
========================

Run a simulation of the Zundel cation using a ML potential based on 
the Point Edge Equivariant model (cf. https://arxiv.org/abs/2305.19302)
trained on a large dataset of CH4 structues, including many high-energy 
configurations (cf. http://doi.org/10.24435/materialscloud:qy-dp).

This is not a highly-optimized model and is only meant for quick testing
and demonstrative purposes. The potential files can be fetched from zenodo 

```bash
wget https://zenodo.org/records/10250171/files/methane.pet_small.zip
unzip methane.pet_small.zip
```

A version of the PET code that is compatible with this example can be 
obtained with 

```bash
pip install git+https://github.com/serfg/pet@83454c6
```

After having obtained the potential file, and a compatible version of PET,
the example can be run as usual

```bash
i-pi input.xml &> log.ipi &
i-pi-py_driver -m pet -o methane.pet_small,ch4.xyz -a pet -u &> log.pet 
```
