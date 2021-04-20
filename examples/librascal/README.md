Run a simulation of the Zundel cation using a ML potential based on 
SOAP features and Gaussian processes/KRR. The potential is obtained
using https://github.com/cosmo-epfl/librascal for fitting, and a 
dedicated interface to compute energy and forces, still using librascal.
The potential can be trained using an example in the librascal source
tree, or fetched as

wget https://www.dropbox.com/s/rk4pddt5uctn0gy/zundel_model.json?dl=0

After having obtained the potential file, the example can be run as usual

$ i-pi input.xml &
$ i-pi-py_driver -m rascal -a zundel -u -o zundel_model.json,h5o2+.extxyz

