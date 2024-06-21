# Downloads the mu-alpha-models
wget -nc https://github.com/venkatkapil24/ML-quantum-vibrational-spectroscopy/raw/main/MACE_dipole_pol.model

# MACE currently doesn't support the extxyz format
ln -s ../3_production_simulations/simulation.pos_0.extxyz dataset.xyz

# Uses the evaluate script to estimate mu and alpha
mace_eval_mu_alpha \
  --configs="dataset.xyz" \
  --model="MACE_dipole_pol.model" \
  --device=cuda \
  --batch_size=50 \
  --output="output.extxyz" \
# Use this flag to also export the spatial derivates of the dipole and polarizability in the extxyz format.
#  --compute_dielectric_derivatives
