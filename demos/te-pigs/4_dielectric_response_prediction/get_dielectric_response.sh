# Downloads the mu-alpha-models
wget -nc https://github.com/venkatkapil24/ML-quantum-vibrational-spectroscopy/raw/main/MACE_dipole_pol.model

# Uses the evaluate script to estimate mu and alpha
mace_eval_mu_alpha \
  --configs="../3_production_simulations/simulation.xc.xyz" \
  --model="MACE_dipole_pol.model" \
  --output="output.extxyz" \
  --device=cuda \
  --batch_size=10 \
#  --compute_dielectric_derivatives

