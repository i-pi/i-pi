#!/bin/bash

set -e

# downloads PET-MAD-xs v1.6.0 from Hugging Face, and converts it to a pet-jax
# checkpoint in pet-mad-xs/ (needs `pip install "pet-jax[convert]"`)
petjax-convert pet-mad-xs --out pet-mad-xs
