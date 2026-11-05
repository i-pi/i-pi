#!/usr/bin/env bash
set -euo pipefail
script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
model_path="${script_dir}/MACE-POLAR-1-M.model"
curl --fail --location --show-error --output "${model_path}.part" \
  "https://github.com/ACEsuit/mace-foundations/releases/download/mace_polar_1/MACE-POLAR-1-M.model"
mv -- "${model_path}.part" "${model_path}"
