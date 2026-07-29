#!/usr/bin/env bash
# Download MACE-POLAR-1-M directly into this example directory.

set -euo pipefail

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
model_path="${script_dir}/MACE-POLAR-1-M.model"
partial_path="${model_path}.part"
model_url="https://github.com/ACEsuit/mace-foundations/releases/download/mace_polar_1/MACE-POLAR-1-M.model"

if [[ -s "${model_path}" ]]; then
    echo "Model already present: ${model_path}"
    exit 0
fi

curl --fail --location --show-error \
    --output "${partial_path}" \
    "${model_url}"
mv -- "${partial_path}" "${model_path}"
echo "Downloaded ${model_path}"
