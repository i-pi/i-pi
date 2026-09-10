#!/usr/bin/env bash

set -euo pipefail

example_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
cd "${example_dir}"
assets_dir="../mace-polar+E-client"

"${assets_dir}/download_mace_polar.sh"

i-pi input.xml > log.i-pi 2>&1 &
ipi_pid=$!
sleep 1

i-pi-py_driver -u -a mace-polar-static-extra -m extmace -o \
  "template=${assets_dir}/start.extxyz,model=${assets_dir}/MACE-POLAR-1-M.model,device=cpu,mace_kwargs=${assets_dir}/mace_kwargs.json,requires_extra=true" \
  > log.driver 2>&1 &
driver_pid=$!

wait "${ipi_pid}"
wait "${driver_pid}"
