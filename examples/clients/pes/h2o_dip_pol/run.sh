#!/usr/bin/env bash
set -euo pipefail

example_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
repo_root=$(cd "$example_dir/../../../.." && pwd)
build_dir=$(mktemp -d)
library="$build_dir/libh2o_dip_pol.so"
trap 'rm -rf "$build_dir"' EXIT

gfortran -O3 -ffree-line-length-none -fPIC -shared \
  "$repo_root/drivers/f90/pes/h2o_dip_pol.f90" -o "$library"

cd "$example_dir"
export IPI_H2O_DIP_POL_LIBRARY="$library"
export PYTHONPATH="$repo_root${PYTHONPATH:+:$PYTHONPATH}"

i-pi input.xml &
sleep 2
i-pi-driver -m qtip4pf -u -a h2o-dip-pol-qtip4pf
