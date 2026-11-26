#!/bin/bash
# Builds two test cases for the upstream dev-spin branch (device.spin = true):
#   A: doubled spinor inputs kron(H, I2) used as they are
#   B: spinless inputs expanded by Quatrex (expand_spinless_inputs = true)
# Usage:  bash devspin_cases.sh [WORKDIR]       (default /tmp/dv)
# Environment overrides: CFG, SP0, REF, POT
set -e
OUT=${1:-/tmp/dv}
CFG=${CFG:-$HOME/quatrex/examples/cp2k/carbon-chain/qtbm/quatrex_config.toml}
SP0=${SP0:-/tmp/nc2/sp0/inputs}        # spinor inputs from spinor_tools.py
REF=${REF:-/tmp/nc2/ref/inputs}        # spinless carbon-chain inputs
POT=${POT:-$REF/potential.npy}         # potential used in both cases

rm -rf "$OUT/A" "$OUT/B"
mkdir -p "$OUT/A" "$OUT/B"
cp -rL "$SP0" "$OUT/A/inputs"
cp -rL "$REF" "$OUT/B/inputs"
cp "$POT" "$OUT/A/inputs/potential.npy"
cp "$POT" "$OUT/B/inputs/potential.npy"

build() {   # name, extra line(s) after spin = true
  f="$OUT/$1/quatrex_config.toml"
  sed "s#^input_dir = .*#input_dir = \"$OUT/$1/inputs\"#; s#^simulation_dir = .*#simulation_dir = \"$OUT/$1\"#" "$CFG" > "$f"
  sed -i "s/^\[device\]\$/[device]\nspin = true$2/" "$f"
  python3 -c "
import sys, tomllib
d = tomllib.load(open(sys.argv[1], 'rb'))['device']
print(sys.argv[1], {k: d[k] for k in d if 'spin' in k}, d['num_orbitals_per_atom'])" "$f"
}
build A ""
build B "\nexpand_spinless_inputs = true"
echo "cases built in $OUT"
