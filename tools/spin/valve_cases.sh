#!/bin/bash
# Builds the spin-valve test cases for the carbon chain and a run script.
# Usage:  bash valve_cases.sh [WORKDIR]      (default WORKDIR = /tmp/valve)
# Needs spinor_tools.py next to this file.
# MODE=mine    (default) branch with DeviceConfig.spinor_basis; 26 orbitals per atom, potential repeated per spin.
# MODE=devspin upstream dev-spin branch: spin = true, 13 orbitals per atom (count per spin channel), original potential.
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
OUT=${1:-/tmp/valve}
EX=$HOME/quatrex/examples/cp2k/carbon-chain
ORIG=$EX/inputs
CFG=$EX/qtbm/quatrex_config.toml
J=0.05
MODE=${MODE:-mine}
S=1300                      # first orbital of the right half: 100 atoms x 13 orbitals
tool="python $HERE/spinor_tools.py"
mkdir -p "$OUT"

make_config() {             # name, spinor(0/1)
  sed "s#^input_dir = .*#input_dir = \"$OUT/$1/inputs\"#; s#^simulation_dir = .*#simulation_dir = \"$OUT/$1\"#" "$CFG" > "$OUT/$1/quatrex_config.toml"
  if [ "$2" = 1 ]; then
    if [ "$MODE" = devspin ]; then
      sed -i 's/^\[device\]$/[device]\nspin = true/' "$OUT/$1/quatrex_config.toml"
    else
      sed -i 's/"C" = 13/"C" = 26/' "$OUT/$1/quatrex_config.toml"
      sed -i 's/^\[device\]$/[device]\nspinor_basis = true/' "$OUT/$1/quatrex_config.toml"
    fi
  fi
}
spinor_case() {             # name, direction of the right half (3 numbers)
  name=$1; shift
  rm -rf "$OUT/$name"; mkdir -p "$OUT/$name"; cp -rL "$ORIG" "$OUT/$name/inputs"
  $tool spinor "$ORIG/hamiltonian.h5" "$OUT/$name/inputs/hamiltonian.h5" --zeeman $J --direction 0 0 1 --split $S --direction-right "$@"
  $tool spinor "$ORIG/overlap.h5" "$OUT/$name/inputs/overlap.h5"
  if [ "$MODE" != devspin ]; then
    $tool potential "$ORIG/potential.npy" "$OUT/$name/inputs/potential.npy"
  fi
  make_config "$name" 1
}
spinless_case() {           # name, shift left, shift right
  name=$1
  rm -rf "$OUT/$name"; mkdir -p "$OUT/$name"; cp -rL "$ORIG" "$OUT/$name/inputs"
  $tool shift "$ORIG/hamiltonian.h5" "$OUT/$name/inputs/hamiltonian.h5" --shift "$2" --shift-right "$3" --split $S
  make_config "$name" 0
}

spinor_case P     0 0 1
spinor_case AP    0 0 -1
spinor_case X     1 0 0
spinor_case Y     0 1 0
spinor_case D45   0.70710678 0 0.70710678
spinor_case D135  0.70710678 0 -0.70710678
spinless_case upAP  $J -$J
spinless_case dnAP -$J  $J

cat > "$OUT/run_qtbm.py" <<'PY'
import sys
from pathlib import Path
from quatrex.cli.main import run as cli_run
cli_run(Path(sys.argv[1]), abort_on_exception=False)
print("DONE", sys.argv[1])
PY
cat > "$OUT/run_all.sh" <<RUN
#!/bin/bash
cd "$OUT"
for c in upAP dnAP P AP X Y D45 D135; do
  rm -f \$c.log
  python run_qtbm.py $OUT/\$c/quatrex_config.toml > \$c.log 2>&1
done
RUN
echo "cases built in $OUT"; grep -c spinor_basis "$OUT"/*/quatrex_config.toml
