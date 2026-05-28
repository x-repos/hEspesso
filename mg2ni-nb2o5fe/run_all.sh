#!/bin/bash
# Run all 7 catalyst-DFT calculations sequentially on the GPU.
# Usage: bash run_all.sh             # run all
#        bash run_all.sh h2 mg2ni    # run a subset
set -e

cd "$(dirname "$0")"

# Load NVHPC + GPU pw.x. env.sh lives at the workspace root
# (../../../env.sh from this script: ../=espresso/, ../../=Mg2Ni-Catalyst-Repo/, ../../../=espresso/ workspace root).
source ../../../env.sh

# Layout (relative to this script):
#   inputs/   QE input files (*.in)        — read-only
#   outputs/  SCF stdout logs (*.out)      — final energies + diagnostics
#   outdir/   QE scratch (wavefns, .save/) — overwritten per run, matches outdir= in *.in
#   logs/     stderr from pw.x             — for debugging failures
mkdir -p outputs outdir logs

ALL_INPUTS=(mg2ni_Nb mg2nih4_Nb mg2ni_NbFe mg2nih4_NbFe mgo nb2o5 ni mg2ni_NbFeMg mg2nih4_NbFeMg)
INPUTS=("$@")
if [ ${#INPUTS[@]} -eq 0 ]; then
    INPUTS=("${ALL_INPUTS[@]}")
fi

echo "GPU: $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | head -1)"
echo "pw.x: $(which pw.x)"
echo "Running: ${INPUTS[@]}"
echo "----"

for name in "${INPUTS[@]}"; do
    if [ ! -f "inputs/$name.in" ]; then
        echo "SKIP $name: inputs/$name.in not found"
        continue
    fi
    echo ""
    echo "==> $name (started $(date '+%H:%M:%S'))"
    t0=$(date +%s)

    pw.x -in "inputs/$name.in" > "outputs/$name.out" 2> "logs/$name.err" || {
        echo "FAILED: $name (see logs/$name.err)"
        tail -20 "logs/$name.err"
        exit 1
    }

    t1=$(date +%s)
    elapsed=$((t1 - t0))
    energy=$(grep '^!' "outputs/$name.out" | tail -1 | awk '{print $5}')
    echo "    done in ${elapsed}s   final E = ${energy} Ry"
done

echo ""
echo "All requested calculations complete."
