#!/bin/bash
# Run all (or selected) pw.x jobs from inputs/, writing outputs/ + logs/.
#   bash run_all.sh                 # every inputs/*.in
#   bash run_all.sh bulk gas_co     # only the named jobs
# Skips any job whose output already contains JOB DONE.
set -e
cd "$(dirname "$0")"
source ../env.sh
# PBE+D3 is GPU-bound; keep OMP at 1 (32 threads oversubscribed the CPU and
# nearly doubled per-iteration time on the slab).
export OMP_NUM_THREADS=1
mkdir -p outputs logs outdir

if [ $# -gt 0 ]; then
    JOBS=("$@")
else
    JOBS=()
    for f in inputs/*.in; do JOBS+=("$(basename "${f%.in}")"); done
fi

echo "pw.x: $(which pw.x)"
echo "Jobs: ${JOBS[*]}"
echo "----"
for name in "${JOBS[@]}"; do
    IN="inputs/$name.in"; OUT="outputs/$name.out"; ERR="logs/$name.err"
    [ -f "$IN" ] || { echo "SKIP $name (no input)"; continue; }
    if [ -f "$OUT" ] && grep -q 'JOB DONE' "$OUT"; then
        echo "==> $name : already done, skipping"; continue
    fi
    echo "==> $name : started $(date '+%H:%M:%S')"; t0=$(date +%s)
    if pw.x -in "$IN" > "$OUT" 2> "$ERR"; then
        e=$(grep '^!' "$OUT" | tail -1 | awk '{print $5}')
        echo "    done in $(( $(date +%s) - t0 ))s   E = $e Ry"
    else
        echo "    FAILED (last 20 lines of $ERR):"; tail -20 "$ERR"; exit 1
    fi
done
echo "All requested jobs complete."
