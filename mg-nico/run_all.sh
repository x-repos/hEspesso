#!/bin/bash
# Run all (or selected) pw.x jobs from inputs/, writing to outputs/ + logs/.
#
# Usage:
#   bash run_all.sh                 # run all 13 jobs in default order
#   bash run_all.sh mg h2 mgh2      # only the cheap baseline
#   bash run_all.sh mgni mgh2ni     # only the Mg:Ni pathway

set -e
cd "$(dirname "$0")"

# GPU pw.x environment (env.sh lives at the workspace root, one level up)
source ../env.sh

mkdir -p outputs logs outdir

# Default order: cheap first so failures are caught early
DEFAULT=(mg h2 mgh2 mgni mgco mg2ni mg2ni_ni mg2ni_co mg2ni_co_nisite
         mg2nih4 mg2nih4_ni mg2nih4_co mg2nih4_co_nisite mgh2ni mgh2co
         mgh2co_spinchk mg2ni_co_spinchk mg2ni_co_nisite_spinchk
         mg2nih4_ni_spinchk mg2nih4_co_spinchk mg2nih4_co_nisite_spinchk)
JOBS=("$@")
if [ ${#JOBS[@]} -eq 0 ]; then
    JOBS=("${DEFAULT[@]}")
fi

GPU=$(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | head -1)
echo "GPU : ${GPU:-none}"
echo "pw.x: $(which pw.x)"
echo "Jobs: ${JOBS[*]}"
echo "----"

for name in "${JOBS[@]}"; do
    IN="inputs/$name.in"
    OUT="outputs/$name.out"
    ERR="logs/$name.err"
    if [ ! -f "$IN" ]; then
        echo "SKIP $name : $IN not found"
        continue
    fi
    if [ -f "$OUT" ] && grep -q 'JOB DONE' "$OUT"; then
        e=$(grep '^!' "$OUT" | tail -1 | awk '{print $5}')
        echo "==> $name : already converged (E = $e Ry), skipping"
        continue
    fi

    echo
    echo "==> $name : started $(date '+%H:%M:%S')"
    t0=$(date +%s)
    if pw.x -in "$IN" > "$OUT" 2> "$ERR"; then
        t1=$(date +%s)
        e=$(grep '^!' "$OUT" | tail -1 | awk '{print $5}')
        echo "    done in $((t1-t0))s   E = $e Ry"
    else
        echo "    FAILED (last 20 lines of $ERR):"
        tail -20 "$ERR"
        exit 1
    fi
done

echo
echo "All requested jobs complete."
