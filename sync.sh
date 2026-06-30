#!/usr/bin/env bash
set -euo pipefail


# SYNC MIO — hEspesso / V2O5 Quantum ESPRESSO study  <->  cluster
REMOTE=mio
REMOTE_DIR="~/Workspace/3-hEspesso"
LOCAL_DIR="$HOME/Workspace/3-hEspesso"




# === push the study UP to the cluster (to run on SLURM) ======================
# Ships the v2o5 study (scripts, structures, SLURM submit, generated inputs)
# plus the full pseudo/ library (~160 MB one-time; rsync -z skips unchanged on
# re-runs).  --relative keeps repo-relative paths so files land in the matching
# remote dir (the `/./` marks where the relative path starts).
#
# EXCLUDED so we never ship GPU-machine cruft or clobber cluster results:
#   outdir/        QE scratch (huge; the cluster makes its own)
#   outputs/ logs/ the GPU run's results -- a DIFFERENT config (1-layer/PBE+D3);
#                  shipping them would make the cluster SKIP jobs that already
#                  show "JOB DONE" and analyze stale numbers
#   docs/          the big reference PDFs the cluster doesn't need
#   __pycache__ / *.pyc
rsync -avz --progress --human-readable --relative \
    --exclude='__pycache__/' \
    --exclude='*.pyc' \
    --exclude='outdir/' \
    --exclude='outputs/' \
    --exclude='logs/' \
    --exclude='pp/' \
    --exclude='docs/' \
    "$LOCAL_DIR/./v2o5/" \
    "$LOCAL_DIR/./pseudo/" \
    "$REMOTE:$REMOTE_DIR/"


# === pull the cluster OUTPUTS back DOWN (after the SLURM run) =================
# Uncomment when the run is done.  Pulls only the QE outputs + post-processing
# (NOT the huge outdir/ scratch) so you can run analyze.py on this machine.
# rsync -avz --progress --human-readable --relative \
#     --exclude='__pycache__/' \
#     "$REMOTE:$REMOTE_DIR/./v2o5/outputs/" \
#     "$REMOTE:$REMOTE_DIR/./v2o5/pp/" \
#     "$LOCAL_DIR/"
