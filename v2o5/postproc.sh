#!/bin/bash
# Post-process the bare slab + each gas's winning adsorption complex with the
# CPU QE build (espresso_cpu). That build is an HPC-X OpenMPI binary, so it must
# be launched through mpirun (a bare ./pp.x hangs in MPI_Init). For each prefix:
#   * valence charge density cube  (plot_num=0)  -> Bader charge transfer + CDD
#   * electrostatic potential      (plot_num=11) -> average.x planar avg -> work fn
#   * projwfc.x                                   -> PDOS + Lowdin charges
#
# Usage:  bash postproc.sh slab ads_co_s0 ads_nh3_s0 ...
set -e
cd "$(dirname "$0")"
MPIRUN=$(ls /opt/nvidia/hpc_sdk/Linux_x86_64/*/comm_libs/*/hpcx/*/ompi/bin/mpirun 2>/dev/null | head -1)
QECPU=/home/x/Programs/espresso_cpu/bin
NP=${NP:-4}
RUN="$MPIRUN -np $NP"
BADER=$(command -v bader 2>/dev/null || echo ./bin/bader)
mkdir -p pp logs

for pre in "$@"; do
    [ -d "outdir/${pre}.save" ] || { echo "SKIP $pre (no outdir/${pre}.save)"; continue; }
    echo "== $pre =="

    # 1) valence charge density (Bader + CDD)
    cat > pp/${pre}_rho.in <<EOF
&INPUTPP
  prefix='${pre}', outdir='./outdir/', plot_num=0, filplot='pp/${pre}.rho'
/
&PLOT
  iflag=3, output_format=6, fileout='pp/${pre}.cube'
/
EOF
    $RUN $QECPU/pp.x -in pp/${pre}_rho.in > logs/${pre}_pp_rho.out 2>&1
    if [ -x "$BADER" ]; then
        "$BADER" pp/${pre}.cube > logs/${pre}_bader.log 2>&1 && mv -f ACF.dat pp/${pre}_ACF.dat || echo "  bader failed (Lowdin fallback in analyze)"
        rm -f AVF.dat BCF.dat
    fi

    # 2) electrostatic potential -> planar average -> work function
    cat > pp/${pre}_pot.in <<EOF
&INPUTPP
  prefix='${pre}', outdir='./outdir/', plot_num=11, filplot='pp/${pre}.pot'
/
&PLOT
  iflag=3, output_format=6, fileout='pp/${pre}.potcube'
/
EOF
    $RUN $QECPU/pp.x -in pp/${pre}_pot.in > logs/${pre}_pp_pot.out 2>&1
    printf "1\npp/%s.pot\n1.0\n2000\n3\n3.0\n" "$pre" | $RUN $QECPU/average.x > pp/avg_${pre}.dat 2> logs/${pre}_avg.out

    # 3) PDOS + Lowdin charges
    cat > pp/${pre}_pdos.in <<EOF
&PROJWFC
  prefix='${pre}', outdir='./outdir/', filpdos='pp/${pre}', DeltaE=0.05, ngauss=0, degauss=0.01
/
EOF
    $RUN $QECPU/projwfc.x -in pp/${pre}_pdos.in > logs/${pre}_pdos.out 2>&1
    echo "  done $pre"
done
echo "post-processing complete."
