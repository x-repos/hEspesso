# GPU Quantum ESPRESSO environment (QE 7.5, NVHPC 25.9 build, RTX 5090).
#
# The NVHPC compilers/lib dir must come FIRST in LD_LIBRARY_PATH so that
# libgomp.so.1 resolves to NVHPC's GOMP-compatibility shim instead of GNU
# libgomp: pw.x links both libnvomp and (via the system libfftw3_omp) a
# libgomp.so.1, and mixing the two OpenMP runtimes aborts with
# "libgomp: TODO" at startup.
export PATH="/home/x/Programs/espresso_gpu/bin:$PATH"
export LD_LIBRARY_PATH="/opt/nvidia/hpc_sdk/Linux_x86_64/25.9/compilers/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export NVCOMPILER_OMP_DISABLE_WARNINGS=true
