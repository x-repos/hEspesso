"""Gate the bulk structure BEFORE any GPU time: confirm structures/V2O5_alpha.cif
is α-V2O5 (orthorhombic, layered Pmmn topology) via composition, the lattice
signature (short ~3.56 A chain axis, long ~11.5 A axis), and V-O bonding (a
short vanadyl bond + ~5-fold V coordination). No spglib dependency.
"""
import sys
from pathlib import Path
import numpy as np
from ase.io import read

HERE = Path(__file__).resolve().parent


def main():
    a = read(HERE / "structures" / "V2O5_alpha.cif")
    e = []
    s = a.get_chemical_symbols()
    if (s.count("V"), s.count("O")) != (4, 10):
        e.append(f"composition {s.count('V')}V{s.count('O')}O, want V4O10")
    L = sorted(a.cell.lengths())
    ang = a.cell.angles()
    if not all(abs(x - 90) < 3 for x in ang):
        e.append(f"not orthorhombic: angles {[round(x, 1) for x in ang]}")
    if not (3.3 < L[0] < 3.9):
        e.append(f"short axis {L[0]:.2f} A, want ~3.56 (chain b)")
    if not (10.8 < L[2] < 12.3):
        e.append(f"long axis {L[2]:.2f} A, want ~11.5 (a)")
    d = a.get_all_distances(mic=True)
    n = len(a)
    if any(d[i][j] < 0.9 for i in range(n) for j in range(i + 1, n)):
        e.append("atom overlap < 0.9 A")
    O = [i for i, x in enumerate(s) if x == "O"]
    for i, x in enumerate(s):
        if x != "V":
            continue
        vo = sorted(d[i][j] for j in O)
        if not (1.50 < vo[0] < 1.75):
            e.append(f"V#{i} shortest V-O {vo[0]:.2f} A, want vanadyl ~1.58")
        nb = sum(v < 2.3 for v in vo)
        if not (4 <= nb <= 6):
            e.append(f"V#{i} coordination {nb}, want ~5")
    if e:
        print("STRUCTURE REJECTED:")
        for x in e:
            print("  -", x)
        sys.exit(1)
    print(f"OK: alpha-V2O5  V4O10, abc={[round(x, 3) for x in L]}, "
          f"angles~90, vanadyl bond present.")


if __name__ == "__main__":
    main()
