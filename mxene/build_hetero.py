#!/usr/bin/env python3
"""
Route B (publication-grade): Mg2Ni(0001) // Ti3C2 (sqrt3 x sqrt3)R30 interface.

Per the referee spec: the monoclinic Mg2NiH4(001) facet does NOT lattice-match
Ti3C2 (<5% impossible); the clean route is the hexagonal metallic
Mg2Ni(0001) 1x1 (a=5.128 A, relaxed) on Ti3C2 sqrt3xsqrt3 (a_relaxed=5.317 A),
locking the shared in-plane lattice to a=5.128 A and straining the thin MXene
by -3.6%. Hydride energetics come from near-surface H (H-vacancy formation +
CI-NEB), not a bulk DH.

This script only BUILDS and INSPECTS geometry (no DFT). Run and eyeball the
atom counts, layer spacings, and strain before generating pw.x inputs.
"""
from pathlib import Path
import numpy as np
from ase.io import read
from ase.build import surface, make_supercell
from ase import Atoms

HERE = Path(__file__).resolve().parent
MGNICO = HERE.parent / "mg-nico" / "outputs"

A_LOCK = 5.128          # shared in-plane lattice (relaxed Mg2Ni a), Angstrom
VAC = 20.0              # vacuum, Angstrom


def load_bulk(name):
    """Relaxed bulk Atoms from a mg-nico pw.x vc-relax output."""
    atoms = read(MGNICO / f"{name}.out", format="espresso-out", index=-1)
    return atoms


def mg2ni_0001_slab(nlayers=7):
    bulk = load_bulk("mg2ni")
    a = np.linalg.norm(bulk.get_cell()[0])
    print(f"  bulk Mg2Ni: {len(bulk)} atoms, a={a:.4f}, "
          f"c={np.linalg.norm(bulk.get_cell()[2]):.4f}")
    slab = surface(bulk, (0, 0, 1), nlayers, vacuum=VAC / 2)
    slab.center(vacuum=VAC / 2, axis=2)
    return slab


def ti3c2_sqrt3_slab():
    """Relaxed Ti3C2 1x1 -> sqrt3xsqrt3 R30 supercell, strained to a=A_LOCK."""
    ti = read(MGNICO.parent.parent / "mxene" / "outputs" / "ti3c2.out",
              format="espresso-out", index=-1)
    P = np.array([[2, 1, 0], [-1, 1, 0], [0, 0, 1]])   # sqrt3 x sqrt3 R30
    sc = make_supercell(ti, P)
    a_now = np.linalg.norm(sc.get_cell()[0])
    scale = A_LOCK / a_now
    cell = sc.get_cell()
    cell[0] *= scale
    cell[1] *= scale
    sc.set_cell(cell, scale_atoms=True)   # strain in-plane, keep z
    print(f"  Ti3C2 sqrt3 slab: {len(sc)} atoms, a {a_now:.4f} -> {A_LOCK:.4f} "
          f"(strain {100*(scale-1):+.2f}%)")
    return sc


def count_layers(atoms, tol=0.6):
    z = np.sort(atoms.get_positions()[:, 2])
    levels = [z[0]]
    for zi in z[1:]:
        if zi - levels[-1] > tol:
            levels.append(zi)
    return len(levels), levels


def stack_hetero(registry=(0.0, 0.0), gap=2.4, vac=20.0):
    """Mg2Ni(0001) slab (relaxed) + Ti3C2 sqrt3 slab (relaxed) stacked along z.
    Both share the in-plane metric (a=5.128, gamma=120); Ti3C2 is mapped into
    the Mg2Ni in-plane cell frame (fractional-preserving) then shifted by
    `registry` (fractional a,b) and lifted `gap` above the Mg2Ni top."""
    mg = read("outputs/mg2ni_slab1.out", format="espresso-out", index=-1)
    ti = read("outputs/ti3c2_sqrt3.out", format="espresso-out", index=-1)
    # map Ti3C2 into Mg2Ni in-plane cell frame (metrics match -> keeps fractional)
    newcell = ti.get_cell().copy()
    newcell[0] = mg.get_cell()[0]
    newcell[1] = mg.get_cell()[1]
    ti.set_cell(newcell, scale_atoms=True)
    mgz = mg.get_positions()[:, 2]
    tp = ti.get_positions().copy()
    shift = registry[0] * mg.get_cell()[0] + registry[1] * mg.get_cell()[1]
    tp[:, 0] += shift[0]; tp[:, 1] += shift[1]
    tp[:, 2] += (mgz.max() + gap) - ti.get_positions()[:, 2].min()
    hetero = mg.copy()
    hetero += Atoms(symbols=ti.get_chemical_symbols(), positions=tp)
    p = hetero.get_positions()
    zmin, zmax = p[:, 2].min(), p[:, 2].max()
    p[:, 2] += 2.0 - zmin
    hetero.set_positions(p)
    cell = hetero.get_cell()
    cell[2] = [0.0, 0.0, (zmax - zmin) + vac]
    hetero.set_cell(cell)
    hetero.wrap()
    return hetero


def min_dist(atoms):
    from ase.geometry import get_distances
    p = atoms.get_positions()
    d = get_distances(p, p, cell=atoms.get_cell(), pbc=[True, True, False])[1]
    np.fill_diagonal(d, 9e9)
    return d.min()


def freeze_z(atoms, frac=0.5):
    """z below which atoms are frozen: bottom `frac` of the slab thickness."""
    z = atoms.get_positions()[:, 2]
    return z.min() + frac * (z.max() - z.min())


def gen_inputs():
    """Write the stage-1/2 clean-slab pw.x inputs (no DFT here)."""
    from build_interface import write_qe
    MAGN = {"Ti": 0.4, "C": 0.0, "Ni": 0.3, "Mg": 0.0, "H": 0.0}

    ti = ti3c2_sqrt3_slab()
    write_qe("inputs/ti3c2_sqrt3.in", ti, prefix="ti3c2_sqrt3", nspin=2,
             magn=MAGN, kpts=(6, 6, 1), degauss=0.02, dipole=False,
             max_seconds=14400)

    for rep, nm in ((1, "mg2ni_slab1"), (2, "mg2ni_slab2")):
        s = mg2ni_0001_slab(rep)
        write_qe(f"inputs/{nm}.in", s, prefix=nm, nspin=2, magn=MAGN,
                 kpts=(6, 6, 1), degauss=0.02, dipole=False,
                 fixed_below=freeze_z(s, 0.5), max_seconds=28800)


if __name__ == "__main__":
    from collections import Counter
    import sys
    if "--gen" in sys.argv:
        gen_inputs()
        sys.exit(0)
    print("== Mg2Ni(0001) slabs (nlayers = unit-cell repeats) ==")
    for n in (1, 2, 3):
        s = mg2ni_0001_slab(n)
        nlay, _ = count_layers(s)
        z = s.get_positions()[:, 2]
        ca = np.linalg.norm(s.get_cell()[0])
        cb = np.linalg.norm(s.get_cell()[1])
        ang = np.degrees(np.arccos(np.dot(s.get_cell()[0], s.get_cell()[1]) /
                                   (ca * cb)))
        print(f"  rep={n}: {len(s)} atoms {dict(Counter(s.get_chemical_symbols()))}, "
              f"{nlay} atomic layers, slab-thick {z.max()-z.min():.2f} A")
        print(f"          in-plane a={ca:.3f} b={cb:.3f} gamma={ang:.1f} deg, "
              f"cell c={s.get_cell()[2][2]:.2f}")
    print("== Ti3C2 sqrt3 slab ==")
    ti = ti3c2_sqrt3_slab()
    ca = np.linalg.norm(ti.get_cell()[0]); cb = np.linalg.norm(ti.get_cell()[1])
    ang = np.degrees(np.arccos(np.dot(ti.get_cell()[0], ti.get_cell()[1]) / (ca*cb)))
    print(f"  composition: {dict(Counter(ti.get_chemical_symbols()))}, "
          f"in-plane a={ca:.3f} b={cb:.3f} gamma={ang:.1f} deg")
