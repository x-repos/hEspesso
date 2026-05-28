"""
Variant A: Fe substituted on Mg site (instead of Ni site).

Rationale: in the original R3 cells, Fe replaced Ni; both are 3d TMs of similar
size and chemistry, so the electronic perturbation was small (shift ≈ 0). Putting
Fe on a Mg site instead places a 3d TM in a position previously occupied by a
2s²2p⁶3s² (alkaline-earth) atom, far closer to the H sublattice in the hydride.
This usually gives a much larger destabilising effect on M-H bonds.

Cells generated:
  mg2ni_NbFeMg     18 atoms  (10 Mg + 6 Ni + 1 Nb + 1 Fe)  - 1 Mg→Nb + 1 Mg→Fe
  mg2nih4_NbFeMg   28 atoms  ( 6 Mg + 4 Ni + 1 Nb + 1 Fe + 16 H)
"""

from pathlib import Path
import numpy as np
from ase.io import read
from ase.io.espresso import write_espresso_in

HERE = Path(__file__).resolve().parent
INDIR = HERE / "inputs"
INDIR.mkdir(exist_ok=True)
CIFDIR = HERE.parent / "mgh2-cif"

PSEUDO_DIR = "../pseudo/"
PSEUDOS = {
    "H":  "H.pbe-rrkjus_psl.1.0.0.UPF",
    "Mg": "Mg.pbe-n-kjpaw_psl.0.3.0.UPF",
    "Ni": "ni_pbe_v1.4.uspp.F.UPF",
    "Nb": "Nb.pbe-spn-kjpaw_psl.0.3.0.UPF",
    "Fe": "Fe.pbe-spn-kjpaw_psl.0.2.1.UPF",
}


def base_input(prefix, calc="relax", nspin=2):
    inp = {
        "control": {
            "calculation": calc, "prefix": prefix,
            "pseudo_dir": PSEUDO_DIR, "outdir": "./outdir/",
            "tprnfor": True, "tstress": True,
            "restart_mode": "from_scratch",
            "etot_conv_thr": 1.4e-4, "forc_conv_thr": 1.0e-4,
        },
        "system": {
            "ecutwfc": 60.0, "ecutrho": 480.0,
            "occupations": "smearing", "smearing": "cold", "degauss": 0.01,
            "nspin": nspin,
        },
        "electrons": {
            "conv_thr": 1.2e-9, "mixing_beta": 0.4,
            "electron_maxstep": 120,
        },
    }
    if calc in ("relax", "vc-relax"):
        inp["ions"] = {"ion_dynamics": "bfgs"}
    return inp


def first_index(atoms, symbol, skip=0):
    """Return the (skip+1)-th index of `symbol`."""
    seen = 0
    for i, a in enumerate(atoms):
        if a.symbol == symbol:
            if seen == skip:
                return i
            seen += 1
    raise ValueError(f"Not enough {symbol} atoms in cell")


def substitute(atoms, idx, new_symbol):
    out = atoms.copy()
    out.symbols[idx] = new_symbol
    return out


def write_qe(path, atoms, prefix, kpts):
    inp = base_input(prefix, calc="relax", nspin=2)
    moms = np.zeros(len(atoms))
    for i, a in enumerate(atoms):
        if a.symbol == "Fe":
            moms[i] = 2.0
    atoms = atoms.copy()
    atoms.set_initial_magnetic_moments(moms)
    pseudos = {sym: PSEUDOS[sym] for sym in set(atoms.get_chemical_symbols())}
    with open(path, "w") as f:
        write_espresso_in(f, atoms, input_data=inp, pseudopotentials=pseudos, kpts=kpts)
    print(f"  wrote {path.name}: {len(atoms)} atoms, kpts={kpts}")


def main():
    mg2ni      = read(CIFDIR / "mg2ni.pwi")
    mg2nih4_28 = read(CIFDIR / "mg2nih4-28.pwi")

    # --- Variant A: Nb on first Mg, Fe on next-available Mg ---
    # Mg2Ni 18a: 12 Mg → 1 Nb + 1 Fe + 10 Mg
    mg2ni_NbFeMg = substitute(mg2ni, first_index(mg2ni, "Mg", skip=0), "Nb")
    mg2ni_NbFeMg = substitute(mg2ni_NbFeMg, first_index(mg2ni_NbFeMg, "Mg", skip=0), "Fe")

    # Mg2NiH4 28a: 8 Mg → 1 Nb + 1 Fe + 6 Mg
    mg2nih4_NbFeMg = substitute(mg2nih4_28, first_index(mg2nih4_28, "Mg", skip=0), "Nb")
    mg2nih4_NbFeMg = substitute(mg2nih4_NbFeMg, first_index(mg2nih4_NbFeMg, "Mg", skip=0), "Fe")

    write_qe(INDIR / "mg2ni_NbFeMg.in",   mg2ni_NbFeMg,   "mg2ni_NbFeMg",   kpts=(8, 8, 4))
    write_qe(INDIR / "mg2nih4_NbFeMg.in", mg2nih4_NbFeMg, "mg2nih4_NbFeMg", kpts=(8, 8, 6))

    print("\nDoping summary (Variant A: Fe on Mg site):")
    from collections import Counter
    for name, atoms in [("mg2ni_NbFeMg", mg2ni_NbFeMg), ("mg2nih4_NbFeMg", mg2nih4_NbFeMg)]:
        print(f"  {name:18s}: {dict(Counter(atoms.get_chemical_symbols()))}")


if __name__ == "__main__":
    main()
