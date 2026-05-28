"""
Generate the doped QE input files using mgh2-cif's relaxed structures as source.

Pristine energies are reused directly from /home/x/Workspace/espresso/mgh2-cif/:
  h2.pwo            -> -2.33322115 Ry  (gas)
  mg.pwo            -> -67.23391048  Ry (hcp Mg)
  mg2ni.pwo         -> -2463.30740842 Ry (18-atom hP18, 6 f.u.)
  mg2nih4-28.pwo    -> -1661.22169488 Ry (28-atom monoclinic, 4 f.u., gives ΔH ≈ -57)
  mgh2.pwo          -> -71.98196008  Ry (rutile)

Only the 4 doped cells + 2 auxiliary references (MgO, Nb2O5) are generated here.
Settings match mgh2-cif: 60/480 Ry, 8x8x4 / 8x8x6 k-points, mixing_beta=0.4.
"""

from pathlib import Path
import numpy as np
from ase import Atoms
from ase.io import read
from ase.io.espresso import write_espresso_in
from ase.spacegroup import crystal

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
    "O":  "O.pbe-n-kjpaw_psl.0.1.UPF",
}


def base_input(prefix, calc="scf", nspin=1):
    inp = {
        "control": {
            "calculation": calc,
            "prefix": prefix,
            "pseudo_dir": PSEUDO_DIR,
            "outdir": "./outdir/",
            "tprnfor": True,
            "tstress": True,
            "restart_mode": "from_scratch",
            "etot_conv_thr": 1.4e-4,
            "forc_conv_thr": 1.0e-4,
        },
        "system": {
            "ecutwfc": 60.0,
            "ecutrho": 480.0,
            "occupations": "smearing",
            "smearing": "cold",
            "degauss": 0.01,
            "nspin": nspin,
        },
        "electrons": {
            "conv_thr": 1.2e-9,
            "mixing_beta": 0.4,
            "electron_maxstep": 120,
        },
    }
    if calc in ("relax", "vc-relax"):
        inp["ions"] = {"ion_dynamics": "bfgs"}
    if calc == "vc-relax":
        inp["cell"] = {"cell_dynamics": "bfgs", "press_conv_thr": 0.5}
    return inp


def first_index(atoms, symbol):
    for i, a in enumerate(atoms):
        if a.symbol == symbol:
            return i
    raise ValueError(f"No {symbol} atom in cell")


def substitute(atoms, idx, new_symbol):
    out = atoms.copy()
    out.symbols[idx] = new_symbol
    return out


def write_qe(path, atoms, prefix, calc, kpts, nspin=1, magnetic_species=None):
    inp = base_input(prefix, calc=calc, nspin=nspin)
    if nspin == 2 and magnetic_species:
        moms = np.zeros(len(atoms))
        for i, a in enumerate(atoms):
            if a.symbol in magnetic_species:
                moms[i] = magnetic_species[a.symbol]
        atoms = atoms.copy()
        atoms.set_initial_magnetic_moments(moms)
    pseudos = {sym: PSEUDOS[sym] for sym in set(atoms.get_chemical_symbols())}
    with open(path, "w") as f:
        write_espresso_in(f, atoms, input_data=inp, pseudopotentials=pseudos, kpts=kpts)
    print(f"  wrote {path.name}: {len(atoms)} atoms, calc={calc}, nspin={nspin}, kpts={kpts}")


def main():
    print("Reading source structures from mgh2-cif/")
    mg2ni     = read(CIFDIR / "mg2ni.pwi")
    mg2nih4_28 = read(CIFDIR / "mg2nih4-28.pwi")
    print(f"  Mg2Ni:      {len(mg2ni)} atoms")
    print(f"  Mg2NiH4-28: {len(mg2nih4_28)} atoms")

    assert len(mg2ni) == 18
    assert len(mg2nih4_28) == 28

    # --- Nb-doped: 1 Mg -> Nb in each ---
    mg2ni_Nb   = substitute(mg2ni,      first_index(mg2ni,      "Mg"), "Nb")
    mg2nih4_Nb = substitute(mg2nih4_28, first_index(mg2nih4_28, "Mg"), "Nb")
    write_qe(INDIR / "mg2ni_Nb.in",   mg2ni_Nb,   "mg2ni_Nb",   calc="relax", kpts=(8, 8, 4))
    write_qe(INDIR / "mg2nih4_Nb.in", mg2nih4_Nb, "mg2nih4_Nb", calc="relax", kpts=(8, 8, 6))

    # --- Nb+Fe co-doped: also 1 Ni -> Fe ---
    mg2ni_NbFe   = substitute(mg2ni_Nb,   first_index(mg2ni_Nb,   "Ni"), "Fe")
    mg2nih4_NbFe = substitute(mg2nih4_Nb, first_index(mg2nih4_Nb, "Ni"), "Fe")
    write_qe(INDIR / "mg2ni_NbFe.in",   mg2ni_NbFe,   "mg2ni_NbFe",   calc="relax",
             kpts=(8, 8, 4), nspin=2, magnetic_species={"Fe": 2.0})
    write_qe(INDIR / "mg2nih4_NbFe.in", mg2nih4_NbFe, "mg2nih4_NbFe", calc="relax",
             kpts=(8, 8, 6), nspin=2, magnetic_species={"Fe": 2.0})

    # --- Auxiliary references for the Nb2O5 reduction reaction ---
    # MgO (rocksalt, primitive 2 atoms)
    mgo = crystal(["Mg", "O"], [(0, 0, 0), (0.5, 0.5, 0.5)],
                  spacegroup=225, cellpar=[4.21, 4.21, 4.21, 90, 90, 90], primitive_cell=True)
    write_qe(INDIR / "mgo.in", mgo, "mgo", calc="scf", kpts=(8, 8, 8))

    # Nb2O5 (synthetic 14-atom orthorhombic cell, vc-relax to settle)
    a, b, c = 6.36, 4.88, 6.20
    nb_frac = [
        [0.20, 0.25, 0.05], [0.70, 0.75, 0.05],
        [0.30, 0.25, 0.55], [0.80, 0.75, 0.55],
    ]
    o_frac = [
        [0.00, 0.25, 0.25], [0.50, 0.75, 0.25],
        [0.50, 0.25, 0.25], [0.00, 0.75, 0.25],
        [0.20, 0.75, 0.05], [0.70, 0.25, 0.05],
        [0.30, 0.75, 0.55], [0.80, 0.25, 0.55],
        [0.20, 0.25, 0.85], [0.70, 0.75, 0.85],
    ]
    nb2o5 = Atoms(symbols=["Nb"] * 4 + ["O"] * 10,
                  scaled_positions=nb_frac + o_frac,
                  cell=[[a, 0, 0], [0, b, 0], [0, 0, c]], pbc=True)
    write_qe(INDIR / "nb2o5.in", nb2o5, "nb2o5", calc="vc-relax", kpts=(4, 4, 4))

    print("\nDoping summary:")
    from collections import Counter
    for name, atoms in [("mg2ni_Nb", mg2ni_Nb), ("mg2nih4_Nb", mg2nih4_Nb),
                        ("mg2ni_NbFe", mg2ni_NbFe), ("mg2nih4_NbFe", mg2nih4_NbFe)]:
        print(f"  {name:18s}: {dict(Counter(atoms.get_chemical_symbols()))}")


if __name__ == "__main__":
    main()
