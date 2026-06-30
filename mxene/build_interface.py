#!/usr/bin/env python3
"""
Build Ti3C2-MXene slabs and adsorbate/interface QE inputs for the
Mg2Ni:Ti3C2 study.

The relaxed BARE Ti3C2 monolayer (vc-relax 2Dxy, nspin=2, PBE+D3-BJ) is the
anchor: a = 3.0699 Ang, FM ~1.85 uB/cell, E = -394.9086 Ry. Its relaxed cell
and fractional coordinates are hard-coded below so downstream slabs inherit the
correct lattice without re-reading the output.

Usage (after the vacuum-unit diagnostics decide cluster-vs-surface):
    python3 build_interface.py            # writes the chosen inputs to inputs/
"""
from pathlib import Path
import re
import numpy as np
from ase import Atoms

HERE = Path(__file__).resolve().parent
PSEUDO = "/home/x/Workspace/3-hEspesso/pseudo"

# ---- relaxed bare Ti3C2 monolayer (Angstrom) -----------------------------
A_TI3C2 = 3.0699                      # relaxed in-plane lattice constant
C_VAC   = 24.641                      # c with ~20 A vacuum (held fixed in 2Dxy)
CELL_1x1 = np.array([
    [ A_TI3C2,            0.0,               0.0   ],
    [-A_TI3C2/2.0, A_TI3C2*np.sqrt(3)/2.0,   0.0   ],
    [ 0.0,               0.0,              C_VAC   ],
])
# relaxed fractional coordinates (Ti-C-Ti-C-Ti, mirror-symmetric about z=0.5)
FRAC_1x1 = [
    ("Ti", (0.000000, 0.000000, 0.5000000)),
    ("Ti", (0.333333, 0.666667, 0.5943115)),
    ("Ti", (0.666667, 0.333333, 0.4056885)),
    ("C",  (0.666667, 0.333333, 0.5527671)),
    ("C",  (0.333333, 0.666667, 0.4472329)),
]

PSEUDOS = {
    "Ti": "ti_pbe_v1.4.uspp.F.UPF",
    "C":  "C.pbe-n-kjpaw_psl.1.0.0.UPF",
    "Mg": "Mg.pbe-n-kjpaw_psl.0.3.0.UPF",
    "Ni": "ni_pbe_v1.4.uspp.F.UPF",
    "H":  "H.pbe-rrkjus_psl.1.0.0.UPF",
}
MASS = {"Ti": 47.867, "C": 12.011, "Mg": 24.305, "Ni": 58.6934, "H": 1.00794}


def ti3c2_slab(nx=3, ny=3):
    """Relaxed Ti3C2 monolayer replicated to an (nx, ny) lateral supercell."""
    syms = [s for s, _ in FRAC_1x1]
    cell1 = Atoms(symbols=syms, scaled_positions=[f for _, f in FRAC_1x1],
                  cell=CELL_1x1, pbc=True)
    slab = cell1.repeat((nx, ny, 1))
    return slab


def top_z(atoms):
    """Highest atomic z (top surface of the slab), Angstrom."""
    return atoms.get_positions()[:, 2].max()


def write_qe(path, atoms, *, calc="relax", prefix, nspin=2, magn=None,
             kpts=(3, 3, 1), degauss=0.02, dipole=False, max_seconds=21600,
             fixed_below=None):
    """Write a pw.x input. magn: dict{species:mag}. fixed_below: z (Ang) below
    which atoms are frozen (if_pos 0 0 0) to mimic a held substrate."""
    species = []
    for s in atoms.get_chemical_symbols():
        if s not in species:
            species.append(s)
    ntyp = len(species)
    nat = len(atoms)

    lines = []
    lines.append("&CONTROL")
    lines.append(f"    calculation      = '{calc}'")
    lines.append(f"    prefix           = '{prefix}'")
    lines.append(f"    pseudo_dir       = '{PSEUDO}'")
    lines.append("    outdir           = './outdir/'")
    lines.append("    restart_mode     = 'from_scratch'")
    lines.append("    tprnfor          = .true.")
    lines.append("    etot_conv_thr    = 1.0d-5")
    lines.append("    forc_conv_thr    = 1.0d-4")
    lines.append("    nstep            = 200")
    lines.append(f"    max_seconds      = {max_seconds}")
    lines.append("/")
    lines.append("")
    lines.append("&SYSTEM")
    lines.append("    ibrav            = 0")
    lines.append(f"    nat              = {nat}")
    lines.append(f"    ntyp             = {ntyp}")
    lines.append("    ecutwfc          = 60.0")
    lines.append("    ecutrho          = 600.0")
    lines.append("    occupations      = 'smearing'")
    lines.append("    smearing         = 'cold'")
    lines.append(f"    degauss          = {degauss}")
    lines.append(f"    nspin            = {nspin}")
    if nspin == 2 and magn:
        for i, s in enumerate(species, start=1):
            lines.append(f"    starting_magnetization({i}) = {magn.get(s, 0.0)}")
    lines.append("    vdw_corr         = 'DFT-D3'")
    lines.append("    dftd3_version    = 4")
    if dipole:
        lines.append("    tefield          = .true.")
        lines.append("    dipfield         = .true.")
        lines.append("    edir             = 3")
        lines.append("    emaxpos          = 0.95")
        lines.append("    eopreg           = 0.05")
    lines.append("/")
    lines.append("")
    lines.append("&ELECTRONS")
    lines.append("    electron_maxstep = 250")
    lines.append("    conv_thr         = 1.0d-7")
    lines.append("    mixing_beta      = 0.2")
    lines.append("    mixing_mode      = 'local-TF'")
    lines.append("    mixing_ndim      = 12")
    lines.append("    diago_david_ndim = 4")
    lines.append("/")
    lines.append("")
    lines.append("&IONS")
    lines.append("    ion_dynamics     = 'bfgs'")
    lines.append("/")
    lines.append("")
    lines.append("ATOMIC_SPECIES")
    for s in species:
        lines.append(f"  {s:2s}  {MASS[s]:9.4f}  {PSEUDOS[s]}")
    lines.append("")
    lines.append("CELL_PARAMETERS angstrom")
    for v in atoms.get_cell():
        lines.append(f"  {v[0]:16.10f} {v[1]:16.10f} {v[2]:16.10f}")
    lines.append("")
    lines.append("ATOMIC_POSITIONS angstrom")
    pos = atoms.get_positions()
    syms = atoms.get_chemical_symbols()
    for s, p in zip(syms, pos):
        if fixed_below is not None and p[2] < fixed_below:
            tail = "  0 0 0"
        else:
            tail = ""
        lines.append(f"  {s:2s} {p[0]:16.10f} {p[1]:16.10f} {p[2]:16.10f}{tail}")
    lines.append("")
    lines.append("K_POINTS automatic")
    lines.append(f"  {kpts[0]} {kpts[1]} {kpts[2]}   0 0 0")
    lines.append("")
    Path(path).write_text("\n".join(lines))
    print(f"wrote {path}  (nat={nat}, ntyp={ntyp})")


def _parse_pos_block(blk):
    syms, pos = [], []
    for line in blk.splitlines():
        p = line.split()
        if len(p) >= 4 and p[0] in MASS:
            syms.append(p[0]); pos.append([float(x) for x in p[1:4]])
    return syms, np.array(pos)


def read_final_coords(prefix):
    """Relaxed coordinates (Angstrom) for `prefix`: prefer the output's final
    coordinates; fall back to the input geometry if the relax never produced a
    relaxed block (e.g. SCF-non-converged run)."""
    out = HERE / "outputs" / f"{prefix}.out"
    if out.exists():
        txt = out.read_text()
        m = re.search(r"Begin final coordinates(.*?)End final coordinates", txt, re.S)
        if m:
            return _parse_pos_block(m.group(1))
        blocks = re.findall(r"ATOMIC_POSITIONS[^\n]*\n(.*?)(?:\n\s*\n|\nEnd|\Z)", txt, re.S)
        if blocks:
            return _parse_pos_block(blocks[-1])
    # fallback: input geometry
    inp = (HERE / "inputs" / f"{prefix}.in").read_text()
    blk = re.search(r"ATOMIC_POSITIONS[^\n]*\n(.*?)(?:\nK_POINTS|\Z)", inp, re.S).group(1)
    print(f"  [warn] {prefix}: no relaxed coords in output, using input geometry")
    return _parse_pos_block(blk)


def place_on_slab(unit_syms, unit_pos, nx=3, ny=3, gap=2.2, vac=12.0):
    """Center a relaxed unit `gap` Ang above the Ti3C2 (nx,ny) slab top, with
    `vac` Ang vacuum on each side along z."""
    slab = ti3c2_slab(nx, ny)
    sp = slab.get_positions()
    tz = sp[:, 2].max()
    topmask = sp[:, 2] > tz - 0.2
    cx, cy = sp[topmask, 0].mean(), sp[topmask, 1].mean()
    u = unit_pos.copy()
    u[:, 0] += cx - u[:, 0].mean()
    u[:, 1] += cy - u[:, 1].mean()
    u[:, 2] += (tz + gap) - u[:, 2].min()
    combined = slab + Atoms(symbols=unit_syms, positions=u)
    pos = combined.get_positions()
    zmin, zmax = pos[:, 2].min(), pos[:, 2].max()
    pos[:, 2] += vac - zmin
    combined.set_positions(pos)
    cell = combined.get_cell()
    cell[2] = [0.0, 0.0, (zmax - zmin) + 2 * vac]
    combined.set_cell(cell)
    return combined


if __name__ == "__main__":
    slab = ti3c2_slab(3, 3)
    print(f"Ti3C2 3x3 slab: {len(slab)} atoms, top z = {top_z(slab):.3f} Ang, "
          f"cell a = {np.linalg.norm(slab.get_cell()[0]):.3f} Ang")

    MAGN = {"Ti": 0.4, "C": 0.0, "Ni": 0.3, "Mg": 0.0, "H": 0.0}
    jobs = [
        ("mg2ni_on_ti3c2",   "outputs/mg2ni_unit.out"),
        ("mg2nih4_on_ti3c2", "outputs/mg2nih4_unit.out"),
    ]
    for name, src in jobs:
        syms, pos = read_final_coords(src)
        atoms = place_on_slab(syms, pos, 3, 3, gap=2.2, vac=12.0)
        zfix = top_z(ti3c2_slab(3, 3)) + 12.0 + 0.8   # freeze the whole MXene slab
        write_qe(f"inputs/{name}.in", atoms, prefix=name, nspin=2, magn=MAGN,
                 kpts=(4, 4, 1), degauss=0.02, dipole=False,
                 fixed_below=zfix, max_seconds=36000)
        print(f"  {name}: {len(atoms)} atoms, c = {atoms.get_cell()[2][2]:.2f} Ang, "
              f"frozen-below z = {zfix:.2f}")
