"""Generate every Quantum ESPRESSO input for the V2O5 gas-sensing study.

Settings are fixed by gas/DESIGN.md: rev-vdW-DF2 (input_dft='vdW-DF2-b86r'),
60/600 Ry, NO Hubbard U, dipole-corrected (001) slabs, per-system spin
(NO2 -> nspin=2). Structure building uses ASE.

Pipeline: build_bulk() + build_molecules() can run immediately. build_slab()
and build_adsorption() need the relaxed bulk (outputs/bulk.out), so main()
emits them only once the bulk vc-relax is done -- re-run after the bulk.
"""
from pathlib import Path
import math
import numpy as np
from ase import Atoms
from ase.io import read
from ase.build import molecule

HERE = Path(__file__).resolve().parent
PSEUDO_DIR = "/home/x/Workspace/3-hEspesso/pseudo"
INPUTS = HERE / "inputs"
for _d in ("inputs", "outputs", "logs", "outdir"):
    (HERE / _d).mkdir(exist_ok=True)

PSEUDOS = {
    "V": "v_pbe_v1.4.uspp.F.UPF",
    "O": "O.pbe-n-kjpaw_psl.0.1.UPF",
    "C": "C.pbe-n-kjpaw_psl.1.0.0.UPF",
    "H": "H.pbe-rrkjus_psl.1.0.0.UPF",
    "N": "N.pbe-n-radius_5.UPF",
    "S": "s_pbe_v1.4.uspp.F.UPF",
}
MASS = {"V": 50.9415, "O": 15.999, "C": 12.011, "H": 1.008, "N": 14.007, "S": 32.06}

GASES = {  # key: (ASE G2 name, nspin)
    "acetone": ("CH3COCH3", 1), "nh3": ("NH3", 1), "c2h4": ("C2H4", 1),
    "h2s": ("H2S", 1), "co": ("CO", 1), "no2": ("NO2", 2),
}
GAS_NAT = {"acetone": 10, "nh3": 4, "c2h4": 6, "h2s": 3, "co": 2, "no2": 3}


def axis_align(atoms):
    """Axis-aligned orthorhombic copy with a=longest, b=shortest, c=middle.
    For alpha-V2O5 the middle-length axis is the interlayer (stacking) axis,
    so this puts the vdW gap along c -- ready for (001) slab construction."""
    L = atoms.cell.lengths()
    o = np.argsort(L)                 # indices: short, mid, long
    order = [o[2], o[0], o[1]]        # -> long(a), short(b), mid(c)
    frac = atoms.get_scaled_positions()
    out = atoms.copy()
    out.set_cell(np.diag(L[order]), scale_atoms=False)
    out.set_scaled_positions(frac[:, order])
    return out


def write_qe(path, atoms, *, calc, prefix, kpts, nspin=1, magn=None,
             fixed_below=None, dipole=False):
    """Write one pw.x input. calc in {scf,relax,vc-relax}; magn={sym:mag};
    fixed_below freezes atoms with z < value (A); dipole adds the z dipole field."""
    syms = atoms.get_chemical_symbols()
    species = list(dict.fromkeys(syms))
    nat, ntyp = len(atoms), len(species)
    L = ["&CONTROL", f"    calculation      = '{calc}'",
         f"    prefix           = '{prefix}'",
         f"    pseudo_dir       = '{PSEUDO_DIR}'",
         "    outdir           = './outdir/'",
         "    restart_mode     = 'from_scratch'",
         "    tprnfor          = .true.",
         f"    tstress          = .{'false' if dipole else 'true'}.",
         "    etot_conv_thr    = 1.0d-5", "    forc_conv_thr    = 1.0d-4",
         "    nstep            = 200", "    max_seconds      = 36000", "/"]
    L += ["&SYSTEM", "    ibrav            = 0", f"    nat              = {nat}",
          f"    ntyp             = {ntyp}", "    ecutwfc          = 60.0",
          "    ecutrho          = 600.0", "    occupations      = 'smearing'",
          "    smearing         = 'cold'", "    degauss          = 0.01",
          f"    nspin            = {nspin}"]
    if nspin == 2:
        for i, s in enumerate(species, 1):
            L.append(f"    starting_magnetization({i}) = {(magn or {}).get(s, 0.0)}")
    # PBE + DFT-D3(BJ): this GPU build does NOT accelerate the rev-vdW-DF2
    # non-local kernel (runs CPU-bound at ~70 s/iter), making 24 slab relaxes
    # infeasible. D3 is GPU-accelerated (~26 s/iter) and the repo-standard
    # dispersion; the sensitivity RANKING is robust to the choice.
    L += ["    vdw_corr         = 'DFT-D3'", "    dftd3_version    = 4"]
    if dipole:  # one-sided slab -> ESM open boundaries along z. This GPU build
                # lacks tefield/dipfield in its &SYSTEM namelist; ESM bc1 (vacuum
                # both sides) removes the spurious inter-image field and gives a
                # correct work function. Requires nk3=1 and tstress off.
        L += ["    assume_isolated  = 'esm'", "    esm_bc           = 'bc1'"]
    L += ["/"]
    L += ["&ELECTRONS", "    electron_maxstep = 250", "    conv_thr         = 1.0d-9",
          "    mixing_beta      = 0.2", "    mixing_mode      = 'local-TF'",
          "    diagonalization  = 'david'", "/"]
    if calc in ("relax", "vc-relax"):
        L += ["&IONS", "    ion_dynamics     = 'bfgs'", "/"]
    if calc == "vc-relax":
        L += ["&CELL", "    cell_dynamics    = 'bfgs'", "    press_conv_thr   = 0.5", "/"]
    L += ["ATOMIC_SPECIES"]
    for s in species:
        L.append(f"  {s:2s}  {MASS[s]:9.4f}  {PSEUDOS[s]}")
    L += ["", "CELL_PARAMETERS angstrom"]
    for v in atoms.get_cell():
        L.append(f"  {v[0]:16.10f} {v[1]:16.10f} {v[2]:16.10f}")
    L += ["", "ATOMIC_POSITIONS angstrom"]
    for s, p in zip(syms, atoms.get_positions()):
        tail = "  0 0 0" if (fixed_below is not None and p[2] < fixed_below) else ""
        L.append(f"  {s:2s} {p[0]:16.10f} {p[1]:16.10f} {p[2]:16.10f}{tail}")
    L += ["", "K_POINTS automatic", f"  {kpts[0]} {kpts[1]} {kpts[2]} 0 0 0", ""]
    Path(path).write_text("\n".join(L))
    print(f"wrote {Path(path).name}  (nat={nat}, ntyp={ntyp}, {calc}, nspin={nspin})")


# --- bulk -------------------------------------------------------------------

def build_bulk():
    atoms = axis_align(read(HERE / "structures" / "V2O5_alpha.cif"))
    write_qe(INPUTS / "bulk.in", atoms, calc="vc-relax", prefix="bulk", kpts=(2, 8, 6))


def relaxed_bulk():
    """Relaxed, axis-aligned alpha-V2O5 from outputs/bulk.out."""
    return axis_align(read(HERE / "outputs" / "bulk.out", index=-1,
                           format="espresso-out"))


# --- slab -------------------------------------------------------------------

def slab_atoms(nlayers=2, nb=3, vacuum=18.0):
    """(001) slab: relaxed bulk replicated nb along b and nlayers along c, with
    `vacuum` A of vacuum along c. alpha-V2O5 layers are vdW-separated and lie
    perpendicular to c, so stacked bulk layers + vacuum IS the (001) surface."""
    bulk = relaxed_bulk()
    slab = bulk.repeat((1, nb, nlayers))
    z = slab.get_positions()[:, 2]
    cell = slab.get_cell()
    cell[2] = [0.0, 0.0, (z.max() - z.min()) + vacuum]
    slab.set_cell(cell, scale_atoms=False)
    # ESM requires atoms within (-Lz/2, +Lz/2], i.e. the slab centered at z=0
    # (not Lz/2), so the open boundaries at the cell edges fall in vacuum.
    zc = slab.get_positions()[:, 2]
    slab.translate([0.0, 0.0, -(zc.min() + zc.max()) / 2.0])
    z = slab.get_positions()[:, 2]
    slab.info["fixed_below"] = float(z.min() + (z.max() - z.min()) / (2 * nlayers))
    return slab


def relaxed_slab_1L():
    """The RELAXED monolayer from outputs/slab_1L.out, re-centered at z=0 for ESM.
    Starting each adsorption complex from the already-relaxed surface (instead of
    the bulk-cut one) cuts ~3x the BFGS steps -- same energy minimum, much faster.
    Freezes the bottom sublayer, as slab_atoms(1) does."""
    s = read(HERE / "outputs" / "slab_1L.out", index=-1, format="espresso-out")
    z = s.get_positions()[:, 2]
    s.translate([0.0, 0.0, -(z.min() + z.max()) / 2.0])
    z = s.get_positions()[:, 2]
    s.info["fixed_below"] = float(z.min() + (z.max() - z.min()) / 2.0)
    return s


def build_slab():
    for name, nl in (("slab", 2), ("slab_1L", 1)):
        s = slab_atoms(nlayers=nl)
        write_qe(INPUTS / f"{name}.in", s, calc="relax", prefix=name,
                 kpts=(2, 2, 1), dipole=True, fixed_below=s.info["fixed_below"])


# --- gas molecules ----------------------------------------------------------

def _h2s():
    """H2S geometry (ASE G2 lacks it): S-H 1.34 A, H-S-H 92.1 deg."""
    r, half = 1.34, math.radians(92.1) / 2
    return Atoms("SHH", positions=[(0, 0, 0),
                 (r * math.sin(half), 0, r * math.cos(half)),
                 (-r * math.sin(half), 0, r * math.cos(half))])


def build_molecules(box=15.0):
    for key, (g2, nspin) in GASES.items():
        m = _h2s() if key == "h2s" else molecule(g2)
        m.set_cell([box, box, box])
        m.center()
        m.pbc = True
        magn = {"N": 0.5} if nspin == 2 else None
        write_qe(INPUTS / f"gas_{key}.in", m, calc="relax", prefix=f"gas_{key}",
                 kpts=(1, 1, 1), nspin=nspin, magn=magn)


# --- adsorption -------------------------------------------------------------

def surface_sites(slab, n=3):
    """Top-surface anchors: highest O (vanadyl), a neighbouring O (bridging),
    and the highest V. Returns up to n (x,y,z) in Angstrom."""
    pos = slab.get_positions()
    syms = slab.get_chemical_symbols()
    ztop = pos[:, 2].max()
    top = [i for i in range(len(slab)) if pos[i, 2] > ztop - 1.6]
    Otop = sorted([i for i in top if syms[i] == "O"], key=lambda i: -pos[i, 2])
    Vtop = sorted([i for i in top if syms[i] == "V"], key=lambda i: -pos[i, 2])
    picks = []
    if Otop:
        picks.append(Otop[0])
    if len(Otop) > 1:
        picks.append(Otop[1])
    if Vtop:
        picks.append(Vtop[0])
    return [tuple(pos[i]) for i in picks[:n]]


def _relaxed_molecule(key):
    return read(HERE / "outputs" / f"gas_{key}.out", index=-1, format="espresso-out")


# Atom each molecule is oriented to present to the surface (lone-pair / donor).
# CO and NO2 bind ambiguously, so both candidate atoms are tried; C2H4 is a
# pi donor placed flat (None = keep relaxed orientation).
ANCHORS = {"acetone": ["O"], "nh3": ["N"], "h2s": ["S"],
           "c2h4": [None], "co": ["C", "O"], "no2": ["N", "O"]}


def _orient_anchor_down(mol, element):
    """Copy of `mol` rotated so the first `element` atom is the lowest point and
    the rest of the molecule extends upward -- that atom faces the surface.
    Returns (oriented_atoms, anchor_index)."""
    m = mol.copy()
    p = m.get_positions()
    idx = next(i for i, s in enumerate(m.get_chemical_symbols()) if s == element)
    v = p.mean(axis=0) - p[idx]                # anchor -> centroid
    if np.linalg.norm(v) > 1e-3:
        m.rotate(v, (0, 0, 1), center=p[idx])  # point that direction along +z
    return m, idx


def build_adsorption(gap=2.2):
    """Place each gas over 3 distinct surface sites (vanadyl O, bridging O, V),
    oriented binding-atom-down; CO and NO2 get both candidate anchors. Indices
    run ads_<key>_s0.. ; analyze.py keeps the lowest-E_ads one per gas."""
    slab = relaxed_slab_1L()   # start from the RELAXED surface -> ~3x fewer steps
    fb = slab.info["fixed_below"]
    sites = surface_sites(slab, 3)
    for key, (_, nspin) in GASES.items():
        mol = _relaxed_molecule(key)
        ci = 0
        for anchor in ANCHORS[key]:
            if anchor is None:
                oriented, aidx = mol.copy(), None
            else:
                oriented, aidx = _orient_anchor_down(mol, anchor)
            for (sx, sy, sz) in sites:
                m = oriented.copy()
                mp = m.get_positions()
                if aidx is None:                       # flat: centroid over site
                    mp[:, 0] += sx - mp[:, 0].mean()
                    mp[:, 1] += sy - mp[:, 1].mean()
                    mp[:, 2] += sz + gap - mp[:, 2].min()
                else:                                  # anchor atom over site
                    mp[:, 0] += sx - mp[aidx, 0]
                    mp[:, 1] += sy - mp[aidx, 1]
                    mp[:, 2] += sz + gap - mp[aidx, 2]
                m.set_positions(mp)
                combined = slab + m
                combined.set_cell(slab.get_cell(), scale_atoms=False)
                magn = {"N": 0.5, "V": 0.0} if nspin == 2 else None
                write_qe(INPUTS / f"ads_{key}_s{ci}.in", combined, calc="relax",
                         prefix=f"ads_{key}_s{ci}", kpts=(2, 2, 1), nspin=nspin,
                         magn=magn, dipole=True, fixed_below=fb)
                ci += 1


def main():
    build_bulk()
    build_molecules()
    bulk_out = HERE / "outputs" / "bulk.out"
    if bulk_out.exists() and "JOB DONE" in bulk_out.read_text(errors="replace"):
        build_slab()
        build_adsorption()
        print("built slab + adsorption inputs (bulk is relaxed)")
    else:
        print("bulk not relaxed yet -> built bulk + molecules only; "
              "re-run after bulk.out has JOB DONE for slab + adsorption")


if __name__ == "__main__":
    main()
