"""
Generate 13 internally-consistent Quantum ESPRESSO inputs in mg-nico/inputs/.

All structures use:
    - PBE + DFT-D3 dispersion (no Hubbard U anywhere -- energies must be
      comparable across all six pathways)
    - ecutwfc = 60 Ry, ecutrho = 480 Ry
    - PAW (Mg, H) + USPP (Ni, Co) pseudopotentials
    - vc-relax (cell + ions) for solids; relax-only for H2 in vacuum
    - cold smearing 0.01 Ry for solids; fixed occupations for the H2 molecule
    - nspin = 1 throughout (see note in write_input)
    - conv_thr = 1e-9, forc_conv_thr = 1e-4 Ry/Bohr

The thirteen structures span six reactions, two hosts x (pure, +Ni, +Co):

    Host Mg (matched 2x2x2 supercells, 16 metal sites, 1 dopant = 6.25%):
      R1  pure      :  Mg       +    H2  ->  MgH2
      R2  Ni-doped  :  Mg15Ni   + 16 H2  ->  Mg15NiH32
      R3  Co-doped  :  Mg15Co   + 16 H2  ->  Mg15CoH32

    Host Mg2Ni (18-atom hexagonal metal cell, 28-atom monoclinic hydride
    cell; 1 catalyst atom on the documented Mg substitution site):
      R4  pure      :  (8/12) Mg12Ni6        + 8 H2  ->  Mg8Ni4H16
      R5  Ni-doped  :  defect shift from Mg11Ni7 / Mg7Ni5H16
      R6  Co-doped  :  defect shift from Mg11CoNi6 / Mg7CoNi4H16

R2/R3 share the same supercell topology, as do R5/R6, so the Ni and Co
columns differ only in dopant chemistry -- no supercell, k-density, or
functional artefacts confound the comparison.
"""

from pathlib import Path

HERE       = Path(__file__).resolve().parent
INPUTS_DIR = HERE / "inputs"
INPUTS_DIR.mkdir(exist_ok=True)

PSEUDO_DIR = "/home/x/Workspace/3-hEspesso/pseudo"
PSEUDO_MG  = "Mg.pbe-n-kjpaw_psl.0.3.0.UPF"
PSEUDO_H   = "H.pbe-rrkjus_psl.1.0.0.UPF"
PSEUDO_NI  = "ni_pbe_v1.4.uspp.F.UPF"
PSEUDO_CO  = "Co_pbe_v1.2.uspp.F.UPF"

# Cell parameters (experimental lattice constants, will be vc-relaxed)
MG_A, MG_C       = 3.20, 5.21       # HCP Mg
MGH2_A, MGH2_C   = 4.5168, 3.0205   # rutile MgH2 (P4_2/mnm)
H2_BOX_BOHR      = 20.0             # 20 Bohr cubic vacuum
H2_BOND_ANG      = 0.7414           # H-H equilibrium bond length

# ---------------------------------------------------------------------------
# Structures: Mg host
# ---------------------------------------------------------------------------

def mg_hcp_primitive():
    """2-atom HCP primitive in fractional (crystal) coords."""
    return [("Mg", 1/3, 2/3, 1/4),
            ("Mg", 2/3, 1/3, 3/4)]

def mgh2_rutile():
    """Rutile MgH2: 6 atoms in tetragonal cell. u = 0.304 internal."""
    u = 0.304
    return [("Mg", 0.0, 0.0, 0.0),
            ("Mg", 0.5, 0.5, 0.5),
            ("H",   u,   u,  0.0),
            ("H",  -u,  -u,  0.0),
            ("H",  0.5+u, 0.5-u, 0.5),
            ("H",  0.5-u, 0.5+u, 0.5)]

def make_mg_supercell_2x2x2():
    """Build 2x2x2 supercell of HCP Mg in CARTESIAN angstrom (16 atoms).
    Hex cell vectors: a1=(a,0,0), a2=(-a/2, a*sqrt(3)/2, 0), a3=(0,0,c).
    """
    import math
    a, c = MG_A, MG_C
    a1 = (a, 0.0, 0.0)
    a2 = (-a/2, a * math.sqrt(3)/2, 0.0)
    a3 = (0.0, 0.0, c)
    basis = [(1/3, 2/3, 1/4), (2/3, 1/3, 3/4)]
    atoms = []
    for ix in range(2):
        for iy in range(2):
            for iz in range(2):
                for fx, fy, fz in basis:
                    cx = (fx + ix) * a1[0] + (fy + iy) * a2[0] + (fz + iz) * a3[0]
                    cy = (fx + ix) * a1[1] + (fy + iy) * a2[1] + (fz + iz) * a3[1]
                    cz = (fx + ix) * a1[2] + (fy + iy) * a2[2] + (fz + iz) * a3[2]
                    atoms.append(("Mg", cx, cy, cz))
    sc1 = tuple(2*x for x in a1)
    sc2 = tuple(2*x for x in a2)
    sc3 = tuple(2*x for x in a3)
    return atoms, (sc1, sc2, sc3)

def make_mgh2_supercell_2x2x2():
    """Build 2x2x2 supercell of rutile MgH2 in CARTESIAN angstrom (48 atoms)."""
    a, c = MGH2_A, MGH2_C
    a1 = (a, 0.0, 0.0)
    a2 = (0.0, a, 0.0)
    a3 = (0.0, 0.0, c)
    basis = mgh2_rutile()
    atoms = []
    for ix in range(2):
        for iy in range(2):
            for iz in range(2):
                for sp, fx, fy, fz in basis:
                    cx = (fx + ix) * a1[0]
                    cy = (fy + iy) * a2[1]
                    cz = (fz + iz) * a3[2]
                    atoms.append((sp, cx, cy, cz))
    return atoms, ((2*a, 0, 0), (0, 2*a, 0), (0, 0, 2*c))

def substitute_first(atoms, src_species, new_species):
    """Replace the first `src_species` atom with `new_species`.

    One substitutional dopant per cell. In HCP Mg and rutile MgH2 every
    Mg site is symmetry-equivalent, so 'first' is the unique choice."""
    out = list(atoms)
    for i, (sp, *_) in enumerate(out):
        if sp == src_species:
            out[i] = (new_species, *out[i][1:])
            return out
    raise ValueError(f"no {src_species} atom found")

# ---------------------------------------------------------------------------
# Structures: Mg2Ni host
#
# Geometries taken from the mg2ni-nb2o5fe study (which in turn used the
# relaxed mgh2-cif structures): 18-atom hexagonal Mg2Ni (P6_222 topology)
# and 28-atom monoclinic LT-Mg2NiH4 (C2/c, Z=4) in a primitive setting.
# The Nb substitution site of that study is restored to Mg and reused as
# the catalyst substitution site here (index 0), so dopant placement is
# identical across studies. Cells are vc-relaxed from these coordinates.
# ---------------------------------------------------------------------------

MG2NI_CELL = (
    (5.1974510419, 0.0000000000, 0.0000000000),
    (-2.5987255209, 4.5011246372, 0.0000000000),
    (0.0000000000, 0.0000000000, 13.2024169613),
)
MG2NI_ATOMS = [   # index 0 = dopant substitution site (the Nb site of the Nb2O5/Fe study)
    ("Mg",  -1.3173663055,   3.7613315492,  11.0020141344),
    ("Mg",   0.0000000000,   1.4795861759,   6.6012084806),
    ("Mg",   3.9160918265,   0.7397930880,  11.0020141344),
    ("Mg",   2.5987255209,   3.0215384613,   6.6012084806),
    ("Mg",   1.3173663055,   3.7613315492,   2.2004028269),
    ("Mg",   1.2813592154,   0.7397930880,   2.2004028269),
    ("Mg",   1.2993627605,   2.2505623186,   9.4710187506),
    ("Mg",   2.5987255209,   0.0000000000,   8.1322038644),
    ("Mg",   1.2993627605,   2.2505623186,  12.5330095181),
    ("Mg",   2.5987255209,   0.0000000000,   5.0702130969),
    ("Mg",   3.8980882814,   2.2505623186,   3.7313982107),
    ("Mg",   3.8980882814,   2.2505623186,   0.6694074432),
    ("Ni",   1.2993627605,   2.2505623186,   4.4008056537),
    ("Ni",   2.5987255209,   0.0000000000,   0.0000000000),
    ("Ni",   3.8980882814,   2.2505623186,   8.8016113075),
    ("Ni",   0.0000000000,   0.0000000000,   4.4008056537),
    ("Ni",   0.0000000000,   0.0000000000,   0.0000000000),
    ("Ni",   0.0000000000,   0.0000000000,   8.8016113075),
]

MG2NIH4_CELL = (
    (7.8519923697, 0.0000000000, 0.0000000000),
    (-5.2504433367, 5.8383755397, 0.0000000000),
    (-2.3428084713, 1.0439429721, 5.9460839518),
)
MG2NIH4_ATOMS = [   # index 0 = dopant substitution site (the Nb site of the Nb2O5/Fe study)
    ("Mg",  -4.2792065963,   5.4955229079,   5.4607179434),
    ("Mg",  -3.1746666110,   4.8234952058,   2.4876759674),
    ("Mg",   3.4334071727,   2.0588233060,   3.4584079844),
    ("Mg",   4.5379471580,   1.3867956039,   0.4853660086),
    ("Mg",   0.7825802604,   6.4825016222,   4.4595629638),
    ("Mg",  -0.5238396986,   0.3998168897,   1.4865209879),
    ("Mg",  -0.5224029966,   3.5538686307,   4.4595629638),
    ("Mg",   0.7811435583,   3.3284498812,   1.4865209879),
    ("Ni",  -1.7226015765,   6.1553205873,   5.4656626820),
    ("Ni",  -1.9557195328,   2.4813316280,   2.4926207061),
    ("Ni",   2.2144600944,   4.4009868838,   3.4534632458),
    ("Ni",   1.9813421383,   0.7269979244,   0.4804212699),
    ("H",   -0.5032506657,   5.1896764017,   5.6320551105),
    ("H",   -0.4223601390,   2.2203827687,   2.6590131347),
    ("H",    0.6811007008,   4.6619357431,   3.2870708173),
    ("H",    0.7619912274,   1.6926421103,   0.3140288413),
    ("H",    5.1967385498,   1.0942273572,   2.5297971576),
    ("H",    4.0914656533,   1.7646102659,   5.5028391336),
    ("H",   -3.8327250917,   5.1177082459,   0.4432448183),
    ("H",   -4.9379979881,   5.7880911547,   3.4162867941),
    ("H",   -2.5682459864,   6.0096553346,   4.1415091751),
    ("H",   -2.4128722450,   3.2075163369,   1.1684671992),
    ("H",    2.6716128067,   3.6748021749,   4.7776167528),
    ("H",    2.8269865480,   0.8726631772,   1.8045747769),
    ("H",    0.5922569907,   1.9608958058,   5.2073775167),
    ("H",    2.7109509201,   3.5648242766,   2.2343355407),
    ("H",   -2.4522103584,   3.3174942352,   3.7117484112),
    ("H",   -0.3335164290,   4.9214227061,   0.7387064353),
]

def substitute_site0(atoms, new_species):
    """Replace the designated substitution site (index 0) with the catalyst."""
    out = list(atoms)
    out[0] = (new_species, *out[0][1:])
    return out

def substitute_first_ni(atoms, new_species):
    """Replace the first Ni atom with the catalyst (Ni-site substitution).

    Yoon et al., J. Magnesium Alloys 12 (2024) 4574 find Co prefers the Ni
    site over the Mg site in stoichiometric Mg2Ni (E_Mg->Ni = -0.53 eV), so
    the Co pathway is computed for both substitution sites."""
    out = list(atoms)
    for i, (sp, *_) in enumerate(out):
        if sp == "Ni":
            out[i] = (new_species, *out[i][1:])
            return out
    raise ValueError("no Ni atom found")

# ---------------------------------------------------------------------------
# Input file emission
# ---------------------------------------------------------------------------

CONTROL_TEMPLATE = """&CONTROL
    calculation      = '{calc}'
    prefix           = '{prefix}'
    pseudo_dir       = '{pseudo_dir}'
    outdir           = './outdir/'
    restart_mode     = 'from_scratch'
    tprnfor          = .true.
    tstress          = .true.
    etot_conv_thr    = 1.0d-5
    forc_conv_thr    = 1.0d-4
    nstep            = 200
    max_seconds      = 21600
/
"""

SYSTEM_TEMPLATE = """&SYSTEM
    ibrav            = 0
    nat              = {nat}
    ntyp             = {ntyp}
    ecutwfc          = 60.0
    ecutrho          = 480.0
    occupations      = '{occ}'
    {smearing_lines}
    {spin_lines}
    {symmetry_lines}
    vdw_corr         = 'DFT-D3'
    dftd3_version    = 4
/
"""

ELECTRONS_TEMPLATE = """&ELECTRONS
    electron_maxstep = {maxstep}
    conv_thr         = 1.0d-9
    mixing_beta      = {mix}
    mixing_mode      = '{mix_mode}'
    mixing_ndim      = {mix_ndim}
    diagonalization  = 'david'
    diago_david_ndim = 4
/
"""

IONS_TEMPLATE = "&IONS\n    ion_dynamics = 'bfgs'\n/\n"
CELL_TEMPLATE = "&CELL\n    cell_dynamics = 'bfgs'\n    press_conv_thr = 0.5\n/\n"


# Cells whose electronic ground state is spin-polarised (Co retains a moment).
# Determined empirically via single-point nspin=2 spin-checks (make_spinchk.py):
# Co is magnetic in the metals and in MgH2, but quenches at the H-coordinated
# Mg/Ni sites of Mg2NiH4 (mg2nih4_co, mg2nih4_co_nisite stay nonmagnetic).
# Ni is nonmagnetic everywhere. See README "Spin treatment".
SPIN2_CELLS = {"mgco", "mgh2co", "mg2ni_co", "mg2ni_co_nisite"}


def write_input(name, atoms, cell_vecs, calc, kpts, nosym=False):
    """Write a complete pw.x input file to inputs/{name}.in.

    atoms: list of (species, x, y, z) in angstrom
    cell_vecs: 3 cell vectors in angstrom (rows of matrix)
    calc: 'scf' | 'relax' | 'vc-relax'
    kpts: (nk1, nk2, nk3)
    nosym: disable symmetry detection (set for cells with a substitutional
           dopant -- the dopant breaks the parent symmetry and vc-relax can
           otherwise trip on "not orthogonal operation" as the cell distorts)
    """
    species = []
    for sp, *_ in atoms:
        if sp not in species:
            species.append(sp)
    ntyp = len(species)
    nat = len(atoms)

    is_h2_only = species == ["H"]
    has_tm = any(sp in ("Ni", "Co") for sp in species)
    if is_h2_only:
        occ = "fixed"
        smearing_lines = ""
        mix = 0.7
        maxstep, mix_mode, mix_ndim = 100, "plain", 8
    else:
        occ = "smearing"
        smearing_lines = "smearing         = 'cold'\n    degauss          = 0.01"
        mix = 0.3 if has_tm else 0.4
        maxstep, mix_mode, mix_ndim = (200, "plain", 12) if has_tm else (120, "plain", 8)

    # Spin: per-cell ground state. Magnetic Co cells get nspin=2 + a starting
    # moment on the TM species; all others stay nspin=1 (Ni nonmagnetic, Co
    # quenched in the Mg2NiH4 hydrides). See SPIN2_CELLS above.
    if name in SPIN2_CELLS:
        mag = "".join(f"\n    starting_magnetization({i+1}) = {0.4 if sp == 'Co' else 0.1}"
                      for i, sp in enumerate(species) if sp in ("Co", "Ni"))
        spin_lines = "nspin            = 2" + mag
    else:
        spin_lines = ""

    symmetry_lines = "nosym            = .true.\n    noinv            = .true." if nosym else ""

    control  = CONTROL_TEMPLATE.format(
        calc=calc, prefix=name, pseudo_dir=PSEUDO_DIR)
    system   = SYSTEM_TEMPLATE.format(
        nat=nat, ntyp=ntyp, occ=occ,
        smearing_lines=smearing_lines, spin_lines=spin_lines,
        symmetry_lines=symmetry_lines)
    electrons = ELECTRONS_TEMPLATE.format(
        mix=mix, maxstep=maxstep, mix_mode=mix_mode, mix_ndim=mix_ndim)

    parts = [control, system, electrons]
    if calc in ("relax", "vc-relax"):
        parts.append(IONS_TEMPLATE)
    if calc == "vc-relax":
        parts.append(CELL_TEMPLATE)

    parts.append("ATOMIC_SPECIES")
    for sp in species:
        if sp == "Mg":   parts.append(f"  Mg  24.305    {PSEUDO_MG}")
        elif sp == "H":  parts.append(f"  H    1.00794  {PSEUDO_H}")
        elif sp == "Ni": parts.append(f"  Ni  58.6934   {PSEUDO_NI}")
        elif sp == "Co": parts.append(f"  Co  58.9332   {PSEUDO_CO}")
    parts.append("")

    parts.append("CELL_PARAMETERS angstrom")
    for v in cell_vecs:
        parts.append(f"  {v[0]:16.10f}  {v[1]:16.10f}  {v[2]:16.10f}")
    parts.append("")

    parts.append("ATOMIC_POSITIONS angstrom")
    for sp, x, y, z in atoms:
        parts.append(f"  {sp:2s}  {x:16.10f}  {y:16.10f}  {z:16.10f}")
    parts.append("")

    parts.append("K_POINTS automatic")
    parts.append(f"  {kpts[0]} {kpts[1]} {kpts[2]} 0 0 0")
    parts.append("")

    (INPUTS_DIR / f"{name}.in").write_text("\n".join(parts))


# ---------------------------------------------------------------------------
# Build all 13 inputs
# ---------------------------------------------------------------------------

def build_mg():
    import math
    a, c = MG_A, MG_C
    cell = ((a, 0, 0), (-a/2, a*math.sqrt(3)/2, 0), (0, 0, c))
    basis = mg_hcp_primitive()
    atoms = []
    for sp, fx, fy, fz in basis:
        cx = fx*cell[0][0] + fy*cell[1][0] + fz*cell[2][0]
        cy = fx*cell[0][1] + fy*cell[1][1] + fz*cell[2][1]
        cz = fx*cell[0][2] + fy*cell[1][2] + fz*cell[2][2]
        atoms.append((sp, cx, cy, cz))
    write_input("mg", atoms, cell, "vc-relax", (12, 12, 8))

def build_h2():
    L_ang = H2_BOX_BOHR * 0.529177
    cell = ((L_ang, 0, 0), (0, L_ang, 0), (0, 0, L_ang))
    atoms = [("H", L_ang/2 - H2_BOND_ANG/2, L_ang/2, L_ang/2),
             ("H", L_ang/2 + H2_BOND_ANG/2, L_ang/2, L_ang/2)]
    write_input("h2", atoms, cell, "relax", (1, 1, 1))

def build_mgh2():
    a, c = MGH2_A, MGH2_C
    cell = ((a, 0, 0), (0, a, 0), (0, 0, c))
    basis = mgh2_rutile()
    atoms = []
    for sp, fx, fy, fz in basis:
        atoms.append((sp, fx*a, fy*a, fz*c))
    write_input("mgh2", atoms, cell, "vc-relax", (8, 8, 12))

# --- Mg host, 1 dopant (6.25% of metal sites) ---

def build_mgni():
    atoms, cell = make_mg_supercell_2x2x2()
    atoms = substitute_first(atoms, "Mg", "Ni")
    write_input("mgni", atoms, cell, "vc-relax", (6, 6, 4), nosym=True)

def build_mgco():
    atoms, cell = make_mg_supercell_2x2x2()
    atoms = substitute_first(atoms, "Mg", "Co")
    write_input("mgco", atoms, cell, "vc-relax", (6, 6, 4), nosym=True)

def build_mgh2ni():
    atoms, cell = make_mgh2_supercell_2x2x2()
    atoms = substitute_first(atoms, "Mg", "Ni")
    write_input("mgh2ni", atoms, cell, "vc-relax", (4, 4, 6), nosym=True)

def build_mgh2co():
    atoms, cell = make_mgh2_supercell_2x2x2()
    atoms = substitute_first(atoms, "Mg", "Co")
    write_input("mgh2co", atoms, cell, "vc-relax", (4, 4, 6), nosym=True)

# --- Mg2Ni host, pristine + 1 catalyst atom on the documented Mg site ---

def build_mg2ni():
    write_input("mg2ni", MG2NI_ATOMS, MG2NI_CELL, "vc-relax", (8, 8, 4))

def build_mg2ni_ni():
    atoms = substitute_site0(MG2NI_ATOMS, "Ni")
    write_input("mg2ni_ni", atoms, MG2NI_CELL, "vc-relax", (8, 8, 4), nosym=True)

def build_mg2ni_co():
    atoms = substitute_site0(MG2NI_ATOMS, "Co")
    write_input("mg2ni_co", atoms, MG2NI_CELL, "vc-relax", (8, 8, 4), nosym=True)

def build_mg2nih4():
    write_input("mg2nih4", MG2NIH4_ATOMS, MG2NIH4_CELL, "vc-relax", (8, 8, 6))

def build_mg2nih4_ni():
    atoms = substitute_site0(MG2NIH4_ATOMS, "Ni")
    write_input("mg2nih4_ni", atoms, MG2NIH4_CELL, "vc-relax", (8, 8, 6), nosym=True)

def build_mg2nih4_co():
    atoms = substitute_site0(MG2NIH4_ATOMS, "Co")
    write_input("mg2nih4_co", atoms, MG2NIH4_CELL, "vc-relax", (8, 8, 6), nosym=True)

def build_mg2ni_co_nisite():
    atoms = substitute_first_ni(MG2NI_ATOMS, "Co")
    write_input("mg2ni_co_nisite", atoms, MG2NI_CELL, "vc-relax", (8, 8, 4), nosym=True)

def build_mg2nih4_co_nisite():
    atoms = substitute_first_ni(MG2NIH4_ATOMS, "Co")
    write_input("mg2nih4_co_nisite", atoms, MG2NIH4_CELL, "vc-relax", (8, 8, 6), nosym=True)


def main():
    builders = [build_mg, build_h2, build_mgh2,
                build_mgni, build_mgco, build_mgh2ni, build_mgh2co,
                build_mg2ni, build_mg2ni_ni, build_mg2ni_co,
                build_mg2nih4, build_mg2nih4_ni, build_mg2nih4_co,
                build_mg2ni_co_nisite, build_mg2nih4_co_nisite]
    for fn in builders:
        fn()
        print(f"  wrote inputs/{fn.__name__[6:]}.in")


if __name__ == "__main__":
    main()
