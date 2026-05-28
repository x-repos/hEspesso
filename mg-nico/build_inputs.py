"""
Generate 7 internally-consistent Quantum ESPRESSO inputs in mg-nico/inputs/.

All structures use:
    - PBE + DFT-D3 dispersion
    - ecutwfc = 60 Ry, ecutrho = 480 Ry
    - PAW (Mg, H) + USPP (Ni, Co) pseudopotentials
    - vc-relax (cell + ions) for solids; relax-only for H2 in vacuum
    - cold smearing 0.01 Ry for metals; fixed occupations for the H2 molecule
    - nspin = 2 with starting_magnetization on TM atoms in TM-containing cells
    - conv_thr = 1e-9, forc_conv_thr = 1e-4 Ry/Bohr

The seven species span three reactions for the catalyst comparison:

    R1  pure         :  Mg          + H2  -> MgH2
    R2  Ni-doped Mg  :  Mg14Ni2     + 16 H2 -> Mg14Ni2H32
    R3  Co-doped Mg  :  Mg14Co2     + 16 H2 -> Mg14Co2H32

R2 and R3 share the same 2x2x2 supercell topology so that ΔH(Ni) and ΔH(Co)
can be compared on equal footing.
"""

from pathlib import Path

HERE       = Path(__file__).resolve().parent
INPUTS_DIR = HERE / "inputs"
INPUTS_DIR.mkdir(exist_ok=True)

PSEUDO_DIR = "/home/x/Workspace/espresso/pseudo"
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
# Structures
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

def substitute_two(atoms, src_species, new_species):
    """Replace first two `src_species` atoms in list with `new_species`.
    Pick the two that are farthest apart (avoid placing dopants as neighbours)."""
    src_idx = [i for i, (sp, *_) in enumerate(atoms) if sp == src_species]
    # pick first; then pick farthest from first
    i0 = src_idx[0]
    def dist(i, j):
        _, x1, y1, z1 = atoms[i]; _, x2, y2, z2 = atoms[j]
        return (x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2
    i1 = max(src_idx[1:], key=lambda j: dist(i0, j))
    out = list(atoms)
    out[i0] = (new_species, *atoms[i0][1:])
    out[i1] = (new_species, *atoms[i1][1:])
    return out

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
    nstep            = 80
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
    vdw_corr         = 'DFT-D3'
    dftd3_version    = 4
/
"""

ELECTRONS_TEMPLATE = """&ELECTRONS
    electron_maxstep = 100
    conv_thr         = 1.0d-9
    mixing_beta      = {mix}
    diagonalization  = 'david'
/
"""

IONS_TEMPLATE = "&IONS\n    ion_dynamics = 'bfgs'\n/\n"
CELL_TEMPLATE = "&CELL\n    cell_dynamics = 'bfgs'\n    press_conv_thr = 0.5\n/\n"


def write_input(name, atoms, cell_vecs, calc, kpts, magnetic_tms=()):
    """Write a complete pw.x input file to inputs/{name}.in.

    atoms: list of (species, x, y, z) in angstrom
    cell_vecs: 3 cell vectors in angstrom (rows of matrix)
    calc: 'scf' | 'relax' | 'vc-relax'
    kpts: (nk1, nk2, nk3)
    magnetic_tms: iterable of TM species names to magnetize (e.g. ('Ni',))
    """
    species = []
    for sp, *_ in atoms:
        if sp not in species:
            species.append(sp)
    ntyp = len(species)
    nat = len(atoms)

    is_h2_only = species == ["H"]
    if is_h2_only:
        occ = "fixed"
        smearing_lines = ""
        mix = 0.7
    else:
        occ = "smearing"
        smearing_lines = "smearing         = 'cold'\n    degauss          = 0.01"
        mix = 0.3 if magnetic_tms else 0.4

    if magnetic_tms:
        spin_block = ["nspin            = 2"]
        for i, sp in enumerate(species, start=1):
            mom = {"Ni": 0.5, "Co": 0.7}.get(sp, 0.0)
            spin_block.append(f"starting_magnetization({i}) = {mom}")
        spin_lines = "\n    ".join(spin_block)
    else:
        spin_lines = ""

    control  = CONTROL_TEMPLATE.format(
        calc=calc, prefix=name, pseudo_dir=PSEUDO_DIR)
    system   = SYSTEM_TEMPLATE.format(
        nat=nat, ntyp=ntyp, occ=occ,
        smearing_lines=smearing_lines, spin_lines=spin_lines)
    electrons = ELECTRONS_TEMPLATE.format(mix=mix)

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
# Build all 7 inputs
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

def build_mgni():
    atoms, cell = make_mg_supercell_2x2x2()
    atoms = substitute_two(atoms, "Mg", "Ni")
    write_input("mgni", atoms, cell, "vc-relax", (6, 6, 4), magnetic_tms=("Ni",))

def build_mgco():
    atoms, cell = make_mg_supercell_2x2x2()
    atoms = substitute_two(atoms, "Mg", "Co")
    write_input("mgco", atoms, cell, "vc-relax", (6, 6, 4), magnetic_tms=("Co",))

def build_mgh2ni():
    atoms, cell = make_mgh2_supercell_2x2x2()
    atoms = substitute_two(atoms, "Mg", "Ni")
    write_input("mgh2ni", atoms, cell, "vc-relax", (4, 4, 6), magnetic_tms=("Ni",))

def build_mgh2co():
    atoms, cell = make_mgh2_supercell_2x2x2()
    atoms = substitute_two(atoms, "Mg", "Co")
    write_input("mgh2co", atoms, cell, "vc-relax", (4, 4, 6), magnetic_tms=("Co",))


def main():
    builders = [build_mg, build_h2, build_mgh2,
                build_mgni, build_mgco, build_mgh2ni, build_mgh2co]
    for fn in builders:
        fn()
        print(f"  wrote inputs/{fn.__name__[6:]}.in")


if __name__ == "__main__":
    main()
