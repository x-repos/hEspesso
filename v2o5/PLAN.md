# V₂O₅ Gas-Sensing DFT — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Compute adsorption energy, charge transfer, and work-function change for six gases on an α-V₂O₅(001) slab in Quantum ESPRESSO, and report whether the descriptors reproduce the experimental sensitivity order acetone > NH₃ > C₂H₄ > H₂S > CO > NO₂.

**Architecture:** A self-contained `v2o5/` sub-study mirroring `mg-nico/`: Python builders generate QE inputs (ASE for structures), a bash loop runs `pw.x`, and an analysis script parses outputs into a ranked descriptor table. One bulk `vc-relax` sets the lattice; everything that touches a gas is a (001) slab. Post-processing (`pp.x`/`projwfc.x`/`average.x`/`bader`) yields charge transfer, work function, CDD, and PDOS.

**Tech Stack:** Quantum ESPRESSO 7.5 (GPU `pw.x`, `pp.x`, `projwfc.x`, `average.x`), Python 3 + ASE + NumPy, Henkelman `bader`, bash. Design spec: `v2o5/DESIGN.md`.

## Global Constraints

Every task implicitly includes these (exact values from `v2o5/DESIGN.md`):

- Functional: rev-vdW-DF2 → `input_dft = 'vdW-DF2-b86r'` (NO `vdw_corr`/DFT-D3, NO Hubbard U at baseline).
- `ecutwfc = 60.0` Ry, `ecutrho = 600.0` Ry.
- `conv_thr = 1.0d-9`, `forc_conv_thr = 1.0d-4`, `etot_conv_thr = 1.0d-5`.
- `pseudo_dir = '/home/x/Workspace/3-hEspesso/pseudo'`. Pseudos (PBE): V `v_pbe_v1.4.uspp.F.UPF`, O `O.pbe-n-kjpaw_psl.0.1.UPF`, C `C.pbe-n-kjpaw_psl.1.0.0.UPF`, H `H.pbe-rrkjus_psl.1.0.0.UPF`, N `N.pbe-n-radius_5.UPF`, S `s_pbe_v1.4.uspp.F.UPF`.
- Spin: `nspin = 1` everywhere EXCEPT NO₂ and its complex → `nspin = 2` with a starting moment on N.
- Slabs: 2 V₂O₅ layers, (1×3) lateral cell, bottom layer frozen (`0 0 0`), ≈18 Å vacuum along c, dipole correction (`tefield`/`dipfield`, `edir = 3`).
- k-points: bulk `4 4 4`; slab `2 3 1`; gas-in-box `1 1 1` (Γ). All `0 0 0` shift.
- Always `source ../env.sh` before any `pw.x`/`pp.x`/`projwfc.x` (GPU libgomp fix).
- Gases (key → ASE-G2 name, spin): acetone→`CH3COCH3` (1), nh3→`NH3` (1), c2h4→`C2H4` (1), h2s→`H2S` (1), co→`CO` (1), no2→`NO2` (2).
- Target experimental order (decreasing sensitivity): acetone > nh3 > c2h4 > h2s > co > no2.
- Energy conversion: 1 Ry = 13.605693 eV.

## File Structure

```
v2o5/
├── DESIGN.md            spec (exists)
├── docs/                3 reference papers (exists)
├── structures/
│   └── V2O5_alpha.cif   verified bulk α-V2O5 (Task 1)
├── verify_structure.py  topology/space-group gate for the bulk (Task 1)
├── build_inputs.py      ASE → all inputs/*.in (Tasks 2–5)
├── check_inputs.py      stoichiometry/settings/pseudo validator (Task 6)
├── run_all.sh           pw.x loop, skip JOB DONE (Task 7)
├── postproc.sh          pp.x/projwfc.x/average.x/bader for winning complexes (Task 9)
├── analyze.py           parse outputs → E_ads, Δq, Δφ, ranking (Tasks 8, 9)
├── test_analyze.py      unit tests for analyze.py parsers (Task 8)
├── inputs/ outputs/ logs/ outdir/
└── README.md            results + caveats (Task 10)
```

---

### Task 1: Acquire and verify the α-V₂O₅ bulk structure

**Files:**
- Create: `v2o5/structures/V2O5_alpha.cif`
- Create: `v2o5/verify_structure.py`

**Interfaces:**
- Produces: `v2o5/structures/V2O5_alpha.cif` — an α-V₂O₅ unit cell (Pmmn, ≈14 atoms) readable by `ase.io.read`. Later tasks read it via `read("structures/V2O5_alpha.cif")`.

- [ ] **Step 1: Write the verification test**

`v2o5/verify_structure.py`:
```python
"""Gate the bulk structure BEFORE any GPU time: confirm it is α-V2O5.

Checks: orthorhombic Pmmn (#59) within symprec, lattice close to the
experimental a=11.512, b=3.564, c=4.368 A (in some axis order), 4 V + 10 O
atoms, every V 5-coordinated to O with a short vanadyl bond ~1.55-1.65 A and
the layers stacked with a van-der-Waals gap (no V-O bond > 2.2 A within a
layer). Exit nonzero on any failure.
"""
import sys
from pathlib import Path
import numpy as np
from ase.io import read
from ase.spacegroup import get_spacegroup

HERE = Path(__file__).resolve().parent
CIF = HERE / "structures" / "V2O5_alpha.cif"
EXP_ABC = sorted([11.512, 3.564, 4.368])

def main():
    atoms = read(CIF)
    errs = []

    sg = get_spacegroup(atoms, symprec=1e-2).no
    if sg != 59:
        errs.append(f"space group #{sg}, expected 59 (Pmmn)")

    abc = sorted(atoms.cell.lengths())
    if not np.allclose(abc, EXP_ABC, rtol=0.06):
        errs.append(f"lattice {abc} vs experimental {EXP_ABC} (>6%)")

    syms = atoms.get_chemical_symbols()
    nV, nO = syms.count("V"), syms.count("O")
    if (nV, nO) != (4, 10):
        errs.append(f"composition V{nV}O{nO}, expected V4O10")

    d = atoms.get_all_distances(mic=True)
    Oidx = [i for i, s in enumerate(syms) if s == "O"]
    for i, s in enumerate(syms):
        if s != "V":
            continue
        vo = sorted(d[i][j] for j in Oidx)
        nbond = sum(1 for x in vo if x < 2.2)
        if not (4 <= nbond <= 6):
            errs.append(f"V#{i} has {nbond} O within 2.2 A (expected ~5)")
        if not (1.50 < vo[0] < 1.70):
            errs.append(f"V#{i} shortest V-O {vo[0]:.3f} A, expected vanadyl ~1.58")

    if errs:
        print("STRUCTURE REJECTED:")
        for e in errs:
            print("  -", e)
        sys.exit(1)
    print(f"OK: alpha-V2O5, Pmmn(#{sg}), abc={abc}, V4O10, vanadyl bond present.")

if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Obtain the structure CIF**

Primary (Materials Project — α-V₂O₅ is entry **mp-25279**):
```bash
cd /home/x/Workspace/3-hEspesso/v2o5 && mkdir -p structures
python3 -c "from mp_api.client import MPRester; import os; \
s=MPRester(os.environ['MP_API_KEY']).get_structure_by_material_id('mp-25279'); \
s.to(filename='structures/V2O5_alpha.cif')"
```
Fallback if no `MP_API_KEY`: place any α-V₂O₅ CIF (COD/ICSD, orthorhombic Pmmn, a≈11.51, b≈3.56, c≈4.37) at `v2o5/structures/V2O5_alpha.cif`. The next step gates it either way. (If the source uses a different axis order, that is fine — `verify_structure.py` sorts the lattice lengths; the slab builder in Task 3 standardises the stacking axis.)

- [ ] **Step 3: Run verification (must pass before proceeding)**

Run: `cd /home/x/Workspace/3-hEspesso/v2o5 && python3 verify_structure.py`
Expected: `OK: alpha-V2O5, Pmmn(#59), abc=[3.56..., 4.36..., 11.51...], V4O10, vanadyl bond present.`
If REJECTED: the CIF is the wrong phase/axis — re-fetch from a cited source and rerun. Do not continue until it passes.

- [ ] **Step 4: Commit**

```bash
cd /home/x/Workspace/3-hEspesso && git add v2o5/structures/V2O5_alpha.cif v2o5/verify_structure.py && git commit -m "gas: add and verify alpha-V2O5 bulk structure"
```

---

### Task 2: Bulk vc-relax — input writer + run + lattice/gap check

**Files:**
- Create: `v2o5/build_inputs.py`
- Create: `v2o5/outputs/` `v2o5/logs/` `v2o5/outdir/` (dirs)

**Interfaces:**
- Produces: `write_qe(path, atoms, *, calc, prefix, kpts, nspin=1, magn=None, fixed_below=None, dipole=False)` — writes one pw.x input. `calc ∈ {'scf','relax','vc-relax'}`; `magn` is `{symbol: float}`; `fixed_below` freezes atoms with z < value (Å); `dipole=True` adds `tefield/dipfield edir=3`.
- Produces: `inputs/bulk.in` (vc-relax) and, after running, `outputs/bulk.out` whose final coordinates feed Task 3.

- [ ] **Step 1: Write `build_inputs.py` with the QE writer and the bulk builder**

`v2o5/build_inputs.py`:
```python
"""Generate every QE input for the V2O5 gas-sensing study into inputs/.

Settings are fixed by v2o5/DESIGN.md: rev-vdW-DF2, 60/600 Ry, no U, dipole-
corrected (001) slabs, per-system spin (NO2 -> nspin=2).
"""
from pathlib import Path
import numpy as np
from ase.io import read
from ase.build import molecule

HERE = Path(__file__).resolve().parent
PSEUDO_DIR = "/home/x/Workspace/3-hEspesso/pseudo"
INPUTS = HERE / "inputs"; INPUTS.mkdir(exist_ok=True)
(HERE / "outputs").mkdir(exist_ok=True)
(HERE / "logs").mkdir(exist_ok=True)
(HERE / "outdir").mkdir(exist_ok=True)

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


def write_qe(path, atoms, *, calc, prefix, kpts, nspin=1, magn=None,
             fixed_below=None, dipole=False):
    syms = atoms.get_chemical_symbols()
    species = list(dict.fromkeys(syms))
    nat, ntyp = len(atoms), len(species)
    L = []
    L += ["&CONTROL", f"    calculation      = '{calc}'",
          f"    prefix           = '{prefix}'",
          f"    pseudo_dir       = '{PSEUDO_DIR}'",
          "    outdir           = './outdir/'",
          "    restart_mode     = 'from_scratch'",
          "    tprnfor          = .true.", "    tstress          = .true.",
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
    L += ["    input_dft        = 'vdW-DF2-b86r'"]
    if dipole:
        L += ["    tefield          = .true.", "    dipfield         = .true.",
              "    edir             = 3", "    emaxpos          = 0.97",
              "    eopreg           = 0.05"]
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
    print(f"wrote {Path(path).name}  (nat={nat}, ntyp={ntyp}, calc={calc}, nspin={nspin})")


def build_bulk():
    atoms = read(HERE / "structures" / "V2O5_alpha.cif")
    write_qe(INPUTS / "bulk.in", atoms, calc="vc-relax", prefix="bulk", kpts=(4, 4, 4))


if __name__ == "__main__":
    build_bulk()
```

- [ ] **Step 2: Generate the bulk input and smoke-check it**

Run: `cd /home/x/Workspace/3-hEspesso/v2o5 && python3 build_inputs.py`
Then: `grep -E "input_dft|ecutrho|calculation|nat " inputs/bulk.in`
Expected: `input_dft = 'vdW-DF2-b86r'`, `ecutrho = 600.0`, `calculation = 'vc-relax'`, `nat = 14`.

- [ ] **Step 3: Run the bulk vc-relax**

Run:
```bash
cd /home/x/Workspace/3-hEspesso/v2o5 && source ../env.sh && pw.x -in inputs/bulk.in > outputs/bulk.out 2> logs/bulk.err
```
Expected: completes with `JOB DONE.` and `bfgs converged` in `outputs/bulk.out` (minutes on the 5090).

- [ ] **Step 4: Verify relaxed lattice and band gap**

Run:
```bash
cd /home/x/Workspace/3-hEspesso/v2o5 && grep -A4 "CELL_PARAMETERS" outputs/bulk.out | tail -4
grep -E "highest occupied, lowest unoccupied|Fermi" outputs/bulk.out | tail -1
```
Expected: lattice lengths within ~6% of 11.51/3.56/4.37 Å; a finite gap (`highest occupied, lowest unoccupied level` split ≈1.8–2.4 eV). If the cell collapsed or no gap appears, the input structure was wrong — return to Task 1.

- [ ] **Step 5: Commit**

```bash
cd /home/x/Workspace/3-hEspesso && git add v2o5/build_inputs.py v2o5/inputs/bulk.in v2o5/outputs/bulk.out && git commit -m "gas: bulk alpha-V2O5 vc-relax (rev-vdW-DF2)"
```

---

### Task 3: (001) slab builder + slab sanity check

**Files:**
- Modify: `v2o5/build_inputs.py` (add `relaxed_bulk()`, `build_slab()`, `slab_atoms()`)

**Interfaces:**
- Consumes: `outputs/bulk.out` (relaxed cell + coordinates).
- Produces: `slab_atoms(nlayers, nb) -> ase.Atoms` (the bare (001) slab, vacuum + frozen mask metadata via `atoms.info['fixed_below']`); `inputs/slab.in` and `inputs/slab_1L.in`.

- [ ] **Step 1: Add the relaxed-bulk reader and slab builder to `build_inputs.py`**

Add to `v2o5/build_inputs.py` (above `if __name__`):
```python
def relaxed_bulk():
    """Relaxed α-V2O5 cell from outputs/bulk.out (final coordinates)."""
    out = HERE / "outputs" / "bulk.out"
    atoms = read(out, index=-1, format="espresso-out")
    # Standardise so the SHORTEST axis-aligned lattice vector is the stacking
    # axis c (the vdW gap direction). For α-V2O5 the layers are perpendicular
    # to the ~4.37 A axis; ensure it is axis 2 (c).
    lengths = atoms.cell.lengths()
    stack = int(np.argmin([abs(l - 4.37) for l in lengths]))
    order = [i for i in range(3) if i != stack] + [stack]
    if order != [0, 1, 2]:
        atoms = atoms[:]  # copy
        atoms.set_cell(atoms.cell[order], scale_atoms=False)
        atoms.set_positions(atoms.get_positions()[:, order])
    return atoms


def slab_atoms(nlayers=2, nb=3, vacuum=18.0):
    """Bare (001) slab: bulk replicated nb along b and nlayers along c, with
    `vacuum` A added along c. α-V2O5 layers are vdW-separated and already lie
    perpendicular to c, so a (001) slab is just stacked bulk layers + vacuum."""
    bulk = relaxed_bulk()
    slab = bulk.repeat((1, nb, nlayers))
    cell = slab.get_cell()
    zmax = slab.get_positions()[:, 2].max()
    zmin = slab.get_positions()[:, 2].min()
    cell[2] = [0.0, 0.0, (zmax - zmin) + vacuum]
    slab.set_cell(cell, scale_atoms=False)
    slab.center(axis=2)
    # freeze the bottom layer (lower ~ one layer thickness)
    z = slab.get_positions()[:, 2]
    fixed_below = z.min() + (z.max() - z.min()) / (2 * nlayers)
    slab.info["fixed_below"] = float(fixed_below)
    return slab


def build_slab():
    for name, nl in (("slab", 2), ("slab_1L", 1)):
        s = slab_atoms(nlayers=nl)
        write_qe(INPUTS / f"{name}.in", s, calc="relax", prefix=name,
                 kpts=(2, 3, 1), dipole=True, fixed_below=s.info["fixed_below"])
```
Update `__main__` to also call `build_slab()`.

- [ ] **Step 2: Generate the slab inputs**

Run: `cd /home/x/Workspace/3-hEspesso/v2o5 && python3 build_inputs.py`
Then: `grep -E "nat |tefield|edir" inputs/slab.in` and `grep -c " 0 0 0$" inputs/slab.in`
Expected: `slab.in` has `nat = 84` (2 layers × (1×3) × 14), `tefield = .true.`, `edir = 3`, and a nonzero count of frozen (`0 0 0`) lines (the bottom layer). `slab_1L.in` has `nat = 42`.

- [ ] **Step 3: Sanity-check slab geometry (no QE yet)**

Run:
```bash
cd /home/x/Workspace/3-hEspesso/v2o5 && python3 -c "
from build_inputs import slab_atoms
s = slab_atoms(2); z = s.get_positions()[:,2]
assert s.get_cell()[2,2] - (z.max()-z.min()) > 15, 'vacuum too small'
assert s.get_chemical_symbols().count('V') == 24, 'expected 24 V in 2L 1x3'
print('slab OK: nat', len(s), 'vacuum', round(s.get_cell()[2,2]-(z.max()-z.min()),1))"
```
Expected: `slab OK: nat 84 vacuum ~18.0`.

- [ ] **Step 4: Run both bare slabs (1-layer and 2-layer convergence check)**

Run:
```bash
cd /home/x/Workspace/3-hEspesso/v2o5 && source ../env.sh
for n in slab slab_1L; do pw.x -in inputs/$n.in > outputs/$n.out 2> logs/$n.err; done
grep -H "^!" outputs/slab.out outputs/slab_1L.out | tail -2
```
Expected: both reach `JOB DONE.`/`bfgs converged`. Record both final energies; Task 8 confirms the surface is converged at 2 layers (top-layer geometry stable vs 1 layer). Use the 2-layer `slab` for production.

- [ ] **Step 5: Commit**

```bash
cd /home/x/Workspace/3-hEspesso && git add v2o5/build_inputs.py v2o5/inputs/slab*.in v2o5/outputs/slab*.out && git commit -m "gas: (001) 2-layer V2O5 slab + 1-vs-2 layer check"
```

---

### Task 4: Gas molecules in a box

**Files:**
- Modify: `v2o5/build_inputs.py` (add `build_molecules()`)

**Interfaces:**
- Produces: `inputs/gas_<key>.in` for each of the six gases; `outputs/gas_<key>.out` feed E_ads in Task 8.

- [ ] **Step 1: Add `build_molecules()` to `build_inputs.py`**

Add to `v2o5/build_inputs.py`:
```python
def _acetone_fallback():
    from ase import Atoms
    return Atoms("OC3H6", positions=[
        (0.000, 0.000, 1.220), (0.000, 0.000, 0.000),
        (1.287, 0.000, -0.745), (-1.287, 0.000, -0.745),
        (1.287, 0.000, -1.840),
        (2.147, 0.000, -0.340), (1.330, 0.880, -0.380),
        (-1.287, 0.000, -1.840),
        (-2.147, 0.000, -0.340), (-1.330, 0.880, -0.380)])

def build_molecules(box=15.0):
    for key, (g2, nspin) in GASES.items():
        try:
            m = molecule(g2)
        except (KeyError, NotImplementedError):
            assert key == "acetone", f"no ASE geometry for {g2}"
            m = _acetone_fallback()
        m.set_cell([box, box, box]); m.center(); m.pbc = True
        magn = {"N": 0.5} if nspin == 2 else None
        write_qe(INPUTS / f"gas_{key}.in", m, calc="relax", prefix=f"gas_{key}",
                 kpts=(1, 1, 1), nspin=nspin, magn=magn)
```
Add `build_molecules()` to `__main__`.

- [ ] **Step 2: Generate and spot-check molecule inputs**

Run: `cd /home/x/Workspace/3-hEspesso/v2o5 && python3 build_inputs.py`
Then: `grep -H "nspin" inputs/gas_no2.in inputs/gas_co.in`
Expected: `gas_no2.in: nspin = 2` (with `starting_magnetization` on N), `gas_co.in: nspin = 1`.

- [ ] **Step 3: Run the six molecules**

Run:
```bash
cd /home/x/Workspace/3-hEspesso/v2o5 && source ../env.sh
for g in acetone nh3 c2h4 h2s co no2; do pw.x -in inputs/gas_$g.in > outputs/gas_$g.out 2> logs/gas_$g.err; done
grep -L "JOB DONE" outputs/gas_*.out || echo "all gas jobs done"
```
Expected: `all gas jobs done`.

- [ ] **Step 4: Commit**

```bash
cd /home/x/Workspace/3-hEspesso && git add v2o5/build_inputs.py v2o5/inputs/gas_*.in v2o5/outputs/gas_*.out && git commit -m "gas: six gas molecules in a box (rev-vdW-DF2)"
```

---

### Task 5: Adsorption configurations (multi-site placement)

**Files:**
- Modify: `v2o5/build_inputs.py` (add `surface_sites()`, `build_adsorption()`)

**Interfaces:**
- Consumes: `slab_atoms()` from Task 3; relaxed molecule geometries from `outputs/gas_<key>.out`.
- Produces: `inputs/ads_<key>_s{0,1,2}.in` (3 sites/gas, 18 inputs); winners selected in Task 8.

- [ ] **Step 1: Add site finder + adsorption builder to `build_inputs.py`**

Add to `v2o5/build_inputs.py`:
```python
def surface_sites(slab, n=3):
    """Top-surface anchor points: the highest O (vanadyl), a neighbouring O
    (bridging), and the highest V. Returns up to n (x,y,z) Angstrom points."""
    pos = slab.get_positions(); syms = slab.get_chemical_symbols()
    ztop = pos[:, 2].max()
    top = [i for i in range(len(slab)) if pos[i, 2] > ztop - 1.6]
    Otop = sorted([i for i in top if syms[i] == "O"], key=lambda i: -pos[i, 2])
    Vtop = sorted([i for i in top if syms[i] == "V"], key=lambda i: -pos[i, 2])
    picks = []
    if Otop: picks.append(Otop[0])            # vanadyl O
    if len(Otop) > 1: picks.append(Otop[1])   # bridging/second O
    if Vtop: picks.append(Vtop[0])            # exposed V
    return [tuple(pos[i]) for i in picks[:n]]

def _relaxed_molecule(key):
    return read(HERE / "outputs" / f"gas_{key}.out", index=-1, format="espresso-out")

def build_adsorption(gap=2.2):
    slab = slab_atoms(2)
    fb = slab.info["fixed_below"]
    sites = surface_sites(slab, 3)
    for key, (_, nspin) in GASES.items():
        mol = _relaxed_molecule(key)
        for si, (sx, sy, sz) in enumerate(sites):
            m = mol.copy()
            mp = m.get_positions()
            mp -= mp.mean(axis=0)                      # centre molecule
            mp[:, 0] += sx; mp[:, 1] += sy
            mp[:, 2] += sz + gap - mp[:, 2].min()      # gap above the site
            m.set_positions(mp)
            combined = slab + m
            combined.set_cell(slab.get_cell(), scale_atoms=False)
            magn = {"N": 0.5, "V": 0.0} if nspin == 2 else None
            write_qe(INPUTS / f"ads_{key}_s{si}.in", combined, calc="relax",
                     prefix=f"ads_{key}_s{si}", kpts=(2, 3, 1), nspin=nspin,
                     magn=magn, dipole=True, fixed_below=fb)
```
Add `build_adsorption()` to `__main__`.

- [ ] **Step 2: Generate adsorption inputs**

Run: `cd /home/x/Workspace/3-hEspesso/v2o5 && python3 build_inputs.py && ls inputs/ads_*.in | wc -l`
Expected: `18` (6 gases × 3 sites). Spot-check: `grep "nat " inputs/ads_co_s0.in` → `nat = 86` (84 slab + 2 CO).

- [ ] **Step 3: Commit (inputs only; runs happen in Task 7)**

```bash
cd /home/x/Workspace/3-hEspesso && git add v2o5/build_inputs.py v2o5/inputs/ads_*.in && git commit -m "gas: multi-site adsorption inputs for six gases"
```

---

### Task 6: `check_inputs.py` validator

**Files:**
- Create: `v2o5/check_inputs.py`

**Interfaces:**
- Consumes: all `inputs/*.in`.
- Produces: a pass/fail report; exit 1 on any failure (the inputs "test suite", like `mg-nico/check_inputs.py`).

- [ ] **Step 1: Write `check_inputs.py`**

`v2o5/check_inputs.py`:
```python
"""Validate every inputs/*.in before burning GPU time.

Checks each file: rev-vdW-DF2 functional, 60/600 cutoffs, NO Hubbard U, NO
DFT-D3, pseudo files exist, nat matches the atom count, NO2 inputs are
nspin=2 and others nspin=1, slab/adsorption inputs carry the dipole field and
at least one frozen atom, and the cell is non-singular.
"""
import re, sys
from pathlib import Path
HERE = Path(__file__).resolve().parent

def check(p):
    t = p.read_text(); e = []
    if "vdW-DF2-b86r" not in t: e.append("missing rev-vdW-DF2 input_dft")
    if "ecutwfc          = 60.0" not in t: e.append("ecutwfc != 60")
    if "ecutrho          = 600.0" not in t: e.append("ecutrho != 600")
    if "HUBBARD" in t or "Hubbard" in t: e.append("has Hubbard U")
    if "vdw_corr" in t or "dftd3" in t: e.append("has DFT-D3")
    nat = int(re.search(r"nat\s*=\s*(\d+)", t).group(1))
    body = t.split("ATOMIC_POSITIONS")[1].splitlines()[1:]
    natoms = sum(1 for ln in body if len(ln.split()) >= 4
                 and ln.split()[0] in ("V", "O", "C", "H", "N", "S"))
    if nat != natoms: e.append(f"nat={nat} but {natoms} atom lines")
    nspin = int(re.search(r"nspin\s*=\s*(\d+)", t).group(1))
    is_no2 = "no2" in p.stem
    if is_no2 and nspin != 2: e.append("NO2 must be nspin=2")
    if not is_no2 and nspin != 1: e.append("non-NO2 must be nspin=1")
    if p.stem.startswith(("slab", "ads_")):
        if "tefield          = .true." not in t: e.append("slab missing dipole field")
        if " 0 0 0" not in t: e.append("slab missing frozen atoms")
    for ps in re.findall(r"^\s*\w+\s+[\d.]+\s+(\S+\.UPF)", t, re.M | re.I):
        if not (Path("/home/x/Workspace/3-hEspesso/pseudo") / ps).exists():
            e.append(f"pseudo {ps} not found")
    return e

def main():
    files = sorted((HERE / "inputs").glob("*.in"))
    bad = 0
    for p in files:
        errs = check(p)
        print(f"  [{'OK ' if not errs else 'FAIL'}] {p.name}")
        for x in errs: print("      -", x); 
        bad += bool(errs)
    if bad: sys.exit(f"\n{bad} input(s) failed.")
    print(f"\nAll {len(files)} inputs pass.")

if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Run the validator**

Run: `cd /home/x/Workspace/3-hEspesso/v2o5 && python3 check_inputs.py`
Expected: every input listed `[OK ]`, ending `All N inputs pass.` (N = bulk + 2 slabs + 6 gas + 18 ads = 27). Fix any FAIL in `build_inputs.py` and regenerate.

- [ ] **Step 3: Commit**

```bash
cd /home/x/Workspace/3-hEspesso && git add v2o5/check_inputs.py && git commit -m "gas: input validator (check_inputs.py)"
```

---

### Task 7: `run_all.sh` + run the adsorption set

**Files:**
- Create: `v2o5/run_all.sh`

**Interfaces:**
- Consumes: `inputs/*.in`. Produces: `outputs/*.out`, `logs/*.err`; skips any job already showing `JOB DONE`.

- [ ] **Step 1: Write `run_all.sh`** (adapted from `mg-nico/run_all.sh`)

`v2o5/run_all.sh`:
```bash
#!/bin/bash
# Run all (or selected) pw.x jobs from inputs/, writing outputs/ + logs/.
#   bash run_all.sh                 # all inputs/*.in
#   bash run_all.sh ads_co_s0 ...   # only the named jobs
set -e
cd "$(dirname "$0")"
source ../env.sh
mkdir -p outputs logs outdir
if [ $# -gt 0 ]; then JOBS=("$@"); else
    JOBS=(); for f in inputs/*.in; do JOBS+=("$(basename "${f%.in}")"); done
fi
echo "pw.x: $(which pw.x)"; echo "Jobs: ${JOBS[*]}"; echo "----"
for name in "${JOBS[@]}"; do
    IN="inputs/$name.in"; OUT="outputs/$name.out"; ERR="logs/$name.err"
    [ -f "$IN" ] || { echo "SKIP $name (no input)"; continue; }
    if [ -f "$OUT" ] && grep -q 'JOB DONE' "$OUT"; then
        echo "==> $name : already done, skipping"; continue; fi
    echo "==> $name : started $(date '+%H:%M:%S')"; t0=$(date +%s)
    if pw.x -in "$IN" > "$OUT" 2> "$ERR"; then
        echo "    done in $(( $(date +%s) - t0 ))s   E=$(grep '^!' "$OUT" | tail -1 | awk '{print $5}') Ry"
    else echo "    FAILED:"; tail -20 "$ERR"; exit 1; fi
done
echo "All requested jobs complete."
```

- [ ] **Step 2: Run the full adsorption set** (bulk/slab/gas already done in Tasks 2–4; this fills the 18 `ads_*`)

Run: `cd /home/x/Workspace/3-hEspesso/v2o5 && bash run_all.sh`
Expected: each `ads_*` job runs to `JOB DONE`; prior jobs report `already done, skipping`. (This is the compute-heavy step — 18 slab relaxations.)

- [ ] **Step 3: Confirm all converged**

Run: `cd /home/x/Workspace/3-hEspesso/v2o5 && grep -L "JOB DONE" outputs/ads_*.out || echo "all adsorption jobs done"`
Expected: `all adsorption jobs done`. For any non-converged job, re-run it; if BFGS oscillates, restart from its best geometry (mirror `mg-nico/restart_input.py` if needed).

- [ ] **Step 4: Commit**

```bash
cd /home/x/Workspace/3-hEspesso && git add v2o5/run_all.sh v2o5/outputs/ads_*.out && git commit -m "gas: run multi-site adsorption set"
```

---

### Task 8: `analyze.py` — E_ads + ranking (TDD)

**Files:**
- Create: `v2o5/analyze.py`, `v2o5/test_analyze.py`

**Interfaces:**
- Produces: `total_energy_ry(path) -> float|None`; `eads_ev(e_complex, e_slab, e_gas) -> float`; `main()` prints the ranked E_ads table vs the target order.
- Consumes: `outputs/slab.out`, `outputs/gas_<key>.out`, `outputs/ads_<key>_s*.out`.

- [ ] **Step 1: Write the failing parser test**

`v2o5/test_analyze.py`:
```python
from analyze import total_energy_ry, eads_ev

SAMPLE = """
     iteration #  5
!    total energy              =     -100.50000000 Ry
     convergence has been achieved
     bfgs converged in 3 scf cycles
     JOB DONE.
"""

def test_total_energy_reads_last_bang():
    p = __import__("pathlib").Path("/tmp/_t.out"); p.write_text(SAMPLE)
    assert abs(total_energy_ry(p) - (-100.5)) < 1e-9

def test_total_energy_none_if_unconverged():
    p = __import__("pathlib").Path("/tmp/_t2.out")
    p.write_text("!    total energy = -1.0 Ry\n")  # no JOB DONE
    assert total_energy_ry(p) is None

def test_eads_ev():
    # E_complex - E_slab - E_gas, Ry -> eV
    assert abs(eads_ev(-110.0, -100.0, -9.5) - (-0.5 * 13.605693)) < 1e-6

if __name__ == "__main__":
    test_total_energy_reads_last_bang()
    test_total_energy_none_if_unconverged()
    test_eads_ev()
    print("analyze tests PASS")
```

- [ ] **Step 2: Run it to confirm it fails**

Run: `cd /home/x/Workspace/3-hEspesso/v2o5 && python3 test_analyze.py`
Expected: FAIL — `ModuleNotFoundError: No module named 'analyze'` (or ImportError).

- [ ] **Step 3: Write `analyze.py`**

`v2o5/analyze.py`:
```python
"""Rank the six gases on V2O5 by adsorption energy and compare to experiment.

E_ads = E(slab+gas) - E(slab) - E(gas)   [eV; more negative = stronger]
For each gas the lowest-E_ads site among ads_<key>_s* is reported.
"""
import re, sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
OUT = HERE / "outputs"
RY_EV = 13.605693
TARGET = ["acetone", "nh3", "c2h4", "h2s", "co", "no2"]

def total_energy_ry(path):
    """Final '!' energy in Ry, or None if the run is not a converged minimum."""
    path = Path(path)
    if not path.exists(): return None
    t = path.read_text(errors="replace")
    if "JOB DONE" not in t: return None
    if "Maximum CPU time exceeded" in t: return None
    if "bfgs converged" not in t: return None   # every production run is a relax
    m = re.findall(r"!\s+total energy\s+=\s+(-?\d+\.\d+)\s+Ry", t)
    return float(m[-1]) if m else None

def eads_ev(e_complex_ry, e_slab_ry, e_gas_ry):
    return (e_complex_ry - e_slab_ry - e_gas_ry) * RY_EV

def best_site(key, e_slab):
    rows = []
    for p in sorted(OUT.glob(f"ads_{key}_s*.out")):
        ec = total_energy_ry(p)
        eg = total_energy_ry(OUT / f"gas_{key}.out")
        if ec is None or eg is None or e_slab is None: continue
        rows.append((eads_ev(ec, e_slab, eg), p.stem))
    return min(rows) if rows else (None, None)

def main():
    e_slab = total_energy_ry(OUT / "slab.out")
    if e_slab is None: sys.exit("slab.out missing/unconverged")
    print(f"{'gas':10s} {'E_ads (eV)':>12s}  {'site':>14s}")
    print("-" * 40)
    results = {}
    for key in TARGET:
        ea, site = best_site(key, e_slab)
        results[key] = ea
        print(f"{key:10s} {ea:12.3f}  {site:>14s}" if ea is not None
              else f"{key:10s} {'(missing)':>12s}")
    ranked = [k for k, v in sorted(results.items(), key=lambda kv: kv[1])
              if v is not None]  # most negative first
    print("\nDFT ranking (strongest binding first):", " > ".join(ranked))
    print("Experimental order (most sensitive first):", " > ".join(TARGET))
    print("MATCH" if ranked == TARGET else "MISMATCH — see DESIGN.md escalation")

if __name__ == "__main__":
    main()
```

- [ ] **Step 4: Run the test to confirm it passes**

Run: `cd /home/x/Workspace/3-hEspesso/v2o5 && python3 test_analyze.py`
Expected: `analyze tests PASS`.

- [ ] **Step 5: Run the real ranking**

Run: `cd /home/x/Workspace/3-hEspesso/v2o5 && python3 analyze.py`
Expected: a table of E_ads per gas, a DFT ranking line, and `MATCH`/`MISMATCH`. (A MISMATCH is a valid scientific result — it triggers Task 9 descriptors and the DESIGN.md escalation, not a code failure.)

- [ ] **Step 6: Commit**

```bash
cd /home/x/Workspace/3-hEspesso && git add v2o5/analyze.py v2o5/test_analyze.py && git commit -m "gas: E_ads ranking vs experimental order (analyze.py)"
```

---

### Task 9: `postproc.sh` + Bader, work function, CDD, PDOS

**Files:**
- Create: `v2o5/postproc.sh`
- Modify: `v2o5/analyze.py` (add `bader_transfer()`, `work_function()`, and a descriptor table)

**Interfaces:**
- Consumes: winning `ads_<key>_s*` prefixes (from Task 8) + their `outdir/<prefix>.save`.
- Produces: per-complex `ACF.dat` (Bader), `avg_<prefix>.dat` (planar potential), CDD cube; descriptor table appended by `analyze.py`.

- [ ] **Step 1: Write `postproc.sh`** (runs on the lowest-E_ads site of each gas; pass the prefixes)

`v2o5/postproc.sh`:
```bash
#!/bin/bash
# Post-process winning adsorption complexes: Bader charge, work function,
# charge-density difference. Usage: bash postproc.sh ads_co_s0 ads_nh3_s1 ...
set -e; cd "$(dirname "$0")"; source ../env.sh
mkdir -p pp
for pre in "$@"; do
  # all-electron valence density for Bader
  cat > pp/${pre}_rho.in <<EOF
&INPUTPP
  prefix='${pre}', outdir='./outdir/', plot_num=0, filplot='pp/${pre}.rho'
/
&PLOT  iflag=3, output_format=6, fileout='pp/${pre}.cube' /
EOF
  pp.x -in pp/${pre}_rho.in > logs/${pre}_pp.out 2>&1
  ( cd pp && bader ${pre}.cube > ${pre}_bader.log 2>&1 && mv ACF.dat ${pre}_ACF.dat )
  # planar-averaged electrostatic potential -> work function
  cat > pp/${pre}_pot.in <<EOF
&INPUTPP
  prefix='${pre}', outdir='./outdir/', plot_num=11, filplot='pp/${pre}.pot'
/
&PLOT iflag=3, output_format=6, fileout='pp/${pre}.potcube' /
EOF
  pp.x -in pp/${pre}_pot.in > logs/${pre}_pot.out 2>&1
  echo "1
pp/${pre}.pot
1.0
2000
3
3.0" | average.x > pp/avg_${pre}.dat 2>> logs/${pre}_pot.out
  echo "post-processed $pre"
done
```
(Also run the same `plot_num=0` density for the bare `slab` and each `gas_<key>` to build CDD cubes; the script can be invoked with those prefixes too.)

- [ ] **Step 2: Add descriptor parsers to `analyze.py`**

Add to `v2o5/analyze.py`:
```python
def work_function(prefix):
    """phi = V_vacuum - E_Fermi (eV). Reads pp/avg_<prefix>.dat (col 3 = planar
    avg potential in Ry vs z) and the run's Fermi energy."""
    avg = HERE / "pp" / f"avg_{prefix}.dat"
    out = OUT / f"{prefix}.out"
    if not (avg.exists() and out.exists()): return None
    pot = [float(l.split()[2]) for l in avg.read_text().splitlines() if len(l.split()) >= 3]
    vvac_ry = max(pot)  # plateau in the vacuum region
    ef = re.findall(r"the Fermi energy is\s+(-?\d+\.\d+)\s+ev", out.read_text())
    if not ef: return None
    return vvac_ry * RY_EV - float(ef[-1])

def bader_transfer(prefix, gas_natoms):
    """Net charge on the adsorbate (e): sum(ZVAL - Bader charge) over the last
    `gas_natoms` atoms in pp/<prefix>_ACF.dat. Positive = molecule donates."""
    acf = HERE / "pp" / f"{prefix}_ACF.dat"
    if not acf.exists(): return None
    rows = [l.split() for l in acf.read_text().splitlines() if l.split() and l.split()[0].isdigit()]
    charges = [float(r[4]) for r in rows][-gas_natoms:]   # Bader populations
    # ZVAL per atom must be read from the .out (pseudo valence); see note below.
    return charges  # finalised in Step 4 against ZVAL from outputs/<prefix>.out
```
Note: in Step 4 the ZVAL per species is parsed from the QE output header (`valence` column) so the net transfer = Σ(ZVAL − Bader population). Keep `bader_transfer` returning the populations until ZVAL wiring is added, so the table can be assembled incrementally.

- [ ] **Step 3: Run post-processing on the winning complexes**

Run (substitute the actual lowest-E_ads prefixes printed by Task 8):
```bash
cd /home/x/Workspace/3-hEspesso/v2o5 && source ../env.sh
bash postproc.sh slab ads_acetone_s0 ads_nh3_s0 ads_c2h4_s0 ads_h2s_s0 ads_co_s0 ads_no2_s0
ls pp/*_ACF.dat pp/avg_*.dat
```
Expected: an `ACF.dat` and an `avg_*.dat` per prefix.

- [ ] **Step 4: Wire ZVAL and print the descriptor table**

Add a `descriptor_table()` to `analyze.py` that, for each gas's winning complex, prints E_ads (Task 8), Bader net transfer Δq (Σ ZVAL − population, ZVAL read from `outputs/<prefix>.out` `valence` per species), and Δφ = `work_function(complex) − work_function('slab')`; then prints the ranking by each descriptor against `TARGET`.
Run: `cd /home/x/Workspace/3-hEspesso/v2o5 && python3 analyze.py`
Expected: a three-column table (E_ads, Δq, Δφ) with a ranking line per descriptor and a MATCH/MISMATCH verdict for each.

- [ ] **Step 5: Commit**

```bash
cd /home/x/Workspace/3-hEspesso && git add v2o5/postproc.sh v2o5/analyze.py && git commit -m "gas: Bader charge transfer + work-function descriptors"
```

---

### Task 10: Results README

**Files:**
- Create: `v2o5/README.md`

**Interfaces:** Consumes the final `analyze.py` table and convergence notes. Produces the study writeup.

- [ ] **Step 1: Write `v2o5/README.md`** in the `mg-nico` style: question; method table (rev-vdW-DF2, 60/600, (001) 2-layer slab, dipole, per-system spin); the descriptor table (E_ads, Δq, Δφ per gas); the DFT-vs-experiment ranking and an honest verdict (which descriptor, if any, reproduces acetone > NH₃ > C₂H₄ > H₂S > CO > NO₂); CDD/PDOS figure callouts; and the caveats from `DESIGN.md` §11 (weak-physisorber near-degeneracy, NO₂ anomaly, single facet, static lattice, no kinetics/humidity). If the order mismatched, document where and which escalation (more sites / U on V / convergence) was applied — `mg-nico/NOTE.md` style.

- [ ] **Step 2: Verify it reflects real numbers**

Run: `cd /home/x/Workspace/3-hEspesso/v2o5 && python3 analyze.py` and confirm every number in `README.md` matches the script output (no hand-edited values).

- [ ] **Step 3: Commit**

```bash
cd /home/x/Workspace/3-hEspesso && git add v2o5/README.md && git commit -m "gas: V2O5 gas-sensing results README"
```

---

## Notes for the implementer

- **MISMATCH is a result, not a bug.** If `analyze.py` prints MISMATCH, follow `DESIGN.md` §10 escalation in order (lean on Δφ/Δq → add sites/orientations → add `&SYSTEM lda_plus_u=.true., Hubbard_U(V)=3.5` only where PDOS shows V⁵⁺→V⁴⁺ → denser k / 3-layer slab for near-degenerate pairs). Re-run only the affected jobs; `run_all.sh` skips the rest.
- **Spin checks.** Before trusting NO₂ and the strongest binder, confirm the spin ground state with a single-point `nspin=2` re-run at the relaxed geometry (mirror `mg-nico/make_spinchk.py`); if the moment quenches to ~0, the `nspin=1` energy stands.
- **Bare-slab convergence.** Compare top-layer V=O bond and surface energy between `slab` (2L) and `slab_1L` (1L); if they differ by < ~0.02 eV/Å² and < ~0.01 Å, 2 layers is converged (it is the production model regardless).
- **`bader` binary.** Task 9 needs the Henkelman `bader` executable on PATH; if absent, fall back to Löwdin charges via `projwfc.x` (`lowdin`-summed per-atom populations) and note the method in the README.
```
