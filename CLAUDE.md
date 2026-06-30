# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

`hEspesso` is a collection of **Quantum ESPRESSO (`pw.x`) DFT studies** testing
whether transition-metal/oxide additives *destabilise* Mg-based hydrides — i.e.
raise ΔH per H₂ toward zero so H₂ release gets easier. It is a research/compute
repository, not an application: the "code" is Python that **generates** QE input
files and **parses** QE output files to compute reaction enthalpies. There is no
build system, package, or test suite.

Each top-level directory (except `pseudo/`) is a **self-contained sub-study**
with its own README, its own settings, and its own `inputs/`→`outputs/` data.
Read the sub-study's README before working in it.

| Directory | Status | What |
|-----------|--------|------|
| `mg-nico/` | complete, canonical | Single Ni/Co dopant on two hosts (Mg, Mg₂Ni). Best-maintained; copy its conventions. |
| `mg2ni-nb2o5fe/` | older, **frozen** | Nb₂O₅+Fe on Mg₂NiH₄. Different settings (see below). Has a reproducibility gap (`mgh2-cif/` refs are gone). |
| `mxene/` | in progress | Mg₂Ni(H₄) on Ti₃C₂ MXene slabs. Uses ASE to build heterostructures; slab/surface settings differ. |
| `v2o5/` | active | V₂O₅(001) gas-sensing DFT — 6 gases, rev-vdW-DF2. See `v2o5/DESIGN.md` + `PLAN.md`. |
| `pseudo/` | shared | UPF pseudopotential library used by all studies (absolute-path referenced). |

## Running calculations

**Always `source env.sh` (or let `run_all.sh` do it) before invoking `pw.x`.**
`env.sh` puts the NVHPC compiler libs *first* on `LD_LIBRARY_PATH` so
`libgomp.so.1` resolves to NVHPC's shim — without it the GPU `pw.x` (QE 7.5,
RTX 5090 build at `/home/x/Programs/espresso_gpu/bin`) aborts at startup with
`libgomp: TODO`. Do not reorder that variable.

The canonical per-study pipeline (run from inside the study directory, e.g. `mg-nico/`):

```bash
python3 build_inputs.py     # (re)generate inputs/*.in from scratch — idempotent, no external deps
python3 check_inputs.py     # sanity-check stoichiometry/settings before burning GPU time
bash run_all.sh             # loop pw.x over inputs/, writing outputs/*.out and logs/*.err
python3 enthalpy.py         # parse outputs/, compute ΔH / ΔΔH and print verdicts
```

`run_all.sh` **skips any job whose `outputs/*.out` already contains `JOB DONE`**,
so re-running is cheap. Run a subset by name (the `inputs/` basename without
`.in`): `bash run_all.sh mgni mgh2ni`. The `DEFAULT=(...)` array in `run_all.sh`
lists every job and its run order (cheap jobs first to fail fast).

Figures and report (where present, e.g. `mg-nico/`):

```bash
python3 render_crystals.py      # ASE → ball-and-stick PNGs into figures/ (needs ase, matplotlib)
cd tex && ./compile.sh report.tex   # 3× pdflatex + bibtex, then cleans aux files
```

`mxene/` uses ASE-based builders (`build_interface.py`, `build_hetero.py`)
instead of a single `build_inputs.py`, and has no `enthalpy.py`/`run_all.sh`
yet — drive `pw.x` manually there for now.

## The cross-study comparability rule (most important correctness constraint)

A DFT total energy is only physically meaningful as a **difference within one
consistent set of settings**. The studies deliberately use *different* settings,
so **energies (`.out` total energies, `.pwo`, ΔH) must never be mixed across
sub-studies.** The same compound lands on different numbers in each:

| | `mg-nico` | `mg2ni-nb2o5fe` |
|---|---|---|
| Hubbard U | **none** | U=9.0 eV Ni-3d, 4.0 eV Fe-3d |
| Geometry | `vc-relax` (cell+ions) | `scf` (frozen cell) |
| Dispersion | DFT-D3 **BJ** (`dftd3_version=4`) | DFT-D3 zero-damping (v3) |
| Spin | per-cell ground state | mixed nspin=1/2 |

`mg2ni-nb2o5fe/NOTE.md` documents this in full. The `mg-nico` study even
re-relaxes the Mg₂Ni geometries it borrows from `mg2ni-nb2o5fe` rather than
reuse its energies. When asked to "compare" or "combine" results, confirm they
came from the same study's settings first.

## How ΔH / ΔΔH work (the scientific core)

For a hydrogenation pathway (reactant metal cell, product hydride cell, n H₂):

```
ΔH_per_H2 = [ E(hydride) − E(metal) − n·E(H2) ] / n   (Ry) × 1312.75 kJ/mol/Ry
```

When doped reactant/product cells have *mismatched* metal content (the Mg₂Ni
host), the doped ΔH is computed in the **dilute-defect limit** so the Mg/dopant
chemical potentials cancel:

```
ΔΔH = [ (E_hyd_doped − E_hyd_pure) − (E_met_doped − E_met_pure) ] / n_H2
```

Verdict convention: **ΔΔH > 0 ⇒ dopant destabilises hydride ⇒ easier release ⇒
"better" catalyst.** Raw cross-host DFT baselines carry per-host GGA bias, so
cross-host claims use experiment-anchored ΔH (`ΔH_exp(pure) + ΔΔH_DFT`), never
raw DFT — see `mg-nico/enthalpy.py` and the README "Scope of validity".

## Conventions baked into the generators/parsers

- **`nosym=.true.`** is set on every doped/substituted cell (the dopant breaks
  parent symmetry and `vc-relax` otherwise trips on "not orthogonal operation").
  Pristine cells keep symmetry. `check_inputs.py` asserts this per job.
- **Spin is per-cell, not blanket.** In `mg-nico`, `SPIN2_CELLS` lists the cells
  run `nspin=2` with a starting moment on the TM (Co is magnetic in metals and
  MgH₂ but quenches at H-coordinated sites in Mg₂NiH₄; Ni is nonmagnetic
  everywhere). `make_spinchk.py` generates a single-point `nspin=2` check to
  confirm a chosen ground state; `make_nspin2.py` re-seeds a cell as an
  `nspin=2` relax.
- **Energy parsing rejects false "JOB DONE".** `enthalpy.py:total_energy()`
  returns `None` (with a warning) if a run hit `max_seconds`, never printed
  `bfgs converged`, or lacks `JOB DONE` — so an unconverged-geometry energy is
  never silently used. Genuinely force-unconverged-but-energy-converged cells
  are whitelisted explicitly in the `ACCEPTED` dict with a documented reason and
  hard-coded global-minimum energy.
- **Restarting a stalled relax:** `restart_input.py <name>` re-seeds
  `inputs/<name>.in` from the *lowest-energy* geometry in its trajectory and
  archives the stale output to `outputs/failed/`, so `run_all.sh` reruns it.
- Pseudopotentials are referenced by **absolute path** to `pseudo/`
  (`PSEUDO_DIR`/`PSEUDO` constant in each builder). Mg/H use PAW, Ni/Co/Ti use
  USPP. `ecutwfc=60 Ry`; `ecutrho` is 480 Ry (bulk) or 600 Ry (mxene slabs).

## Data hygiene / git

`.gitignore` excludes QE scratch (`outdir/`, `*.wfc*`, `*.save/`, `CRASH`),
`logs/`, and LaTeX aux files. **`outputs/*.out` and `inputs/*.in` ARE committed**
— they are the study's data of record. `outdir/` is pure scratch and safe to
delete. Note `*.out` is intentionally *not* gitignored despite being a LaTeX
aux extension, because QE output files use it.
