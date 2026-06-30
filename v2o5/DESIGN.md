# Design spec — V₂O₅ gas-sensing DFT study

_Date: 2026-06-28. Status: design, awaiting review._

## 1. Goal

Compute, in Quantum ESPRESSO, the DFT descriptors of six gases adsorbed on an
α-V₂O₅ surface and show they reproduce the **experimentally measured sensitivity
ranking** (decreasing):

    acetone  >  NH₃  >  C₂H₄  >  H₂S  >  CO  >  NO₂

The DFT plays the same supporting/explanatory role it does in the three
reference papers in `v2o5/docs/` (CuFe₂O₄, SnS₂, CoFe₂O₄): explain the measured
selectivity through **adsorption energy + charge transfer (+ work function)**,
not predict absolute sensor response.

**Success criterion.** A descriptor — primarily |E_ads|, supported by Bader
charge transfer and work-function change Δφ — ranks the six gases in the
experimental order, or as close as DFT credibly allows, with honest
documentation of any disagreement.

## 2. What we are matching

- **Experiment:** a V₂O₅-based sensor study (Vietnamese; to be supplied by the
  user) measuring the six-gas sensitivity order above. The *order* is the firm
  input; exact response values and the measured morphology/operating
  temperature will refine the surface choice and the "match" target if the
  paper is provided.
- **Genre:** the three `v2o5/docs/` papers are each experiment + supporting DFT
  on a surface slab/monolayer (E_ads, charge transfer, CDD, DOS, sometimes work
  function). Two used VASP with non-local vdW-DF functionals (rev-vdW-DF2,
  optPBE-vdW); one used plane-wave PBE. This study reproduces that method in QE
  — the hEspesso house code — for V₂O₅.

## 3. Scope

**In scope:** six gas molecules; one V₂O₅ facet ((001)); static-lattice DFT
adsorption thermodynamics and electronic descriptors.

**Out of scope (documented as caveats, not modelled):** reaction kinetics and
response/recovery dynamics, operating-temperature effects, humidity, multiple
facets / grain boundaries, oxygen-vacancy chemistry (unless a reduction trigger
forces it), and coverage effects beyond the dilute single-molecule limit.

## 4. Material and surface model

- **Phase:** orthorhombic α-V₂O₅ (space group *Pmmn*, No. 59), experimental
  cell a≈11.512, b≈3.564, c≈4.368 Å, Z=2 (V₄O₁₀, 14 atoms/cell). Layered,
  stacked along **c**, van-der-Waals-bound between layers.
- **Surface:** the **(001)** basal cleavage plane — the dominant exposed face in
  V₂O₅ sensors — which presents the three chemically distinct oxygens (terminal
  **vanadyl** O=V, **bridging** O, **chain** O) and exposed V⁵⁺ sites.
- **Slab thickness:** **2 V₂O₅ layers**; bottom layer frozen at bulk geometry,
  top layer + adsorbate relaxed. Validated against a **1-layer (monolayer)**
  model (α-V₂O₅ layers are vdW-decoupled, so a monolayer is a standard, legitimate
  V₂O₅(001) model); if the 1-vs-2-layer check agrees, the monolayer is used for
  the production runs to cut cost.
- **Lateral cell:** **(1×3)** supercell — a≈11.5 Å is already large; b is short
  (≈3.56 Å) so it is tripled to ≈10.7 Å — giving ≥~10 Å between an adsorbate and
  its periodic images (dilute limit). Resulting size ≈42 atoms (1 layer) to ≈84
  atoms (2 layers), plus the molecule.
- **Vacuum:** ≈18–20 Å along c, with a **dipole correction** along z. One-sided
  adsorption creates a net slab dipole; the correction is required for accurate
  E_ads *and* for a meaningful work function.

## 5. DFT method

| Setting | Value | Rationale |
|---|---|---|
| Code | `pw.x` (GPU build via `../env.sh`), QE 7.5 | repo standard |
| Functional | **rev-vdW-DF2**, `input_dft='vdW-DF2-b86r'` | adsorption of the weak physisorbers is dispersion-controlled, where vdW-DF ≫ GGA+D3; also matches the SnS₂ reference paper → comparable to the genre |
| Hubbard U | **none at baseline**; conditional U≈3–4 eV on V-3d | V₂O₅ is formally V⁵⁺ (d⁰) → U is irrelevant to the clean surface; add it only where PDOS shows a gas reducing V⁵⁺→V⁴⁺ |
| Spin | per-system ground state; nspin=2 for NO₂ and any V⁴⁺-reduced complex | NO₂ is an open-shell doublet; a reduced V⁴⁺ carries a moment. Single-point nspin=2 spin-checks confirm each choice (mg-nico practice) |
| Pseudopotentials | V `v_pbe_v1.4.uspp.F.UPF`, O `O.pbe-n-kjpaw_psl.0.1.UPF`, C `C.pbe-n-kjpaw_psl.1.0.0.UPF`, H `H.pbe-rrkjus_psl.1.0.0.UPF`, N `N.pbe-n-radius_5.UPF`, S `s_pbe_v1.4.uspp.F.UPF` | all in `pseudo/`; PBE pseudos are correct for vdW-DF (non-local correlation added on top) |
| Cutoffs | ecutwfc=60 Ry, ecutrho=600 Ry | 600 (10×) for the ultrasoft V/S pseudos |
| Convergence | conv_thr=1e-9, forc_conv_thr=1e-4 Ry/Bohr, etot_conv_thr=1e-5 | repo standard |
| k-points | slab ≈2×3×1 Γ-centered (converged); gas box Γ only | large cell → sparse k |
| Smearing | cold, degauss=0.01 Ry (slab); molecules insulating in box (NO₂ smeared + spin-polarized) | |

**Geometry protocol:** (1) `vc-relax` bulk α-V₂O₅ under rev-vdW-DF2; (2) cut the
(001) slab from the relaxed cell; (3) relax ions only (fixed cell, bottom layer
frozen); molecules relaxed separately in a box.

## 6. Gas references and adsorption configurations

- **Gas geometries** from ASE's G2 set — `CH3COCH3` (acetone), `NH3`, `C2H4`,
  `H2S`, `CO`, `NO2` — each relaxed in a ≈15 Å cubic box with the same
  functional/cutoffs and the correct spin (NO₂ doublet; others singlet).
- **Adsorption:** for each gas, **3–4 starting configurations** (molecule over
  the vanadyl O, over a bridging O, and over a V site; 1–2 orientations), relax
  each, and report the **lowest-energy** complex.

## 7. Descriptors (for the winning complex of each gas)

1. **E_ads = E(slab+gas) − E(slab) − E(gas)** [eV]. Primary ranking quantity;
   more negative ⇒ stronger binding ⇒ (hypothesis) higher sensitivity.
2. **Bader charge transfer Δq** (`pp.x` all-electron density → Henkelman `bader`).
   Sign distinguishes donor (reducing) vs acceptor (oxidizing). Löwdin charges
   (`projwfc.x`) as a cross-check.
3. **Work-function change Δφ = φ(slab+gas) − φ(slab)** (`pp.x` plot_num=11
   electrostatic potential → `average.x` planar average; φ = E_vac − E_F). The
   measured resistance change is driven by exactly this band-bending/charge
   transfer, so Δφ is arguably the best sensitivity proxy.
4. **CDD** Δρ = ρ(slab+gas) − ρ(slab) − ρ(gas) and **PDOS** (`projwfc.x`):
   the explanatory figures, and the check for V⁵⁺→V⁴⁺ reduction (the U trigger).

## 8. Workflow / architecture (mirrors `mg-nico/`)

```
v2o5/
├── docs/            the 3 reference papers (already present)
├── DESIGN.md        this spec
├── build_inputs.py  ASE: relaxed α-V2O5 → cut (001) slab → place each molecule
│                    (multi-config) → emit inputs/*.in
├── check_inputs.py  stoichiometry / settings / pseudo-existence sanity check
├── run_all.sh       loop pw.x (sources ../env.sh), skip JOB DONE, subset selection
├── analyze.py       parse outputs → E_ads, Bader Δq, Δφ; rank gases; compare to
│                    the experimental order
├── postproc.sh      pp.x / projwfc.x / average.x / bader for the winning
│                    complexes (CDD, PDOS, Δφ)
├── inputs/ outputs/ logs/ outdir/
└── README.md        written at the end: results + caveats, mg-nico style
```

Structure building uses ASE (as in `mxene/build_interface.py`); the slab
freeze + dipole-correction conventions reuse the `mxene` patterns.

**Job budget:** 1 bulk vc-relax + 1–2 bare slabs + 6 gas-in-box + ~18
adsorption relaxations (6 gases × ~3 configs) + post-processing ≈ **25–30 jobs**,
each a small-to-medium slab — comfortable on the RTX 5090.

## 9. Validation plan

1. **Bulk:** lattice constants within a few % of experiment; band gap ≈2.0–2.3 eV
   (rev-vdW-DF2 typically ≈2 eV) — a sanity check, not an exactness claim.
2. **1-vs-2 layer** bare-slab check: top-layer geometry / surface energy
   converged → fixes the production slab thickness.
3. **k-point + vacuum** spot check on one adsorption system (e.g. NH₃).
4. **Spin checks** on NO₂ and the most strongly bound complex.

## 10. Matching strategy and escalation

Start from the E_ads ranking. If it does not reproduce the experimental order:

1. Lean on **Δφ** and **Bader Δq** (the resistance-relevant descriptors) and
   report which descriptor tracks the experiment.
2. Add more **sites/orientations** for the mis-ranked gas.
3. Add **U on V-3d** where reduction is detected.
4. Escalate **convergence** (denser k, thicker slab) for near-degenerate pairs.

Document honestly where DFT matches and where it does not — in the
`mg-nico/NOTE.md` style.

## 11. Risks / caveats

- The weak physisorbers (CO, C₂H₄, H₂S) have small, near-degenerate E_ads, so
  the ranking *among* them is delicate.
- **NO₂ ranked least sensitive is unusual** — NO₂ normally binds oxides strongly.
  This is the sharpest test of the model and may signal that the experimental
  ordering is governed by kinetics / operating temperature rather than by
  adsorption thermodynamics alone.
- Single facet (001), dilute coverage, static lattice, no humidity, no kinetics —
  all genre-standard simplifications, to be stated plainly.
- The vdW-DF2 + U combination, if triggered, needs careful U calibration.

## 12. Deliverables

- A ranked descriptor table (E_ads, Δq, Δφ per gas) against the experimental order.
- CDD + PDOS figures per gas.
- `v2o5/README.md` with results and honest caveats.
