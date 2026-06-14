# mg-nico: Ni & Co catalyst effect on MgH2 and Mg2NiH4 thermodynamics

Self-contained study showing whether a single substitutional **Ni** or **Co**
catalyst atom destabilises the hydride of two hosts — **Mg** and **Mg2Ni** —
i.e. makes ΔH per H2 less negative so that H2 release becomes easier (the
requirement for a "better" hydrogen-storage catalyst).

## Reactions

Two hosts x (pure, +Ni, +Co), one catalyst atom per cell:

| ID | Host  | Reaction                                  | What it tests          |
|----|-------|-------------------------------------------|------------------------|
| R1 | Mg    | Mg + H2 -> MgH2                           | baseline (pure)        |
| R2 | Mg    | Mg15Ni + 16 H2 -> Mg15NiH32               | 6.25% Ni @ Mg site     |
| R3 | Mg    | Mg15Co + 16 H2 -> Mg15CoH32               | 6.25% Co @ Mg site     |
| R4 | Mg2Ni | (8/12) Mg12Ni6 + 8 H2 -> Mg8Ni4H16        | baseline (pure)        |
| R5 | Mg2Ni | defect shift: Mg11Ni7 / Mg7Ni5H16         | +1 Ni @ Mg site        |
| R6 | Mg2Ni | defect shift: Mg11CoNi6 / Mg7CoNi4H16     | +1 Co @ Mg site        |
| R6'| Mg2Ni | defect shift: Mg12CoNi5 / Mg8CoNi3H16     | +1 Co @ Ni site        |

R6' exists because Co prefers the Ni site over the Mg site in stoichiometric
Mg2Ni (E_Mg->Ni = -0.53 eV: Yoon et al., J. Magnesium Alloys 12 (2024)
4574), and the doping effect there is much weaker; both sites are computed.
For R5, excess Ni in Ni-rich Mg2Ni naturally occupies the Mg sublattice
(same reference), so the Ni @ Mg site model is the physical one.

R2/R3 share the **same 2x2x2 supercell topology** (16 metal sites, 1
substituted with TM, 32 H sites), and R5/R6 share the same Mg2Ni cells and
substitution site. Within each host the only difference between the Ni and
Co columns is the chemistry of the dopant — no supercell, k-density, or
functional artefacts confound the comparison.

For the Mg host, doped reactant and product supercells contain identical
metal content, so ΔH follows from a single balanced reaction. For Mg2Ni the
doped 18-atom metal and 28-atom hydride cells do not match in metal content,
so the doped ΔH is computed in the dilute-defect limit:

    ΔΔH = [ (E_hyd_doped - E_hyd_pure) - (E_met_doped - E_met_pure) ] / 8

where the Mg and dopant chemical potentials cancel exactly, and ΔH(doped) =
ΔH(pure) + ΔΔH.

The Mg2Ni / Mg2NiH4 geometries (18-atom hexagonal P6_222 cell; 28-atom
monoclinic LT-Mg2NiH4, C2/c, Z=4) are taken from the sibling
`mg2ni-nb2o5fe` study with its Nb substitution site restored to Mg; that
same site hosts the catalyst atom here, so dopant placement is identical
across studies. All cells are re-relaxed (vc-relax) under this study's
settings — the old study's energies are NOT reused (they carried a Hubbard
U on Ni-3d and fixed geometries, both incompatible with this framework).

## Method

All fifteen structures use identical settings:

| Setting            | Value                                        |
|--------------------|----------------------------------------------|
| Functional         | PBE (no Hubbard U anywhere)                  |
| Dispersion         | DFT-D3 (Grimme, BJ damping, v4)              |
| ecutwfc / ecutrho  | 60 Ry / 480 Ry                               |
| Pseudopotentials   | PAW for Mg/H, USPP for Ni/Co (PSlibrary 1.x) |
| Smearing           | cold, degauss = 0.01 Ry (solids)             |
| Occupations (H2)   | fixed (insulating molecule)                  |
| Spin               | per-cell ground state (see Spin treatment)   |
| conv_thr           | 1e-9                                         |
| forc_conv_thr      | 1e-4 Ry/Bohr                                 |
| Cell optimisation  | vc-relax for solids; relax-only for H2/box   |

### Spin treatment (per-cell ground state)

Each cell is computed in whichever spin state is its electronic ground
state, verified by single-point nspin=2 spin-checks (`make_spinchk.py`):

* **Ni is nonmagnetic everywhere** — pure cells, Ni-doped cells, and the
  antisite-Ni metal `mg2ni_ni` (nspin=2 gives m=0). nspin=1 used.
* **Co is magnetic** with an environment-dependent moment:
  `mgco` 1.84, `mg2ni_co` 2.83, `mg2ni_co_nisite` 1.45 µB (metals);
  `mgh2co` 1.18 µB (MgH2 hydride). These cells are relaxed with nspin=2
  and a starting moment on Co (`make_nspin2.py`).
* **Co quenches** specifically at the Mg site of Mg2NiH4 (`mg2nih4_co`
  0.37 µB) and at the Ni site of Mg2NiH4 (`mg2nih4_co_nisite` 0.00 µB) —
  octahedral H coordination kills the moment. nspin=1 is correct there.

Treating each cell at its true ground state matters: the nspin=2 correction
shifts ΔΔH(Mg:Co) from +4.70 to +4.09 and ΔΔH(Mg2Ni:Co@Mg) from +13.28 to
+15.03 (the asymmetry — magnetic metal, nonmagnetic hydride — makes the
Mg2Ni:Co case *more* destabilising). Neither shift flips a verdict.

K-points: 12x12x8 (2-atom HCP Mg), 8x8x12 (6-atom MgH2), 6x6x4 (16-atom Mg
supercell), 4x4x6 (48-atom MgH2 supercell), 8x8x4 (18-atom Mg2Ni), 8x8x6
(28-atom Mg2NiH4). Doped cells use the same grid as their pristine parent.

## Layout

```
mg-nico/
├── README.md          this file
├── build_inputs.py    generates inputs/*.in from scratch (no external deps)
├── check_inputs.py    stoichiometry/settings sanity check of inputs/*.in
├── restart_input.py   re-seed a stalled relax from its best geometry
├── make_spinchk.py    single-point nspin=2 spin-stability check of a cell
├── make_nspin2.py     re-seed a magnetic cell as nspin=2 relax
├── run_all.sh         loops pw.x over inputs/, writes outputs/ and logs/
├── enthalpy.py        parses outputs/, computes the 6 ΔH values + verdicts
├── inputs/            15 QE input files (.in)
├── outputs/           15 SCF/vc-relax outputs (.out) -- from run_all.sh
│   ├── archive-2dopant/  superseded 12.5%-dopant (2 TM atoms) run
│   └── failed/           retry trail (timeouts, nspin1 pre-correction)
├── logs/              per-job stderr
└── outdir/            QE scratch (.wfc, .save, etc.) -- safe to delete
```

## Reproduce

```bash
# (1) regenerate inputs (idempotent) and sanity-check them
python3 build_inputs.py
python3 check_inputs.py

# (2) run all 15 DFT calculations (sources ../env.sh for the GPU pw.x)
bash run_all.sh

# (3) compute the six ΔH and print the catalyst verdicts
python3 enthalpy.py
```

`run_all.sh` skips any job whose output already contains `JOB DONE`, so
re-running it is cheap. You can also select a subset:
`bash run_all.sh mgni mgh2ni` to only do the Mg:Ni pathway.

## How the comparison works

For pathway X (reactant `met_X`, product `hyd_X`, n H2 absorbed):

    ΔH_per_H2 = [ E(hyd_X) - E(met_X) - n * E(H2) ] / n     [Ry]

multiplied by 1312.75 kJ/mol/Ry (for Mg2Ni pathways via the defect-shift
formula above). The verdict is:

    ΔΔH = ΔH(dopant) - ΔH(pure)
    if ΔΔH > 0   -> dopant DESTABILISES hydride -> easier release -> BETTER
    if ΔΔH < 0   -> dopant STABILISES   hydride -> harder release -> WORSE

Experimental ΔH references (verified against primary sources):
* MgH2:    -74.5 kJ/mol H2 (Bogdanovic 1999, -74.513 at 298 K; Stampfer 1960, -17790 cal/mol = -74.4)
* Mg2NiH4: -64.4 kJ/mol H2 (Reilly & Wiswall 1968, -15.4 kcal/mol H2)

**Scope of validity.** GGA's systematic error differs per host, as
documented in the literature (all values static-lattice, kJ/mol H2):

| System  | Experiment | Published GGA                          | This work (PBE+D3) |
|---------|-----------|-----------------------------------------|--------------------|
| MgH2    | -74.5 (Bogdanovic 1999, 298 K), -74.4 (Stampfer 1960) | -63.7 PW91 (Er, PRB 79, 024105 (2009), explicitly "underestimated"); -45.3 PBE+ZPE / DMC -82 (Pozzo & Alfe, PRB 77, 104103 (2008)) | -63.75 |
| Mg2NiH4 | -64.4 (Reilly-Wiswall 1968, -15.4 kcal/mol H2) | (literature ~-63.7 to -68; not re-verified here) | -66.97 |

GGA under-binds MgH2 by 10-20 kJ/mol H2 (bare PBE is worst; benchmark
DMC+ZPE gives -82, Pozzo & Alfe, PRB 77, 104103 (2008)) yet lands within a
few kJ/mol for Mg2NiH4 -- different bonding character (ionic vs covalent
Ni-H). The same inversion appears in the large VASP-PBE survey of Yoon et
al., J. Magnesium Alloys 12 (2024) 4574 (MgH2 -52.65 vs Mg2NiH4 -57.06,
their Table 3), who state explicitly that the objective of such
calculations "was not to provide the absolute plateau pressure and enthalpy
values but rather to predict the relative change" -- and whose
relative-change predictions reproduce experimental plateau-pressure trends
for ~10 solutes in MgH2 (their Table 5) and ~14 in Mg2NiH4 (their Table 8),
including the correct magnitudes. That survey is the methodological
precedent for this study's ΔΔH framing. Raw DFT baselines therefore do NOT reproduce the experimental
cross-host ordering, and this study makes no cross-host claims from raw ΔH.
The quantitative results are (a) the within-host shifts ΔΔH, where the
host's GGA bias cancels, and (b) experiment-anchored estimates
ΔH(doped) ≈ ΔH_exp(pure host) + ΔΔH_DFT, which inherit the correct
cross-host ordering from experiment. (ZPE, not included here, destabilises
both hydrides by ~+10 kJ/mol H2: +9 for MgH2 (Pozzo & Alfe), +10.6 for
Mg2NiH4 (van Setten) -- again mostly cancelling in ΔΔH.)

## Results

Converged on RTX 5090 (PBE + DFT-D3, vc-relax, per-cell spin ground state).
ΔΔH > 0 means the dopant raises ΔH per H2 toward zero → easier H2 release →
better catalyst.

| Pathway                       | ΔH per H₂ (kJ/mol) | ΔΔH vs pure host | Verdict |
|-------------------------------|--------------------|------------------|---------|
| **R1** pure Mg                | **−63.75**         | —                | baseline |
| **R2** Mg : Ni (6.25%)        | **−51.89**         | **+11.86**       | ✅ destabilises → easier release |
| **R3** Mg : Co (6.25%)        | **−59.66**         | **+4.09**        | ✅ destabilises → easier release |
| **R4** pure Mg₂Ni             | **−66.97**         | —                | baseline |
| **R5** Mg₂Ni : Ni             | **−52.54**         | **+14.43**       | ✅ destabilises → easier release |
| **R6** Mg₂Ni : Co (Mg site)   | **−51.94**         | **+15.03**       | ✅ destabilises → easier release |
| **R6′** Mg₂Ni : Co (Ni site)  | **−68.71**         | **−1.74**        | ⚠️ slightly stabilises (≈ neutral) |

**Verdict.** Ni destabilises both hosts strongly; Co destabilises both at
the Mg site. The one exception is Co at its *thermodynamically preferred*
Ni site in Mg₂Ni (Yoon et al. 2024, E_Mg→Ni = −0.53 eV), where the effect
is essentially neutral (−1.74 — an order of magnitude below the Mg-site
shifts, comparable to the ~1.4 kJ/mol site-to-site spread in Yoon et al.,
and well below the ~10 kJ/mol ZPE correction). **Best overall catalyst: Ni** — strong on both
hosts, and Mg:Ni keeps the high 7.1 wt% H₂ capacity of Mg.

Experiment-anchored ΔH and estimated 1-bar desorption temperature
(ΔH_exp(pure host) + ΔΔH_DFT; T_des via van 't Hoff, ΔS = 130.5 J/mol/K):

| System                | ΔH anchored (kJ/mol) | T_des @1 bar | wt% H₂ |
|-----------------------|----------------------|--------------|--------|
| Mg : Ni               | −62.64               | 207 °C       | 7.08   |
| Mg₂Ni : Co (Mg site)  | −49.37               | 105 °C       | 3.36   |
| Mg₂Ni : Ni            | −49.97               | 110 °C       | 3.36   |
| Mg : Co               | −70.41               | 266 °C       | 7.08   |
| Mg pure               | −74.50               | 298 °C       | 7.66   |
| Mg₂Ni pure            | −64.40               | 220 °C       | 3.62   |
| Mg₂Ni : Co (Ni site)  | −66.14               | 234 °C       | 3.62   |

Convergence quality: all cells BFGS-converged (|F| < 5·10⁻⁴ Ry/Bohr)
**except** three dilute-defect metallic cells with sub-mRy-flat energy
landscapes — `mg2ni_ni`, `mg2nih4_ni` (and historically `mgco`, since
fixed). These oscillate under vc-relax; `mgco` converged cleanly once
switched to ions-only relax at the optimised cell, while `mg2ni_ni` and
`mg2nih4_ni` are accepted at their global-minimum energy with documented
force-unconvergence (|F| ~ 1·10⁻² and 2·10⁻³ Ry/Bohr; energy converged to
the same minimum across 4 independent attempts — see `ACCEPTED` in
`enthalpy.py`). The residual largely cancels in ΔΔH.

An earlier revision used **two** TM atoms per Mg supercell (12.5%);
archived in `outputs/archive-2dopant/`, it gave ΔΔH +14.95 (Ni) / +8.51
(Co) for MgH2 — same direction as the present 6.25% results, larger
magnitude (sublinear with concentration).

## Caveats

- **Spin: per-cell ground state (not blanket nspin=1).** Co carries a
  magnetic moment in the metals and in MgH2 but quenches at the H-coordinated
  Mg/Ni sites of Mg2NiH4; Ni is nonmagnetic throughout. See the Spin
  treatment section under Method — each cell is run at its verified ground
  state, with single-point nspin=2 spin-checks confirming the choice. (The
  earlier blanket-nspin=1 assumption was valid at 12.5% Co but wrong for an
  isolated Co atom, which retains its moment.)
- **No ZPE.** Static-lattice DFT only. ZPE typically shifts ΔH by
  +6 to +10 kJ/mol H2 *uniformly* for these systems, so the inter-pathway
  comparison is robust without it.
- **Single dopant site per sublattice in Mg2Ni.** In HCP Mg and rutile MgH2
  every Mg site is symmetry-equivalent, so 1-dopant placement is unique
  there. Mg2Ni and LT-Mg2NiH4 have inequivalent Mg sites (Mg1/Mg2/Mg3 in
  the hydride); this study uses the single (documented) substitution site
  of the mg2ni-nb2o5fe study. Yoon et al. (2024) report Co at the three
  inequivalent Mg sites of Mg2NiH4 at -54.58/-53.64/-55.02 kJ/mol H2 (Table
  7), a ~1.4 kJ/mol H2 site-to-site spread -- the scale of the uncertainty
  on our single-site Mg2Ni shifts. (Ni-site substitution is unique up to the
  two Ni Wyckoff positions; we use the first.)
- **Dilute-defect asymmetry (Mg2Ni host).** The dopant sits in an 18-atom
  metal cell (1/18 of metal sites) but a 28-atom hydride cell (1/12 of
  metal sites); the residual is a finite-size/finite-concentration error
  of the dilute approximation.
- **Cross-host dopant concentrations differ.** The Mg host carries 1
  dopant per 16 H2 (6.25% of metal sites), the Mg2Ni host 1 per 8 H2
  (8.3% of Mg sites). In the dilute regime ΔΔH per H2 scales roughly
  linearly with dopant-per-H2 content, so comparing dopant *strength*
  across hosts must account for this ~2x lever arm. Within-host
  comparisons (the study's main claims) are unaffected.
- **Doped metals are forced solid solutions.** Substitutional Ni/Co in HCP
  Mg is metastable -- experimentally Ni segregates to Mg2Ni rather than
  dissolving in Mg, and measured plateau enthalpies of catalyzed MgH2 are
  essentially unchanged from pure MgH2. The computed ΔΔH answers "what if
  the dopant sits substitutionally in the lattice", a model of intimately
  mixed dopants (e.g. ball-milled, thin-film or nanoconfined systems),
  not of phase-separated catalyst particles.
- **Bias cancellation in ΔΔH is incomplete at the dopant site.** A Ni/Co
  atom in MgH2 introduces locally the covalent TM-H bonding whose GGA
  error differs from the ionic host's; the host bias cancels in ΔΔH but
  this local contribution (one TM per cell) does not. It is bounded by
  the per-TM-H-bond error difference between the two host classes.
- **No kinetics.** ΔH speaks only to thermodynamics. A catalyst can also
  lower activation barriers; that's not addressed here. (Experimentally,
  the dominant effect of Ni/Co additives on MgH2 is kinetic.)
