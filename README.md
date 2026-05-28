# hEspesso — DFT studies of Mg-based hydrogen-storage catalysts

Quantum ESPRESSO (PBE + DFT-D3) workflows that test whether transition-metal
and oxide additives **destabilise** Mg-based hydrides — i.e. raise ΔH per H₂
toward zero so that H₂ release becomes easier under practical conditions.

![Mg2NiH4 supercells: pristine, Nb-doped, and Nb,Fe-doped](mg2ni-nb2o5fe/tex/mg2nih4.png)

*Supercells studied in the Mg₂NiH₄ branch: (a) pristine, (b) Nb-substituted,
(c) Nb,Fe co-substituted. One Mg site at the cell edge is replaced by Nb
(blue) and the neighbouring Mg by Fe (red-brown).*

## Sub-studies

| Directory | Question |
|-----------|----------|
| [`mg-nico/`](mg-nico/README.md) | Do Ni and Co substitutional dopants in MgH₂ lower ΔH per H₂? Baseline + 12.5% TM dopant at fixed 2×2×2 supercell topology. |
| [`mg2ni-nb2o5fe/`](mg2ni-nb2o5fe/README.md) | Do Nb (from Nb₂O₅) and Fe additives destabilise Mg₂NiH₄? Pristine vs. Nb-doped vs. Nb,Fe co-doped. |
| [`pseudo/`](pseudo/) | Shared UPF pseudopotential library (PSlibrary 1.x, ONCV, GBRV) used by both studies. |

Each sub-study is self-contained: `build_inputs.py` → `run_all.sh` →
`enthalpy.py`. See the per-directory README for the reactions, settings,
and reproduction steps.

## Method (shared across studies)

| Setting            | Value                                            |
|--------------------|--------------------------------------------------|
| Code               | Quantum ESPRESSO (`pw.x`, GPU build)             |
| Functional         | PBE                                              |
| Dispersion         | DFT-D3 (Grimme, BJ damping, v4)                  |
| ecutwfc / ecutrho  | 60 Ry / 480 Ry                                   |
| Pseudopotentials   | PAW for Mg/H, USPP/ONCV for TMs (PSlibrary 1.x)  |
| Smearing           | cold, degauss = 0.01 Ry (metallic cells)         |
| Spin               | nspin = 2 with magnetic init for TM-containing cells |
| conv_thr           | 1e-9                                             |
| forc_conv_thr      | 1e-4 Ry/Bohr                                     |
| k-density          | ~0.03 Å⁻¹                                        |

## Repository layout

```
hEspesso/
├── README.md             this file
├── mg-nico/              MgH2 + Ni/Co dopant study
├── mg2ni-nb2o5fe/        Mg2NiH4 + Nb2O5/Fe catalyst study
│   └── tex/              LaTeX report + figures
└── pseudo/               shared UPF pseudopotential library
```

## Reproduce

Each sub-study has its own `run_all.sh` / `enthalpy.py`. From a fresh checkout:

```bash
# pick a study
cd mg2ni-nb2o5fe        # or: cd mg-nico

python3 build_inputs.py # (re)generate QE input files
bash run_all.sh         # run pw.x over inputs/
python3 enthalpy.py     # parse outputs/, compute ΔH and verdict
```

`run_all.sh` skips any job whose output already contains `JOB DONE`, so
re-running it after a partial run is cheap.

## Caveats

- **Static-lattice DFT.** No ZPE, no finite-T phonon contributions. ZPE
  typically shifts ΔH by +6 to +10 kJ/mol H₂ roughly uniformly for these
  systems, so cross-pathway comparisons are robust without it.
- **Single dopant configuration.** TM atoms are placed to minimise
  dopant-dopant interaction; a full study would average over
  symmetry-inequivalent placements.
- **Thermodynamics only.** ΔH does not capture kinetic barriers.
