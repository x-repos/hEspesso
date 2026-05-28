# mg-nico: Ni & Co catalyst effect on MgH2 thermodynamics

Self-contained study showing whether substitutional **Ni** and **Co** dopants
in magnesium destabilise the resulting magnesium hydride (i.e. lower the
ΔH per H2 of formation, making H2 release easier — the requirement for a
"better" hydrogen-storage catalyst).

## Reactions

| ID | Reaction                                          | What it tests        |
|----|---------------------------------------------------|----------------------|
| R1 | Mg          +    H2  ->  MgH2                     | baseline (pure)      |
| R2 | Mg14Ni2     + 16 H2  ->  Mg14Ni2H32               | 12.5% Ni dopant      |
| R3 | Mg14Co2     + 16 H2  ->  Mg14Co2H32               | 12.5% Co dopant      |

R2 and R3 share the **same 2x2x2 supercell topology** (16 metal sites,
2 substituted with TM, 32 H sites). That guarantees the only difference
between the Ni and Co columns is the chemistry of the dopant — no
supercell, k-density, or Hubbard artefacts confound the comparison.

## Method

All seven structures use identical settings:

| Setting            | Value                                        |
|--------------------|----------------------------------------------|
| Functional         | PBE                                          |
| Dispersion         | DFT-D3 (Grimme, BJ damping, v4)              |
| ecutwfc / ecutrho  | 60 Ry / 480 Ry                               |
| Pseudopotentials   | PAW for Mg/H, USPP for Ni/Co (PSlibrary 1.x) |
| Smearing           | cold, degauss = 0.01 Ry (metals)             |
| Occupations (H2)   | fixed (insulating molecule)                  |
| Spin               | nspin=2 for TM-containing cells, mag init    |
| conv_thr           | 1e-9                                         |
| forc_conv_thr      | 1e-4 Ry/Bohr                                 |
| Cell optimisation  | vc-relax for solids; relax-only for H2/box   |

K-points per cell are chosen for ~0.03 Å-1 density (e.g. 12×12×8 for the
2-atom HCP primitive, 6×6×4 for the 16-atom supercell, 4×4×6 for the
48-atom hydride supercell).

## Layout

```
mg-nico/
├── README.md          this file
├── build_inputs.py    generates inputs/*.in from scratch (no external deps)
├── run_all.sh         loops pw.x over inputs/, writes outputs/ and logs/
├── enthalpy.py        parses outputs/, computes the 3 ΔH values + verdict
├── inputs/            7 QE input files (.in)
├── outputs/           7 SCF/vc-relax outputs (.out) -- populated by run_all.sh
├── logs/              per-job stderr
└── outdir/            QE scratch (.wfc, .save, etc.) -- safe to delete
```

## Reproduce

```bash
# (1) regenerate inputs (idempotent)
python3 build_inputs.py

# (2) run all 7 DFT calculations (sources env.sh for the GPU pw.x)
bash run_all.sh

# (3) compute the three ΔH and print the catalyst verdict
python3 enthalpy.py
```

`run_all.sh` skips any job whose output already contains `JOB DONE`, so
re-running it is cheap. You can also select a subset:
`bash run_all.sh mgni mgh2ni` to only do the Ni pathway.

## How the comparison works

For pathway X (reactant `met_X`, product `hyd_X`, n H2 absorbed):

    ΔH_per_H2 = [ E(hyd_X) - E(met_X) - n * E(H2) ] / n     [Ry]

multiplied by 1312.75 kJ/mol/Ry. The verdict is:

    ΔΔH = ΔH(dopant) - ΔH(pure)
    if ΔΔH > 0   -> dopant DESTABILISES hydride -> easier release -> BETTER
    if ΔΔH < 0   -> dopant STABILISES   hydride -> harder release -> WORSE

Experimental ΔH references (for sanity):
* MgH2 desorption:    -75   kJ/mol H2 (Bogdanovic et al., 1999)
* Mg2NiH4 desorption: -64.5 kJ/mol H2 (Reilly & Wiswall, 1968)

## Caveats

- **No ZPE.** Static-lattice DFT only. ZPE typically shifts ΔH by
  +6 to +10 kJ/mol H2 *uniformly* for these systems, so the inter-pathway
  comparison is robust without it.
- **Single dopant configuration.** The two TM atoms are placed at the
  far ends of the supercell to minimise dopant-dopant interaction.
  A full study would average over symmetry-inequivalent placements.
- **No kinetics.** ΔH speaks only to thermodynamics. A catalyst can also
  lower activation barriers; that's not addressed here.
