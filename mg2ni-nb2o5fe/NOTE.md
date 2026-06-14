# NOTE — scope, method, and comparability of this study

_Last reviewed: 2026-06-13._

This folder is a **separate study** from `../mg-nico/`. Read this before reusing
any number from here.

## What this study is

- **System:** Nb₂O₅ + Fe additives on Mg₂NiH₄ (Nb and Fe co-doping, plus the
  in-situ Nb₂O₅ ball-milling reduction).
- **Not** the Ni/Co single-atom catalyst study — that is `../mg-nico/`.

## Method used here (and why it is not comparable to `../mg-nico/`)

| Setting | this study (`mg2ni-nb2o5fe`) | `../mg-nico/` |
|---|---|---|
| Hubbard U | **U = 9.0 eV on Ni-3d, 4.0 eV on Fe-3d** | none |
| Geometry | `calculation='scf'` (single point, frozen cell) | `vc-relax` (full cell+ion) |
| Dispersion | `vdw_corr='dft-d3'` (default **zero-damping**, D3 v3) | D3-**BJ** (`dftd3_version=4`) |
| Spin | `nspin=2` for the Fe-doped variants and elemental Ni; `nspin=1` for pristine and Nb-only cells | per-cell spin ground state |
| Pristine refs | imported from external `mgh2-cif/` | computed in-project |

**Do not mix energies across the two projects.** A DFT total energy is only
meaningful as a difference within one consistent set of settings. The same
compound lands on different numbers here vs `mg-nico` because of the +U,
dispersion damping, relaxation, and reference-state differences:

| ΔH per H₂ (kJ/mol) | here | `mg-nico` |
|---|---|---|
| MgH₂ | −67.95 | −63.75 |
| Mg₂NiH₄ (pristine) | −64.41 | −66.97 |

Note: the +U here only touches the **Ni-containing** cells (the `HUBBARD` block
declares `U Ni-3d` only); MgH₂ has no Ni, so its difference vs `mg-nico` comes
purely from the D3-damping / relaxation / reference differences, not from U.

The pristine Mg₂NiH₄ = −64.41 matching experiment (~−64) is **tuned/coincidental**
(this specific U with a frozen geometry, ZPE omitted), not evidence it is "more
correct." For the Ni/Co paper, use `../mg-nico/` exclusively.

## Reproducibility gap

The pristine reference energies in `enthalpy.py` are read from a sibling
`mgh2-cif/` directory (`mg2nih4-28.pwo`, `mg2ni.pwo`, `mg.pwo`, `h2.pwo`).
**That directory no longer exists on disk** (`/home/x/Workspace/espresso/mgh2-cif/`
is gone). Only the input `inputs/mg2nih4-28.pwi` survives. So the pristine
branch of the ΔH chain here is currently **not reproducible** and would need
re-running before the absolute values are trustworthy.

## "Mg₂NiH₄ cell looks smaller than Mg₂Ni" — not an error

The two cells hold a different number of formula units (Z):

| | structure | atoms | Z | volume | per f.u. |
|---|---|---|---|---|---|
| Mg₂Ni | hexagonal P6₂22 | 18 (Mg₁₂Ni₆) | 6 | ~296 Å³ | 49.3 Å³ |
| Mg₂NiH₄ | monoclinic C2/c (primitive) | 28 (Mg₈Ni₄H₁₆) | 4 | ~260 Å³ | 64.9 Å³ |

- Per **whole cell**: the hydride is ~12% smaller — only because it packs 4 f.u.
  against 6 (smaller box, fewer repeat units).
- Per **formula unit**: the hydride is **+32% larger** — exactly the documented
  Mg₂Ni → Mg₂NiH₄ hydrogenation expansion.

Never compare these two cells whole-to-whole. The figure showing a smaller
hydride box is a Z / crystal-system artifact, not a computational mistake.
(Volumes above are the `mg-nico` relaxed cells; the frozen-geometry cells used
here are ~273/309 Å³ but give the same per-f.u. conclusion.)
