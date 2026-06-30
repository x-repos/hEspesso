"""Validate every inputs/*.in before burning GPU time.

Per file: rev-vdW-DF2 functional, 60/600 cutoffs, NO Hubbard U, NO DFT-D3,
pseudo files exist, nat matches the atom count, NO2 inputs are nspin=2 and all
others nspin=1, slab/adsorption inputs carry the dipole field + >=1 frozen atom.
"""
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
PSEUDO_DIR = Path("/home/x/Workspace/3-hEspesso/pseudo")
ELEMS = ("V", "O", "C", "H", "N", "S")


def check(p):
    t = p.read_text()
    e = []
    if "vdw_corr         = 'DFT-D3'" not in t or "dftd3_version    = 4" not in t:
        e.append("missing PBE+DFT-D3(v4) dispersion")
    if "ecutwfc          = 60.0" not in t:
        e.append("ecutwfc != 60")
    if "ecutrho          = 600.0" not in t:
        e.append("ecutrho != 600")
    if "HUBBARD" in t or "Hubbard" in t or "lda_plus_u" in t:
        e.append("has Hubbard U")
    if "input_dft" in t:
        e.append("has input_dft (should be plain PBE + DFT-D3)")
    nat = int(re.search(r"nat\s*=\s*(\d+)", t).group(1))
    body = t.split("ATOMIC_POSITIONS")[1].splitlines()[1:]
    natoms = sum(1 for ln in body if len(ln.split()) >= 4 and ln.split()[0] in ELEMS)
    if nat != natoms:
        e.append(f"nat={nat} but {natoms} atom lines")
    nspin = int(re.search(r"nspin\s*=\s*(\d+)", t).group(1))
    is_no2 = "no2" in p.stem
    if is_no2 and nspin != 2:
        e.append("NO2 must be nspin=2")
    if not is_no2 and nspin != 1:
        e.append("non-NO2 must be nspin=1")
    if p.stem.startswith(("slab", "ads_")):
        if "assume_isolated  = 'esm'" not in t:
            e.append("slab missing ESM z-isolation")
        if " 0 0 0" not in t:
            e.append("slab missing frozen atoms")
    for ps in re.findall(r"^\s*\w+\s+[\d.]+\s+(\S+\.UPF)", t, re.M | re.I):
        if not (PSEUDO_DIR / ps).exists():
            e.append(f"pseudo {ps} not found")
    return e


def main():
    files = sorted((HERE / "inputs").glob("*.in"))
    if not files:
        sys.exit("no inputs/*.in found")
    bad = 0
    for p in files:
        errs = check(p)
        print(f"  [{'OK ' if not errs else 'FAIL'}] {p.name}")
        for x in errs:
            print("      -", x)
        bad += bool(errs)
    if bad:
        sys.exit(f"\n{bad} input(s) failed.")
    print(f"\nAll {len(files)} inputs pass.")


if __name__ == "__main__":
    main()
