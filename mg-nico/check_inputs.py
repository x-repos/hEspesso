"""
Sanity-check the 13 generated QE inputs before burning GPU time.

Verifies for every inputs/*.in:
  - composition matches the intended stoichiometry exactly
  - nat / ntyp are consistent with the atom list
  - every pseudopotential file exists in pseudo_dir
  - no Nb remnants (structures were recovered from the Nb-doped study)
  - nosym set on doped cells only (pristine cells keep symmetry)
  - cell matrix is non-singular and K_POINTS block present
"""

from pathlib import Path
from collections import Counter
import re
import sys

HERE = Path(__file__).resolve().parent

# name -> (composition, calc, nosym)
EXPECT = {
    "mg":          ({"Mg": 2},                          "vc-relax", False),
    "h2":          ({"H": 2},                           "relax",    False),
    "mgh2":        ({"Mg": 2, "H": 4},                  "vc-relax", False),
    "mgni":        ({"Mg": 15, "Ni": 1},                "vc-relax", True),
    "mgco":        ({"Mg": 15, "Co": 1},                "vc-relax", True),
    "mgh2ni":      ({"Mg": 15, "Ni": 1, "H": 32},       "vc-relax", True),
    "mgh2co":      ({"Mg": 15, "Co": 1, "H": 32},       "vc-relax", True),
    "mg2ni":       ({"Mg": 12, "Ni": 6},                "vc-relax", False),
    "mg2ni_ni":    ({"Mg": 11, "Ni": 7},                "vc-relax", True),
    "mg2ni_co":    ({"Mg": 11, "Co": 1, "Ni": 6},       "vc-relax", True),
    "mg2nih4":     ({"Mg": 8, "Ni": 4, "H": 16},        "vc-relax", False),
    "mg2nih4_ni":  ({"Mg": 7, "Ni": 5, "H": 16},        "vc-relax", True),
    "mg2nih4_co":  ({"Mg": 7, "Co": 1, "Ni": 4, "H": 16}, "vc-relax", True),
    "mg2ni_co_nisite":   ({"Mg": 12, "Co": 1, "Ni": 5}, "vc-relax", True),
    "mg2nih4_co_nisite": ({"Mg": 8, "Co": 1, "Ni": 3, "H": 16}, "vc-relax", True),
}


def det3(m):
    return (m[0][0]*(m[1][1]*m[2][2]-m[1][2]*m[2][1])
          - m[0][1]*(m[1][0]*m[2][2]-m[1][2]*m[2][0])
          + m[0][2]*(m[1][0]*m[2][1]-m[1][1]*m[2][0]))


def check(name, expect_comp, expect_calc, expect_nosym):
    path = HERE / "inputs" / f"{name}.in"
    errors = []
    if not path.exists():
        return [f"{name}: missing input file"]
    text = path.read_text()
    lines = text.splitlines()

    if "Nb" in text:
        errors.append("contains 'Nb'")

    calc = re.search(r"calculation\s*=\s*'(\S+)'", text).group(1)
    if calc != expect_calc:
        errors.append(f"calculation={calc}, expected {expect_calc}")

    nat = int(re.search(r"nat\s*=\s*(\d+)", text).group(1))
    ntyp = int(re.search(r"ntyp\s*=\s*(\d+)", text).group(1))
    has_nosym = "nosym" in text
    if has_nosym != expect_nosym:
        errors.append(f"nosym={'set' if has_nosym else 'absent'}, "
                      f"expected {'set' if expect_nosym else 'absent'}")
    if "HUBBARD" in text or re.search(r"\bU \w+-\dd", text):
        errors.append("contains a Hubbard U block")

    # atoms
    j = lines.index("ATOMIC_POSITIONS angstrom")
    atoms = []
    for line in lines[j+1:]:
        p = line.split()
        if len(p) != 4:
            break
        atoms.append(p[0])
    comp = dict(Counter(atoms))
    if comp != expect_comp:
        errors.append(f"composition {comp}, expected {expect_comp}")
    if nat != len(atoms):
        errors.append(f"nat={nat} but {len(atoms)} atom lines")
    if ntyp != len(set(atoms)):
        errors.append(f"ntyp={ntyp} but {len(set(atoms))} species present")

    # species block must list exactly the species present, with existing pseudos
    pseudo_dir = Path(re.search(r"pseudo_dir\s*=\s*'([^']+)'", text).group(1))
    i = lines.index("ATOMIC_SPECIES")
    declared = []
    for line in lines[i+1:]:
        p = line.split()
        if len(p) != 3:
            break
        declared.append(p[0])
        if not (pseudo_dir / p[2]).exists():
            errors.append(f"pseudo {p[2]} not found in {pseudo_dir}")
    if set(declared) != set(atoms) or len(declared) != ntyp:
        errors.append(f"ATOMIC_SPECIES {declared} inconsistent with atoms")

    # cell + k-points
    k = lines.index("CELL_PARAMETERS angstrom")
    cell = [[float(x) for x in lines[k+r+1].split()] for r in range(3)]
    if abs(det3(cell)) < 1.0:
        errors.append(f"degenerate cell, |det|={abs(det3(cell)):.3f}")
    if "K_POINTS automatic" not in text:
        errors.append("missing K_POINTS block")

    return [f"{name}: {e}" for e in errors]


def main():
    all_errors = []
    for name, (comp, calc, nosym) in EXPECT.items():
        errs = check(name, comp, calc, nosym)
        status = "OK " if not errs else "FAIL"
        print(f"  [{status}] {name:12s} {'+'.join(f'{v}{k}' for k, v in comp.items())}")
        all_errors.extend(errs)
    if all_errors:
        print("\nERRORS:")
        for e in all_errors:
            print(f"  - {e}")
        sys.exit(1)
    print(f"\nAll {len(EXPECT)} inputs pass.")


if __name__ == "__main__":
    main()
