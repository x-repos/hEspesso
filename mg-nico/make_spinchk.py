"""
Generate inputs/<name>_spinchk.in: the converged final geometry of
outputs/<name>.out rerun as vc-relax with nspin=2 and a starting moment on
the TM species.

Purpose: verify that the nspin=1 ground state used in the pathways is
spin-stable. If the moment converges to ~0, the nspin=1 energy stands; if a
finite moment survives with lower energy, the spin-polarized energy must
replace it (see the Co bistability incident: nonmagnetic Co d-states at
E_F of a metal can even be electronically bistable).

Usage: python3 make_spinchk.py <name> [mag]
"""

import re
import sys
from pathlib import Path
import importlib.util

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("ri", HERE / "restart_input.py")
ri = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ri)


def main():
    name = sys.argv[1]
    mag = float(sys.argv[2]) if len(sys.argv) > 2 else 0.4
    out_path = HERE / "outputs" / f"{name}.out"
    text_out = out_path.read_text(errors="replace")
    if "bfgs converged" not in text_out:
        sys.exit(f"{name} is not converged -- spinchk needs a converged geometry")
    energies, geoms = ri.parse_trajectory(out_path)
    cell, pos = geoms[-1]
    assert cell and pos, "missing final geometry"

    text = (HERE / "inputs" / f"{name}.in").read_text()
    # magnetize every TM species present (index in ATOMIC_SPECIES order)
    species = re.findall(r"^\s*(\w+)\s+[\d.]+\s+\S+\.UPF\s*$",
                         text.split("ATOMIC_SPECIES")[1], re.M)
    mag_lines = "".join(f"\n    starting_magnetization({i+1}) = {mag}"
                        for i, sp in enumerate(species) if sp in ("Ni", "Co"))
    if not mag_lines:
        sys.exit(f"{name} has no TM species")
    text = text.replace(f"prefix           = '{name}'",
                        f"prefix           = '{name}_spinchk'")
    # Single-point SCF at the converged geometry: we only need to know whether
    # a finite moment survives (and lowers the energy), not to re-relax.
    text = text.replace("calculation      = 'vc-relax'",
                        "calculation      = 'scf'")
    text = text.replace("degauss          = 0.01",
                        "degauss          = 0.01\n    nspin            = 2" + mag_lines)
    text = re.sub(r"(CELL_PARAMETERS angstrom\n)(?:.+\n){3}",
                  lambda m: m.group(1) + "\n".join(cell) + "\n", text)
    text = re.sub(r"(ATOMIC_POSITIONS angstrom\n)(?:.+\n)+?\n",
                  lambda m: m.group(1) + "\n".join(pos) + "\n\n", text)
    (HERE / "inputs" / f"{name}_spinchk.in").write_text(text)
    print(f"wrote inputs/{name}_spinchk.in  (E_nspin1 = {energies[-1]:.8f} Ry, "
          f"magnetized: {[sp for sp in species if sp in ('Ni', 'Co')]})")


if __name__ == "__main__":
    main()
