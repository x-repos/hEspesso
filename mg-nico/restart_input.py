"""
Re-seed inputs/<name>.in from a geometry visited by a (possibly failed)
vc-relax in outputs/<name>.out, so a fresh BFGS can continue the
minimisation.

By default picks the LOWEST-energy geometry of the trajectory, not the last:
a failed BFGS often ends on an uphill oscillation branch. Optionally tightens
the BFGS trust radius (optimizer knob only -- the physics is untouched).

The stale output (and stderr log) are archived to outputs/failed/ so
run_all.sh reruns the job.

Usage:
    python3 restart_input.py <name> [--trust 0.2]
"""

import argparse
import re
import shutil
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent


def parse_trajectory(out_path):
    """Return (energies, geoms) where geoms[k] = (cell_lines, pos_lines) is
    the geometry written AFTER the SCF that produced energies[k]."""
    lines = out_path.read_text().splitlines()
    energies, geoms = [], []
    cell, pos = None, None
    i = 0
    while i < len(lines):
        line = lines[i]
        m = re.match(r"!\s+total energy\s+=\s+(-?\d+\.\d+)\s+Ry", line)
        if m:
            energies.append(float(m.group(1)))
            cell, pos = None, None
        elif line.startswith("CELL_PARAMETERS"):
            assert "angstrom" in line, f"unexpected units: {line}"
            cell = lines[i+1:i+4]
            i += 3
        elif line.startswith("ATOMIC_POSITIONS"):
            assert "angstrom" in line, f"unexpected units: {line}"
            pos = []
            j = i + 1
            while j < len(lines):
                p = lines[j].split()
                if len(p) not in (4, 7):     # optional 0/1 relax flags
                    break
                try:
                    [float(v) for v in p[1:4]]
                except ValueError:
                    break
                pos.append(lines[j])
                j += 1
            geoms.append((cell, pos))      # geometry after current energy
            i = j - 1
        i += 1
    return energies, geoms


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("name")
    ap.add_argument("--trust", type=float, default=None,
                    help="set trust_radius_max in &IONS (bohr)")
    args = ap.parse_args()

    out_path = HERE / "outputs" / f"{args.name}.out"
    in_path = HERE / "inputs" / f"{args.name}.in"
    energies, geoms = parse_trajectory(out_path)
    if not energies:
        sys.exit(f"no energies found in {out_path}")

    best = min(range(len(energies)), key=lambda k: energies[k])
    print(f"{args.name}: {len(energies)} SCF energies, "
          f"best #{best+1} = {energies[best]:.8f} Ry "
          f"(last = {energies[-1]:.8f} Ry)")

    # geometry evaluated at energies[best] is the one written after
    # energies[best-1]; for best == 0 the input geometry already was best
    if best == 0:
        sys.exit("initial geometry already lowest -- restart would not move; "
                 "inspect the run instead")
    cell, pos = geoms[best - 1]
    if cell is None:
        # ions-only step (no cell change printed); reuse the most recent cell
        for k in range(best - 2, -1, -1):
            if geoms[k][0] is not None:
                cell = geoms[k][0]
                break
    assert cell and pos, "could not recover a full geometry"

    text = in_path.read_text()
    text = re.sub(r"(CELL_PARAMETERS angstrom\n)(?:.+\n){3}",
                  lambda m: m.group(1) + "\n".join(cell) + "\n", text)
    pos_block = "\n".join(pos) + "\n"
    text = re.sub(r"(ATOMIC_POSITIONS angstrom\n)(?:.+\n)+?\n",
                  lambda m: m.group(1) + pos_block + "\n", text)
    if args.trust is not None:
        if "trust_radius_max" in text:
            text = re.sub(r"trust_radius_max\s*=\s*\S+",
                          f"trust_radius_max = {args.trust}", text)
        else:
            text = text.replace("ion_dynamics = 'bfgs'",
                                f"ion_dynamics = 'bfgs'\n    trust_radius_max = {args.trust}")
    in_path.write_text(text)
    print(f"re-seeded {in_path} from geometry #{best}"
          + (f", trust_radius_max={args.trust}" if args.trust else ""))

    failed_dir = HERE / "outputs" / "failed"
    failed_dir.mkdir(exist_ok=True)
    n = len(list(failed_dir.glob(f"{args.name}.attempt*.out"))) + 1
    shutil.move(out_path, failed_dir / f"{args.name}.attempt{n}.out")
    err = HERE / "logs" / f"{args.name}.err"
    if err.exists():
        shutil.move(err, failed_dir / f"{args.name}.attempt{n}.err")
    print(f"archived stale output to outputs/failed/{args.name}.attempt{n}.out")


if __name__ == "__main__":
    main()
