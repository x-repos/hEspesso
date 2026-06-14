import re, sys, importlib.util
from pathlib import Path
HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("ri", HERE/"restart_input.py")
ri = importlib.util.module_from_spec(spec); spec.loader.exec_module(ri)

name = sys.argv[1]
src = HERE/"outputs"/f"{name}.out"
if not src.exists():                       # fall back to best archived attempt
    cands = sorted((HERE/"outputs"/"failed").glob(f"{name}.attempt*.out"))
    src = cands[-1] if cands else None
txt = src.read_text(errors="replace")
energies, geoms = ri.parse_trajectory(src)
cell, pos = geoms[-1]
base = (HERE/"inputs"/f"{name}.in").read_text()
# species order for starting_magnetization indices
species = re.findall(r"^\s*(\w+)\s+[\d.]+\s+\S+\.UPF", base.split("ATOMIC_SPECIES")[1], re.M)
mag = "".join(f"\n    starting_magnetization({i+1}) = {0.4 if sp=='Co' else 0.1}"
              for i,sp in enumerate(species) if sp in ("Co","Ni"))
base = base.replace("calculation      = 'vc-relax'","calculation      = 'relax'")
if "nspin" not in base:
    base = base.replace("degauss          = 0.01","degauss          = 0.01\n    nspin            = 2"+mag)
base = re.sub(r"(CELL_PARAMETERS angstrom\n)(?:.+\n){3}", lambda m:m.group(1)+"\n".join(cell)+"\n", base)
base = re.sub(r"(ATOMIC_POSITIONS angstrom\n)(?:.+\n)+?\n", lambda m:m.group(1)+"\n".join(pos)+"\n\n", base)
(HERE/"inputs"/f"{name}.in").write_text(base)
print(f"{name}: nspin=2 relax seeded from {src.name} (E_nspin1={energies[-1]:.6f}), magnetized {[s for s in species if s in ('Co','Ni')]}")
