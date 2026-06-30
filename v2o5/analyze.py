"""Rank the six gases on V2O5 by adsorption energy (and, after postproc.sh, by
Bader charge transfer and work-function change), and compare to the experimental
sensitivity order acetone > nh3 > c2h4 > h2s > co > no2.

E_ads = E(slab+gas) - E(slab) - E(gas)   [eV; more negative = stronger binding]
For each gas the lowest-E_ads configuration among ads_<key>_s* is reported.
"""
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUT = HERE / "outputs"
RY_EV = 13.605693
TARGET = ["acetone", "nh3", "c2h4", "h2s", "co", "no2"]
GAS_NAT = {"acetone": 10, "nh3": 4, "c2h4": 6, "h2s": 3, "co": 2, "no2": 3}


def total_energy_ry(path):
    """Final '!' total energy (Ry), or None if the run is not a converged min."""
    path = Path(path)
    if not path.exists():
        return None
    t = path.read_text(errors="replace")
    if "JOB DONE" not in t:
        return None
    if "Maximum CPU time exceeded" in t:
        return None
    if "bfgs converged" not in t:
        return None
    m = re.findall(r"!\s+total energy\s+=\s+(-?\d+\.\d+)\s+Ry", t)
    return float(m[-1]) if m else None


def eads_ev(e_complex_ry, e_slab_ry, e_gas_ry):
    return (e_complex_ry - e_slab_ry - e_gas_ry) * RY_EV


def best_site(key, e_slab):
    """(E_ads_eV, prefix) of the lowest-E_ads converged config for `key`."""
    eg = total_energy_ry(OUT / f"gas_{key}.out")
    rows = []
    for p in sorted(OUT.glob(f"ads_{key}_s*.out")):
        ec = total_energy_ry(p)
        if None in (ec, eg, e_slab):
            continue
        rows.append((eads_ev(ec, e_slab, eg), p.stem))
    return min(rows) if rows else (None, None)


def _zval_map(out_path):
    """ZVAL per species from the QE 'atomic species valence mass' header table."""
    t = Path(out_path).read_text(errors="replace")
    return {sp: float(v) for sp, v in
            re.findall(r"^\s+([A-Z][a-z]?)\s+([\d.]+)\s+[\d.]+\s+\w", t, re.M)}


def _input_species(prefix):
    """Atom species in ATOMIC_POSITIONS order from inputs/<prefix>.in."""
    body = (HERE / "inputs" / f"{prefix}.in").read_text().split("ATOMIC_POSITIONS")[1]
    return [ln.split()[0] for ln in body.splitlines()[1:]
            if len(ln.split()) >= 4 and ln.split()[0] in ("V", "O", "C", "H", "N", "S")]


def work_function(prefix):
    """phi = V_vacuum - E_Fermi (eV). Reads pp/avg_<prefix>.dat (col 3 = planar
    average potential, Ry, vs z) + the run's Fermi energy. None if not yet made."""
    avg = HERE / "pp" / f"avg_{prefix}.dat"
    out = OUT / f"{prefix}.out"
    if not (avg.exists() and out.exists()):
        return None
    pot = [float(l.split()[2]) for l in avg.read_text().splitlines()
           if len(l.split()) >= 3 and l.split()[0].lstrip("-").replace(".", "").isdigit()]
    ef = re.findall(r"the Fermi energy is\s+(-?\d+\.\d+)\s+ev", out.read_text())
    if not pot or not ef:
        return None
    return max(pot) * RY_EV - float(ef[-1])


def bader_transfer(prefix, key):
    """Net charge donated by the adsorbate (e); + = molecule -> surface. Sum over
    the molecule's atoms (last GAS_NAT[key]) of (ZVAL - Bader population)."""
    acf = HERE / "pp" / f"{prefix}_ACF.dat"
    out = OUT / f"{prefix}.out"
    if not (acf.exists() and out.exists()):
        return None
    pops = [float(l.split()[4]) for l in acf.read_text().splitlines()
            if l.split() and l.split()[0].isdigit()]
    zval = _zval_map(out)
    species = _input_species(prefix)
    n = GAS_NAT[key]
    if len(pops) < n or any(s not in zval for s in species[-n:]):
        return None
    return sum(zval[s] - q for s, q in zip(species[-n:], pops[-n:]))


def descriptor_table():
    """E_ads, Bader dq, dphi per gas's winning complex; rank by each descriptor."""
    e_slab = total_energy_ry(OUT / "slab.out")
    phi_slab = work_function("slab")
    cols = {"eads": {}, "dq": {}, "dphi": {}}
    print(f"\n{'gas':10s} {'E_ads(eV)':>10s} {'dq(e)':>8s} {'dphi(eV)':>9s}  prefix")
    print("-" * 58)
    for key in TARGET:
        ea, pre = best_site(key, e_slab)
        if pre is None:
            print(f"{key:10s} {'(no data)':>10s}")
            continue
        dq = bader_transfer(pre, key)
        phi = work_function(pre)
        dphi = (phi - phi_slab) if None not in (phi, phi_slab) else None
        cols["eads"][key], cols["dq"][key], cols["dphi"][key] = ea, dq, dphi
        g = lambda x, w, d: (f"{x:{w}.{d}f}" if x is not None else f"{'--':>{w}}")
        print(f"{key:10s} {g(ea,10,3)} {g(dq,8,3)} {g(dphi,9,3)}  {pre}")
    fns = {"E_ads": lambda k: cols["eads"].get(k),
           "|dq| ": lambda k: (-abs(cols["dq"][k]) if cols["dq"].get(k) is not None else None),
           "dphi ": lambda k: cols["dphi"].get(k)}
    print()
    for label, fn in fns.items():
        avail = [k for k in TARGET if fn(k) is not None]
        rank = sorted(avail, key=fn)
        tag = "MATCH" if rank == TARGET else "MISMATCH"
        print(f"  rank by {label}: {' > '.join(rank) or '(none)'}  [{tag}]")


def main():
    e_slab = total_energy_ry(OUT / "slab.out")
    print(f"{'gas':10s} {'E_ads (eV)':>12s}  {'best config':>16s}")
    print("-" * 42)
    res = {}
    for key in TARGET:
        ea, site = best_site(key, e_slab)
        res[key] = ea
        print(f"{key:10s} {ea:12.3f}  {site:>16s}" if ea is not None
              else f"{key:10s} {'(pending)':>12s}")
    ranked = [k for k, v in sorted(res.items(),
              key=lambda kv: (kv[1] if kv[1] is not None else 1e9)) if v is not None]
    print("\nDFT E_ads ranking (strongest first):", " > ".join(ranked) or "(pending)")
    print("Experimental order                 :", " > ".join(TARGET))
    if ranked == TARGET:
        print("=> E_ads MATCHES experiment")
    elif ranked:
        print("=> E_ads differs; consult descriptor table + DESIGN.md escalation")
    if (HERE / "pp").exists() and any((HERE / "pp").glob("*_ACF.dat")):
        descriptor_table()


if __name__ == "__main__":
    main()
