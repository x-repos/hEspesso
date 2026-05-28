"""
Compute and compare hydrogenation enthalpies (kJ/mol H2) for three pathways:

    R1  pure         :  Mg          +    H2  ->  MgH2
    R2  Ni dopant 12%:  Mg14Ni2     + 16 H2  ->  Mg14Ni2H32
    R3  Co dopant 12%:  Mg14Co2     + 16 H2  ->  Mg14Co2H32

All seven SCF energies are read from outputs/ in this folder. The doped metals
and hydrides use the *same* 2x2x2 supercell topology, so the comparison
between Ni and Co is apples-to-apples (only the dopant species changes).

The catalyst is "better" (for hydrogen storage) if it makes ΔH per H2
less negative -- i.e. easier to release H2 (lower desorption temperature).
"""

from pathlib import Path
import re
import sys

HERE = Path(__file__).resolve().parent
OUT  = HERE / "outputs"

RY_TO_KJ = 1312.75


def total_energy(out_path):
    """Final converged total energy in Ry, or None."""
    if not out_path.exists():
        return None
    energies = []
    with out_path.open() as f:
        for line in f:
            m = re.match(r"!\s+total energy\s+=\s+(-?\d+\.\d+)\s+Ry", line)
            if m:
                energies.append(float(m.group(1)))
    return energies[-1] if energies else None


def main():
    species = ("mg", "h2", "mgh2", "mgni", "mgco", "mgh2ni", "mgh2co")
    E = {s: total_energy(OUT / f"{s}.out") for s in species}

    # report what's available
    print("=" * 74)
    print("Total energies (Ry) from outputs/")
    print("=" * 74)
    print(f"  {'species':10s}  {'cell content':24s}  {'E (Ry)':>16s}")
    print("  " + "-" * 56)
    info = {
        "mg":     ("2 Mg              (1 unit Mg HCP)",      "R1 reactant"),
        "h2":     ("2 H              (1 H2 molecule)",      "all reactions"),
        "mgh2":   ("2 Mg + 4 H       (2 fu rutile MgH2)",   "R1 product"),
        "mgni":   ("14 Mg + 2 Ni     (12.5% Ni @ Mg site)", "R2 reactant"),
        "mgco":   ("14 Mg + 2 Co     (12.5% Co @ Mg site)", "R3 reactant"),
        "mgh2ni": ("14 Mg + 2 Ni + 32 H  (Ni-doped MgH2)",  "R2 product"),
        "mgh2co": ("14 Mg + 2 Co + 32 H  (Co-doped MgH2)",  "R3 product"),
    }
    for s in species:
        e = E[s]
        e_str = f"{e:16.6f}" if e is not None else "       (missing)"
        print(f"  {s:10s}  {info[s][0]:24s}  {e_str}    [{info[s][1]}]")

    missing = [s for s in species if E[s] is None]
    if missing:
        print(f"\n(skipping reactions with missing outputs: {missing})")

    # --- ΔH per H2 ---------------------------------------------------------
    # Cell counts:
    #   mg.out          : 2 Mg cell (1 "Mg2" unit; 2 Mg per cell)
    #   mgh2.out        : 6 atoms (2 fu MgH2 = 2 Mg + 4 H = 2 H2 absorbed)
    #   mgni/mgco.out   : 16 metal sites (14 Mg + 2 TM)
    #   mgh2ni/mgh2co   : 14 Mg + 2 TM + 32 H = 16 H2 absorbed
    pathways = [
        ("R1  pure         Mg     +    H2  -> MgH2",
         "mg", "mgh2", 2),
        ("R2  Ni dopant 12% Mg14Ni2 + 16 H2 -> Mg14Ni2H32",
         "mgni", "mgh2ni", 16),
        ("R3  Co dopant 12% Mg14Co2 + 16 H2 -> Mg14Co2H32",
         "mgco", "mgh2co", 16),
    ]

    print()
    print("=" * 74)
    print("ΔH of hydrogenation per H2 (static lattice, PBE+D3)")
    print("=" * 74)
    print(f"  {'pathway':50s}  {'ΔH (kJ/mol H2)':>14s}")
    print("  " + "-" * 66)
    results = []
    for label, met_key, hyd_key, n_h2 in pathways:
        e_m, e_h, e_h2 = E[met_key], E[hyd_key], E["h2"]
        if None in (e_m, e_h, e_h2):
            print(f"  {label:50s}  (incomplete)")
            results.append((label, None))
            continue
        dE_ry = (e_h - e_m - n_h2 * e_h2) / n_h2
        dH_kj = dE_ry * RY_TO_KJ
        results.append((label, dH_kj))
        print(f"  {label:50s}  {dH_kj:+14.2f}")

    # --- Catalyst verdict --------------------------------------------------
    print()
    print("=" * 74)
    print("Catalyst effect (less-negative ΔH means easier H2 release = BETTER)")
    print("=" * 74)
    dH_r1 = results[0][1]
    if dH_r1 is None:
        print("  R1 unavailable; cannot compute catalyst shift.")
        return
    print(f"  R1 pure        ΔH = {dH_r1:+7.2f} kJ/mol H2   (baseline)")
    for label, dH in results[1:]:
        if dH is None:
            tm = "Ni" if "Ni" in label else "Co"
            print(f"  {tm} pathway not yet computed.")
            continue
        shift = dH - dH_r1
        tm = "Ni" if "Ni" in label else "Co"
        verdict = "BETTER (destabilises hydride)" if shift > 0 else "WORSE (stabilises hydride)"
        print(f"  {tm} dopant    ΔH = {dH:+7.2f} kJ/mol H2   "
              f"ΔΔH = {shift:+6.2f}   {verdict}")

    # --- Reference values --------------------------------------------------
    print()
    print("Reference (experimental, for sanity):")
    print("  MgH2  desorption     ΔH ≈ -75   kJ/mol H2  (Bogdanovic 1999)")
    print("  Mg2NiH4 desorption   ΔH ≈ -64.5 kJ/mol H2  (Reilly-Wiswall 1968)")


if __name__ == "__main__":
    main()
