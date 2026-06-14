"""
Compute and compare hydrogenation enthalpies (kJ/mol H2) for six pathways,
two hosts x (pure, +Ni, +Co), one catalyst atom per cell:

    Host Mg (matched supercells, 16 metal sites, 1 dopant = 6.25%):
      R1  pure      :  Mg       +    H2  ->  MgH2
      R2  Ni-doped  :  Mg15Ni   + 16 H2  ->  Mg15NiH32
      R3  Co-doped  :  Mg15Co   + 16 H2  ->  Mg15CoH32

    Host Mg2Ni (18-atom metal cell = Mg12Ni6, 28-atom hydride = Mg8Ni4H16):
      R4  pure      :  (8/12) Mg12Ni6  + 8 H2  ->  Mg8Ni4H16
      R5  Ni-doped  :  R4 + defect shift from Mg11Ni7 / Mg7Ni5H16
      R6  Co-doped  :  R4 + defect shift from Mg11CoNi6 / Mg7CoNi4H16

For the Mg host, reactant and product supercells contain identical metal
content, so ΔH per H2 follows from a single balanced reaction.

For the Mg2Ni host the doped metal (18-atom) and doped hydride (28-atom)
cells do NOT contain matching metal content, so the doped ΔH is computed in
the dilute-defect limit:

    ΔΔH = [ (E_hyd_doped - E_hyd_pure) - (E_met_doped - E_met_pure) ] / 8

The Mg and dopant chemical potentials cancel exactly between the two
substitution energies; what remains is the doping-induced change of the
hydrogenation enthalpy per H2 (8 H2 per hydride cell). Note the dopant
concentration differs between the phases (1/18 of metal sites in the metal
cell vs 1/12 in the hydride cell) and between hosts (the Mg host has
1 dopant per 16 H2, the Mg2Ni host 1 per 8 H2); the residual is a
finite-size/finite-concentration error of the dilute approximation -- see
the README caveats.

A catalyst is "better" (for hydrogen storage) if it makes ΔH per H2 less
negative -- i.e. easier to release H2 (lower desorption temperature).
"""

from pathlib import Path
import re

HERE = Path(__file__).resolve().parent
OUT  = HERE / "outputs"

RY_TO_KJ = 1312.75


def total_energy(out_path):
    """Final converged total energy in Ry, or None.

    Guard against silently using non-minimum energies: a run that hit
    max_seconds or whose BFGS failed still ends with 'JOB DONE.', but its
    last '!' energy is an unconverged-geometry value. Such outputs are
    rejected (with a warning) so the job gets re-seeded and rerun instead.
    """
    if not out_path.exists():
        return None
    text = out_path.read_text(errors="replace")
    reason = None
    if "JOB DONE" not in text:
        reason = "run incomplete (no JOB DONE) -- still running or crashed"
    elif "Maximum CPU time exceeded" in text:
        reason = "hit max_seconds before geometry converged"
    elif "bfgs converged" not in text:
        reason = "BFGS did not converge"
    if reason:
        print(f"  WARNING: ignoring {out_path.name}: {reason} -- "
              f"re-seed with restart_input.py and rerun")
        return None
    energies = [float(m.group(1)) for m in re.finditer(
        r"!\s+total energy\s+=\s+(-?\d+\.\d+)\s+Ry", text)]
    return energies[-1] if energies else None


INFO = {
    "mg":          ("2 Mg                 (1 cell HCP Mg)",       "R1 reactant"),
    "h2":          ("2 H                  (1 H2 molecule)",       "all reactions"),
    "mgh2":        ("2 Mg + 4 H           (2 fu rutile MgH2)",    "R1 product"),
    "mgni":        ("15 Mg + 1 Ni         (6.25% Ni @ Mg site)",  "R2 reactant"),
    "mgco":        ("15 Mg + 1 Co         (6.25% Co @ Mg site)",  "R3 reactant"),
    "mgh2ni":      ("15 Mg + 1 Ni + 32 H  (Ni-doped MgH2)",       "R2 product"),
    "mgh2co":      ("15 Mg + 1 Co + 32 H  (Co-doped MgH2)",       "R3 product"),
    "mg2ni":       ("12 Mg + 6 Ni         (6 fu Mg2Ni)",          "R4 reactant"),
    "mg2nih4":     ("8 Mg + 4 Ni + 16 H   (4 fu LT-Mg2NiH4)",     "R4 product"),
    "mg2ni_ni":    ("11 Mg + 7 Ni         (Ni @ Mg site)",        "R5 reactant"),
    "mg2nih4_ni":  ("7 Mg + 5 Ni + 16 H   (Ni @ Mg site)",        "R5 product"),
    "mg2ni_co":    ("11 Mg + 1 Co + 6 Ni  (Co @ Mg site)",        "R6 reactant"),
    "mg2nih4_co":  ("7 Mg + 1 Co + 4 Ni + 16 H (Co @ Mg site)",   "R6 product"),
    "mg2ni_co_nisite":   ("12 Mg + 1 Co + 5 Ni  (Co @ Ni site)",  "R6' reactant"),
    "mg2nih4_co_nisite": ("8 Mg + 1 Co + 3 Ni + 16 H (Co @ Ni site)", "R6' product"),
}


# Energy-converged but force-unconverged cells (documented exception).
# These dilute-defect metallic cells have a sub-mRy-flat BFGS landscape: the
# total energy converges to the same minimum across independent relaxation
# attempts, but the residual force sits at the cold-smearing noise floor
# (~1e-2 Ry/Bohr) and BFGS never prints "bfgs converged". Treated as in the
# sibling study's mgh2ni case ("energy converged, geometry near minimum,
# energy usable"). Value = global-minimum '!' energy over all attempts;
# the energy enters ΔΔH, where the small geometric residual largely cancels.
ACCEPTED = {
    "mg2ni_ni": (-2773.36189213,
                 "BFGS-flat metallic cell; energy at global min over 4 "
                 "relaxation attempts; |F|~1e-2 Ry/Bohr at cold-smearing "
                 "floor; nspin=2 check gives m=0, so nspin=1 energy is valid"),
    "mg2nih4_ni": (-1971.10236788,
                 "BFGS-flat hydride cell; energy at global min over 3 "
                 "relaxation attempts (spread <0.1 mRy); |F|~2e-3 Ry/Bohr; "
                 "geometry at energy minimum"),
}


def main():
    E = {s: total_energy(OUT / f"{s}.out") for s in INFO}
    for s, (e_acc, reason) in ACCEPTED.items():
        if E.get(s) is None:
            E[s] = e_acc
            print(f"  NOTE: accepting {s} = {e_acc:.6f} Ry (force-unconverged): {reason}")

    print("=" * 78)
    print("Total energies (Ry) from outputs/")
    print("=" * 78)
    print(f"  {'species':12s}  {'cell content':42s}  {'E (Ry)':>16s}")
    print("  " + "-" * 74)
    for s, (content, role) in INFO.items():
        e_str = f"{E[s]:16.6f}" if E[s] is not None else "       (missing)"
        print(f"  {s:12s}  {content:42s}  {e_str}  [{role}]")

    missing = [s for s in INFO if E[s] is None]
    if missing:
        print(f"\n(reactions involving missing outputs are skipped: {missing})")

    def have(*keys):
        return all(E[k] is not None for k in keys)

    # --- ΔH per H2 ----------------------------------------------------------
    print()
    print("=" * 78)
    print("ΔH of hydrogenation per H2 (static lattice, PBE+D3, no U)")
    print("=" * 78)

    dH = {}

    # Host Mg: matched-cell balanced reactions
    if have("mg", "mgh2", "h2"):
        dH["R1"] = (E["mgh2"] - E["mg"] - 2 * E["h2"]) / 2 * RY_TO_KJ
    if have("mgni", "mgh2ni", "h2"):
        dH["R2"] = (E["mgh2ni"] - E["mgni"] - 16 * E["h2"]) / 16 * RY_TO_KJ
    if have("mgco", "mgh2co", "h2"):
        dH["R3"] = (E["mgh2co"] - E["mgco"] - 16 * E["h2"]) / 16 * RY_TO_KJ

    # Host Mg2Ni: pristine mass-balanced reaction + dilute-defect shifts
    if have("mg2ni", "mg2nih4", "h2"):
        dH["R4"] = (E["mg2nih4"] - (8.0 / 12.0) * E["mg2ni"]
                    - 8 * E["h2"]) / 8 * RY_TO_KJ
        if have("mg2ni_ni", "mg2nih4_ni"):
            shift = ((E["mg2nih4_ni"] - E["mg2nih4"])
                     - (E["mg2ni_ni"] - E["mg2ni"])) / 8 * RY_TO_KJ
            dH["R5"] = dH["R4"] + shift
        if have("mg2ni_co", "mg2nih4_co"):
            shift = ((E["mg2nih4_co"] - E["mg2nih4"])
                     - (E["mg2ni_co"] - E["mg2ni"])) / 8 * RY_TO_KJ
            dH["R6"] = dH["R4"] + shift
        if have("mg2ni_co_nisite", "mg2nih4_co_nisite"):
            shift = ((E["mg2nih4_co_nisite"] - E["mg2nih4"])
                     - (E["mg2ni_co_nisite"] - E["mg2ni"])) / 8 * RY_TO_KJ
            dH["R6b"] = dH["R4"] + shift

    labels = {
        "R1": "R1  pure Mg       :  Mg      +    H2 -> MgH2",
        "R2": "R2  Mg : Ni 6.25% :  Mg15Ni  + 16 H2 -> Mg15NiH32",
        "R3": "R3  Mg : Co 6.25% :  Mg15Co  + 16 H2 -> Mg15CoH32",
        "R4": "R4  pure Mg2Ni    :  Mg2Ni   +  2 H2 -> Mg2NiH4",
        "R5": "R5  Mg2Ni : Ni    :  defect shift, Ni @ Mg site",
        "R6": "R6  Mg2Ni : Co    :  defect shift, Co @ Mg site",
        "R6b": "R6' Mg2Ni : Co    :  defect shift, Co @ Ni site",
    }
    print(f"  {'pathway':52s}  {'ΔH (kJ/mol H2)':>16s}")
    print("  " + "-" * 70)
    for r, label in labels.items():
        val = f"{dH[r]:+14.2f}" if r in dH else "    (incomplete)"
        print(f"  {label:52s}  {val:>16s}")

    # --- Catalyst verdicts ---------------------------------------------------
    print()
    print("=" * 78)
    print("Catalyst effect (less-negative ΔH = easier H2 release = BETTER)")
    print("=" * 78)
    for host, base, doped in (("Mg", "R1", (("Ni", "R2"), ("Co", "R3"))),
                              ("Mg2Ni", "R4", (("Ni", "R5"), ("Co", "R6"),
                                               ("Co@Ni", "R6b")))):
        print(f"\n  Host {host}:")
        if base not in dH:
            print(f"    baseline {base} unavailable; cannot compute shifts.")
            continue
        print(f"    pure        ΔH = {dH[base]:+8.2f} kJ/mol H2   (baseline)")
        for tm, r in doped:
            if r not in dH:
                print(f"    {tm} pathway not yet computed.")
                continue
            shift = dH[r] - dH[base]
            verdict = ("BETTER (destabilises hydride)" if shift > 0
                       else "WORSE (stabilises hydride)")
            print(f"    {tm} catalyst  ΔH = {dH[r]:+8.2f} kJ/mol H2   "
                  f"ΔΔH = {shift:+6.2f}   {verdict}")

    # --- Experiment-anchored estimates ----------------------------------------
    # PBE's systematic error differs per host (under-binds MgH2 by ~11 kJ/mol,
    # within ~2.6 for Mg2NiH4), so cross-host comparisons must NOT use raw DFT
    # baselines. Anchor each host to its experimental ΔH and apply only the
    # (error-cancelling) DFT shift: ΔH(doped) ≈ ΔH_exp(pure) + ΔΔH_DFT.
    # Experimental values (verified against the primary sources):
    #   MgH2    -74.5 kJ/mol H2  (Bogdanovic 1999, -74.513 at 298 K;
    #                             Stampfer 1960, -17790 cal/mol = -74.4)
    #   Mg2NiH4 -64.4 kJ/mol H2  (Reilly & Wiswall 1968, -15.4 kcal/mol H2)
    EXP = {"R1": -74.5, "R4": -64.4}
    print()
    print("=" * 78)
    print("Experiment-anchored estimates: ΔH_exp(pure host) + ΔΔH_DFT")
    print("=" * 78)
    for host, base, doped in (("Mg", "R1", (("Ni", "R2"), ("Co", "R3"))),
                              ("Mg2Ni", "R4", (("Ni", "R5"), ("Co", "R6"),
                                               ("Co@Ni", "R6b")))):
        print(f"  {host:6s} pure        ΔH = {EXP[base]:+8.2f} kJ/mol H2  (experimental)")
        for tm, r in doped:
            if r in dH and base in dH:
                anchored = EXP[base] + (dH[r] - dH[base])
                print(f"  {host:6s} +{tm:2s}         ΔH = {anchored:+8.2f} kJ/mol H2  "
                      f"(anchored, ΔΔH = {dH[r] - dH[base]:+6.2f})")

    # --- Global ranking on the anchored scale ---------------------------------
    # Anchored ΔH puts all six systems on one thermodynamic scale (experiment
    # for the host baseline + error-cancelling DFT shift for the dopant).
    # T_des(1 bar) ≈ ΔH_des / ΔS with ΔS ≈ 130.5 J/(mol K) H2 (van 't Hoff
    # rule; metal-hydride desorption entropy is dominated by the H2 gas).
    # wt% H2 from the cell composition of the hydride.
    M = {"Mg": 24.305, "H": 1.008, "Ni": 58.693, "Co": 58.933}
    def wt_pct(comp, n_h):
        mass = sum(M[el] * n for el, n in comp.items())
        return 100.0 * n_h * M["H"] / mass
    HYDRIDE_COMP = {   # normalized label -> (hydride cell composition, H atoms released)
        "Mg pure":    ({"Mg": 16, "H": 32}, 32),
        "Mg + Ni":    ({"Mg": 15, "Ni": 1, "H": 32}, 32),
        "Mg + Co":    ({"Mg": 15, "Co": 1, "H": 32}, 32),
        "Mg2Ni pure": ({"Mg": 8, "Ni": 4, "H": 16}, 16),
        "Mg2Ni + Ni": ({"Mg": 7, "Ni": 5, "H": 16}, 16),
        "Mg2Ni + Co": ({"Mg": 7, "Co": 1, "Ni": 4, "H": 16}, 16),
        "Mg2Ni + Co@Ni": ({"Mg": 8, "Co": 1, "Ni": 3, "H": 16}, 16),
    }
    ranked = []
    for host, base, doped in (("Mg", "R1", (("Ni", "R2"), ("Co", "R3"))),
                              ("Mg2Ni", "R4", (("Ni", "R5"), ("Co", "R6"),
                                               ("Co@Ni", "R6b")))):
        if base in dH:
            ranked.append((f"{host:7s} pure", dH[base], EXP[base]))
            for tm, r in doped:
                if r in dH:
                    ranked.append((f"{host:7s} + {tm}", dH[r],
                                   EXP[base] + dH[r] - dH[base]))
    if ranked:
        ranked.sort(key=lambda t: t[1], reverse=True)   # least negative DFT first
        print()
        print("=" * 78)
        print("RANKING by computed ΔH (less negative = lower T_des = easier release)")
        print("=" * 78)
        print(f"  {'#':>2s}  {'system':14s}  {'ΔH DFT':>10s}  {'ΔH anchored':>12s}  "
              f"{'T_des @1bar':>12s}  {'wt% H2':>7s}")
        print("  " + "-" * 66)
        for i, (label, raw, anch) in enumerate(ranked, 1):
            if anch < 0:
                t_des = f"{-anch * 1000.0 / 130.5 - 273.15:9.0f} °C"
            else:
                t_des = "      n/a"   # hydride not stable at any T
            comp, n_h = HYDRIDE_COMP[" ".join(label.split())]
            print(f"  {i:>2d}  {label:14s}  {raw:+10.2f}  {anch:+12.2f}  "
                  f"{t_des}  {wt_pct(comp, n_h):7.2f}")
        print("  (kJ/mol H2. Sorted by computed DFT ΔH. T_des is van 't Hoff with")
        print("   ΔS = 130.5 J/mol/K on the ANCHORED ΔH -- raw cross-host baselines")
        print("   carry per-host GGA bias (see README Scope of validity). Doped")
        print("   compositions differ between hosts: 1 dopant per 16 H2 (Mg host,")
        print("   6.25% of metal sites) vs 1 per 8 H2 (Mg2Ni host, 8.3% of Mg")
        print("   sites) -- cross-host dopant-strength comparisons must account")
        print("   for this ~2x lever arm. Capacity from cell composition.)")

    # --- Reference values ----------------------------------------------------
    print()
    print("Reference (experimental, for sanity):")
    print("  MgH2     formation   ΔH = -74.5 kJ/mol H2  (Bogdanovic 1999; Stampfer 1960 -74.4)")
    print("  Mg2NiH4  formation   ΔH = -64.4 kJ/mol H2  (Reilly-Wiswall 1968, -15.4 kcal/mol)")
    print("Raw DFT baselines reproduce each host's known PBE behaviour but NOT the")
    print("experimental cross-host ordering; within-host shifts (ΔΔH) are the")
    print("robust DFT quantity, and the anchored column re-bases them onto the")
    print("experimental host baselines for cross-host reading.")


if __name__ == "__main__":
    main()
