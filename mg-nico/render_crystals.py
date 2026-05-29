"""
Render ball-and-stick PNGs of the relaxed crystal structures into figures/.

Reads the final geometry from each outputs/*.out and writes a single PNG
per structure using matplotlib (ase.io.write backend='matplotlib').

By default produces 7 figures (mg, mgh2, mgni, mgco, mgh2ni, mgh2co, h2),
but only the 3 hydrides are wired into the LaTeX report.
"""

from pathlib import Path
from ase.io import read, write
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
OUTPUTS = HERE / "outputs"
FIGDIR  = HERE / "figures"
FIGDIR.mkdir(exist_ok=True)

# species -> (descriptive title, rotation for camera)
JOBS = [
    ("mg",     "Mg HCP",                        "-72x,-12y,0z"),
    ("h2",     "H$_2$ molecule",                "0x,0y,0z"),
    ("mgh2",   "Rutile MgH$_2$",                "-75x,-10y,0z"),
    ("mgni",   "Mg$_{14}$Ni$_2$ (Ni 12.5%)",    "-72x,-12y,0z"),
    ("mgco",   "Mg$_{14}$Co$_2$ (Co 12.5%)",    "-72x,-12y,0z"),
    ("mgh2ni", "Mg$_{14}$Ni$_2$H$_{32}$ (Ni 12.5%)", "-75x,-10y,0z"),
    ("mgh2co", "Mg$_{14}$Co$_2$H$_{32}$ (Co 12.5%)", "-75x,-10y,0z"),
]

# Element colours (CPK-ish, with Ni/Co distinct)
COLORS = {
    "Mg": "#9bd86d",
    "H":  "#f8f8f8",
    "Ni": "#3c8bd8",
    "Co": "#d84f8e",
}


def render(name, title, rotation):
    """Render a single .out -> figures/<name>.png ball-and-stick."""
    src = OUTPUTS / f"{name}.out"
    if not src.exists():
        print(f"  skip {name}: {src} missing")
        return
    atoms = read(src, format="espresso-out", index=-1)
    fig, ax = plt.subplots(figsize=(4, 4), dpi=200)
    colors = [COLORS.get(sp, "#aaaaaa") for sp in atoms.get_chemical_symbols()]
    radii  = [0.55 if sp == "H" else 0.95 for sp in atoms.get_chemical_symbols()]
    from ase.visualize.plot import plot_atoms
    plot_atoms(
        atoms, ax,
        rotation=rotation,
        radii=radii,
        colors=colors,
        show_unit_cell=2,
    )
    ax.set_axis_off()
    ax.set_title(title, fontsize=14, pad=30)
    out = FIGDIR / f"{name}.png"
    fig.savefig(out, bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)
    print(f"  wrote {out.relative_to(HERE)}")


def render_legend():
    """A horizontal strip of CPK-coloured circles + element labels."""
    elements = [("Mg", 0.95), ("H", 0.55), ("Ni", 0.95), ("Co", 0.95)]
    fig, ax = plt.subplots(figsize=(4.0, 0.7), dpi=200)
    for i, (sp, r) in enumerate(elements):
        ax.scatter(i, 0, s=(r * 220)**1.1, c=COLORS[sp],
                   edgecolors="black", linewidths=0.7, zorder=3)
        ax.annotate(sp, (i, -0.55), ha="center", va="top", fontsize=11)
    ax.set_xlim(-0.6, len(elements) - 0.4)
    ax.set_ylim(-1.0, 0.6)
    ax.set_aspect("equal")
    ax.set_axis_off()
    out = FIGDIR / "legend.png"
    fig.savefig(out, bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)
    print(f"  wrote {out.relative_to(HERE)}")


def render_combined(out_name="mgnico_nolegend.png",
                    panels=("mgh2", "mgh2ni", "mgh2co"),
                    width_ratios=(0.55, 1.0, 1.0)):
    """Render several structures side by side in one figure, with titles
    but NO legend. width_ratios shrinks individual panels (the small
    rutile MgH2 cell gets a narrower panel by default)."""
    from ase.visualize.plot import plot_atoms
    titlemap = {n: t for n, t, _ in JOBS}
    rotmap   = {n: r for n, _, r in JOBS}
    if width_ratios is None or len(width_ratios) != len(panels):
        width_ratios = [1.0] * len(panels)
    fig, axes = plt.subplots(1, len(panels), figsize=(4 * sum(width_ratios), 4),
                             dpi=200, gridspec_kw={"width_ratios": list(width_ratios)})
    if len(panels) == 1:
        axes = [axes]
    for ax, name in zip(axes, panels):
        src = OUTPUTS / f"{name}.out"
        if not src.exists():
            print(f"  skip {name}: {src} missing")
            ax.set_axis_off()
            continue
        atoms = read(src, format="espresso-out", index=-1)
        colors = [COLORS.get(sp, "#aaaaaa") for sp in atoms.get_chemical_symbols()]
        radii  = [0.55 if sp == "H" else 0.95 for sp in atoms.get_chemical_symbols()]
        plot_atoms(atoms, ax, rotation=rotmap[name], radii=radii,
                   colors=colors, show_unit_cell=2)
        ax.set_axis_off()
        # Pin the (aspect-shrunk) axes box to the top of its cell so every
        # panel's title sits on the same line.
        ax.set_anchor("N")
        ax.set_title(titlemap[name], fontsize=14, pad=20)
    out = FIGDIR / out_name
    fig.savefig(out, bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)
    print(f"  wrote {out.relative_to(HERE)}")


def main():
    for name, title, rot in JOBS:
        render(name, title, rot)
    render_legend()
    render_combined()


if __name__ == "__main__":
    main()
