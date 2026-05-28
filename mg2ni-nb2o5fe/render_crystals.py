from ase.io import read, write
import matplotlib.pyplot as plt
from ase.visualize.plot import plot_atoms

# --- CONFIGURATION ---
# Change this string to rotate the crystal structures in the output images.
# Format is 'Xx,Yy,Zz' in degrees (e.g., '30x,-70y,15z').
VIEW_ROTATION = '15x,-30y,0z'
# ---------------------

import numpy as np

def draw_axes_triad(ax, rotation_str, origin_frac=(0.1, 0.1), length_frac=0.1):
    R = np.eye(3)
    if rotation_str:
        for r in rotation_str.split(','):
            axis = r[-1]
            angle = float(r[:-1]) * np.pi / 180.0
            c, s = np.cos(angle), np.sin(angle)
            if axis == 'x': r_mat = np.array([[1,0,0], [0,c,s], [0,-s,c]])
            elif axis == 'y': r_mat = np.array([[c,0,-s], [0,1,0], [s,0,c]])
            elif axis == 'z': r_mat = np.array([[c,s,0], [-s,c,0], [0,0,1]])
            R = np.dot(R, r_mat)
            
    # ASE's projection aligns Z out of the screen. We project X,Y,Z vectors
    # We will place the triad dynamically
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    origin = np.array([xlim[0] + origin_frac[0]*(xlim[1]-xlim[0]), ylim[0] + origin_frac[1]*(ylim[1]-ylim[0])])
    length = length_frac * (xlim[1] - xlim[0])
    
    # Project basis vectors
    basis = np.eye(3) * length
    proj = np.dot(basis, R)[:, :2]  # Keep only X and Y coordinates on screen
    
    colors = ['r', 'g', 'b']
    labels = ['X', 'Y', 'Z']
    for i in range(3):
        ax.annotate('', xy=origin + proj[i], xytext=origin,
                    arrowprops=dict(arrowstyle="->", color=colors[i], lw=2))
        ax.text(origin[0] + proj[i,0]*1.2, origin[1] + proj[i,1]*1.2, labels[i],
                color=colors[i], fontsize=12, fontweight='bold', ha='center', va='center')

def get_custom_colors(atoms):
    from ase.data.colors import jmol_colors
    colors_array = []
    for atom in atoms:
        if atom.symbol == 'Mg':
            colors_array.append((1.0, 0.7, 0.2)) # Orange for Mg
        else:
            colors_array.append(jmol_colors[atom.number])
    return np.array(colors_array)

def render(filepath, outpath, title, format_type, rotation=VIEW_ROTATION):
    atoms = read(filepath, format=format_type)
    
    fig, ax = plt.subplots(figsize=(6, 6))
    custom_colors = get_custom_colors(atoms)
    plot_atoms(atoms, ax, radii=0.8, rotation=rotation, colors=custom_colors)
    ax.set_axis_off()
    
    # Draw the coordinate axes
    draw_axes_triad(ax, rotation)
    ax.set_title(title, fontsize=16, pad=15)
    
    plt.tight_layout()
    plt.savefig(outpath, dpi=300, bbox_inches='tight', transparent=True)
    plt.close()
    print(f"Rendered {outpath}")

def render_all_in_one(outpath, rotation=VIEW_ROTATION):
    files = [
        ('mg2nih4-28.pwi', '(a) Pristine Mg$_2$NiH$_4$', 'espresso-in'),
        ('mg2nih4_Nb.in', '(b) Mg$_2$NiH$_4$:Nb', 'espresso-in'),
        ('mg2nih4_NbFe.in', '(c) Mg$_2$NiH$_4$:Nb,Fe', 'espresso-in')
    ]
    
    fig, axes = plt.subplots(2, 2, figsize=(10, 10))
            
    unique_symbols = set()
    
    # Render the 3 crystal structures
    for i, (filepath, title, format_type) in enumerate(files):
        ax = axes.flat[i]
        atoms = read(filepath, format=format_type)
        unique_symbols.update(atoms.get_chemical_symbols())
        custom_colors = get_custom_colors(atoms)
        plot_atoms(atoms, ax, radii=0.8, rotation=rotation, colors=custom_colors)
        ax.set_axis_off()
        ax.set_title(title, fontsize=18, pad=2)
        
    # Use the bottom right (4th) subplot for axes triad and legend
    ax_br = axes.flat[3]
    ax_br.set_axis_off()
    ax_br.set_aspect('equal')
    ax_br.set_xlim(0, 1)
    ax_br.set_ylim(0, 1)
    
    # Draw triad nicely in the left-center of this quadrant
    draw_axes_triad(ax_br, rotation, origin_frac=(0.20, 0.5), length_frac=0.2)
    
    # Add legend to this quadrant
    from ase.data.colors import jmol_colors
    from ase.data import atomic_numbers
    from matplotlib.lines import Line2D
    
    handles = []
    for sym in sorted(unique_symbols):
        if sym == 'Mg':
            c = (1.0, 0.7, 0.2)
        else:
            c = jmol_colors[atomic_numbers[sym]]
        handles.append(Line2D([0], [0], marker='o', color='w', 
               markerfacecolor=c, 
               markeredgecolor='k', markersize=14, label=sym))

    # Place legend in the right-center of the bottom right plot
    ax_br.legend(handles=handles, loc='center right', ncol=1, 
                 fontsize=18, frameon=False, bbox_to_anchor=(0.95, 0.5))
        
    # Squeeze the plots together heavily vertically
    plt.subplots_adjust(hspace=-0.45)
    
    plt.savefig(outpath, dpi=300, bbox_inches='tight', pad_inches=0.01, transparent=True)
    plt.close()
    print(f"Rendered {outpath}")


if __name__ == '__main__':
    # System 0: Pristine Mg2NiH4
    render('mg2nih4-28.pwi', 'tex/struct_sys0.png', '(a) Pristine Mg$_2$NiH$_4$', 'espresso-in')
    # System 1: Doped Mg2NiH4 (Nb)
    render('mg2nih4_Nb.in', 'tex/struct_sys1.png', '(b) Mg$_2$NiH$_4$:Nb', 'espresso-in')

    # System 2: Doped Mg2NiH4 (Nb, Fe) - using .in because fast_sweep didn't finish an .out file completely for this one yet
    render('mg2nih4_NbFe.in', 'tex/struct_sys2.png', '(c) Mg$_2$NiH$_4$:Nb,Fe', 'espresso-in')

    # Combined figure
    render_all_in_one('tex/struct_combined.png')
