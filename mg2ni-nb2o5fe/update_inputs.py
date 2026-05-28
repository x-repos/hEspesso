import glob
import os

def update_file(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
        
    out_lines = []
    in_system = False
    has_ni = False
    has_fe = False
    
    for line in lines:
        if line.strip().upper() == '&SYSTEM':
            in_system = True
            out_lines.append(line)
            continue
            
        if in_system:
            if 'vdw_corr' in line.lower() or 'input_dft' in line.lower():
                continue # Strip old params
            if line.strip() == '/':
                in_system = False
                out_lines.append("   vdw_corr         = 'dft-d3'\n")
        
        # Check for Ni or Fe in ATOMIC_SPECIES block
        if line.strip().startswith('Ni '):
            has_ni = True
        if line.strip().startswith('Fe '):
            has_fe = True
            
        out_lines.append(line)
        
    content = "".join(out_lines)
    
    # Strip existing HUBBARD block if it exists
    if "HUBBARD" in content:
        # Just a naive strip for this specific workspace
        content = content.split("HUBBARD")[0]
        
    if has_ni or has_fe:
        content += "\nHUBBARD (ortho-atomic)\n"
        if has_ni:
            content += "   U Ni-3d 9.0\n"
        if has_fe:
            content += "   U Fe-3d 4.0\n"
            
    with open(filename, 'w') as f:
        f.write(content)
        
    print(f"Updated {filename}: Ni={has_ni}, Fe={has_fe}")

if __name__ == '__main__':
    here = os.path.dirname(os.path.abspath(__file__))
    for f in sorted(glob.glob(os.path.join(here, 'inputs', '*.in'))):
        update_file(f)
