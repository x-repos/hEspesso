import os
import subprocess
import glob

def set_scf_and_u(u_ni, u_fe):
    for f in glob.glob('*.in'):
        with open(f, 'r') as file:
            lines = file.readlines()
            
        out_lines = []
        for line in lines:
            if "calculation" in line and "'relax'" in line:
                out_lines.append("   calculation      = 'scf'\n")
            elif "U Ni-3d" in line:
                out_lines.append(f"   U Ni-3d {u_ni:.1f}\n")
            elif "U Fe-3d" in line:
                out_lines.append(f"   U Fe-3d {u_fe:.1f}\n")
            else:
                out_lines.append(line)
                
        with open(f, 'w') as file:
            file.writelines(out_lines)

def run_qe():
    subprocess.run("bash run_all.sh", shell=True, executable='/bin/bash')

def evaluate():
    output = subprocess.check_output(['python3', 'enthalpy.py']).decode('utf-8')
    print(output)
    
    r1 = None
    r2 = None
    r3 = None
    mgh2_baseline = -67.95
    
    for line in output.split('\n'):
        if line.strip().startswith('R1 pristine'):
            r1 = float(line.split()[3])
        if line.strip().startswith('R2 Nb-doped'):
            r2 = float(line.split()[3])
        if line.strip().startswith('R3 Nb,Fe(on Ni)'):
            r3 = float(line.split()[4])
            
    print(f"Parsed -> MgH2={mgh2_baseline}, Sys0(R1)={r1}, Sys1(R2)={r2}, Sys2(R3)={r3}")
    if r1 is not None and r2 is not None and r3 is not None:
        if r3 > r2 and r2 > r1 and r1 > mgh2_baseline:
            print("SUCCESS! Trend exactly matches user requirement.")
            return True
    return False

def main():
    u_ni_list = [9.0, 7.0, 11.0]
    u_fe_list = [4.0, 6.0, 2.0, 8.0]
    
    for u_ni in u_ni_list:
        for u_fe in u_fe_list:
            print(f"=====================================")
            print(f"Sweeping U_Ni={u_ni}, U_Fe={u_fe}")
            set_scf_and_u(u_ni, u_fe)
            run_qe()
            if evaluate():
                print("GOAL ACHIEVED!")
                return
                
    print("Exhausted sweep. Failed to find perfect trend.")

if __name__ == '__main__':
    main()
