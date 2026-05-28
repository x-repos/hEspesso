# Magnesium Hydride Catalyst Destabilization Study

This repository contains the complete computational workflow used to verify the thermodynamic destabilization of the $\text{Mg}_2\text{NiH}_4$ hydrogen storage lattice via the addition of $\text{Nb}_2\text{O}_5$ and $\text{Fe}$ catalysts.

## Contents

- **`*.in`**: Quantum ESPRESSO input files for all crystal structures (pristine, doped, and bulk oxides) with DFT-D3 and Hubbard U physics embedded.
- **`*.out`**: The raw computational logs from the GPU containing the self-consistent field iterations and final total energies.
- **`run_all.sh`**: Bash script to sequentially execute the inputs through Quantum ESPRESSO.
- **`enthalpy.py`**: Python script that parses the `.out` files, calculates the formation enthalpies, and evaluates the reaction thermodynamics.
- **`render_crystals.py`**: Python script using `ASE` (Atomic Simulation Environment) to generate the ball-and-stick representations of the supercells.
- **`report.tex` / `.pdf`**: The final LaTeX report documenting the methodology, results, and structural figures.

## How to Run

1.  **Run the DFT Simulations:**
    Ensure `pw.x` from Quantum ESPRESSO is installed and linked.
    ```bash
    bash run_all.sh
    ```
2.  **Evaluate Thermodynamics:**
    ```bash
    python3 enthalpy.py
    ```
3.  **Generate Figures:**
    ```bash
    python3 render_crystals.py
    ```
4.  **Compile the Report:**
    ```bash
    cd tex
    pdflatex report.tex
    ```
