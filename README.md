# Grotthuss-type-proton-transport
Grotthuss-type proton transport through hydrogen molecules H2 in hydrogen hydrate under extreme conditions.

This repository contains the MATLAB analysis scripts used to identify and characterize the Grotthuss-type proton relay mechanism in compressed hydrogen-water mixtures (C3 phase) using First-Principles Molecular Dynamics (FPMD) trajectories (VASP `XDATCAR` format).
## ⚙️ Requirements

- **MATLAB** 
- **Statistics and Machine Learning Toolbox** (for `histcounts`, `fspecial`, `imfilter`)
- **Image Processing Toolbox** (for `fspecial`, `imfilter`)
## 🚀 How to Use

1. **Prepare the trajectory**: Ensure you have a VASP `XDATCAR` file in the same directory. Update the `filename` (or `fpath`) variable at the top of each script to match your file name.
2. **Set system parameters**: Modify `n_H` and `n_O` to match the number of atoms in your system.
3. **Run the script**: Execute each script sequentially or independently in MATLAB.
   - For example, run `cross_section_total.m` to generate `Anisotropic_Cross_Section.png` and `.csv` files.
4. **Adjust thresholds**: The scripts use a **hysteresis criterion** (formation/dissociation thresholds) for defining \(H_2\) and \(H_3^+\). You can adjust `r_HH_on`, `r_HH_off`, `r_H3_on`, and `r_H3_off` at the top of each script.
