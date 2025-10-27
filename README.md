# molplane - Molecular Planarity Analyzer - script series for Gaussian 16 output analysis (n.2)
**Python** tool for automated molecular planarity analysis from Gaussian 16 output files.
Calculates molecular dimensions, planarity (best-fit plane distances), and Hirshfeld charge/spin sums for selected atoms.

## Features

- **3D Interactive GUI**: Visualize molecules and select atoms in real-time with matplotlib 3D rendering
- **Best-Fit Plane Calculation**: 
  - Plane determined via SVD (Singular Value Decomposition) on selected atoms
  - Planarity measured as average absolute distance from the plane
- **Molecular Dimensions**: Calculate extensions above and below the plane along the normal vector
- **Hirshfeld Analysis**: Sum Q-H charges, S-H spin densities, and Q-CM5 charges for the same atoms used in plane fitting
- **Multi-format Export**: CSV files with all metrics + XYZ files with optimized geometries

## Input Requirements
- Gaussian 16 output files (.out, .log) from **geometry optimization** calculations
- Must contain `Geom=AllCheck` section with optimized coordinates
- Hirshfeld charge section must be present for charge analysis (pop=hirshfeld)

## How to run it on Windows
Detailed instructions on how to run the script are reported at the beginning of the script.

## Output files
**XYZ file**

**CSV file** having as columns:
- **Filename**: Input Gaussian output file name
- **AvgDistance** (Å): Average distance of selected atoms from best-fit plane  
- **DimensionAbove** (Å): Maximum molecular extension above the plane
- **DimensionBelow** (Å): Maximum molecular extension below the plane
- **sum_Spin-Hirshfeld**: Sum of Hirshfeld spin densities (S-H) for selected atoms
- **sum_Charge-Hirshfeld**: Sum of Hirshfeld charges (Q-H) for selected atoms  
- **sum_Charge-CM5**: Sum of CM5 charges (Q-CM5) for selected atoms

## Citation
J.G. Vitillo, molplane, https://github.com/jennygvitillo/molplane (2025).




