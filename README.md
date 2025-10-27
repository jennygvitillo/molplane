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

## How to run it
(on Windows, for dummies) Detailed instructions on how to run the script are reported at the beginning of the script.
for the other users: like all the other Python tools.

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

## Italiano
Strumento **Python** per analizzare file output di Gaussian 16 da calcoli di ottimizzazione geometria. Calcola planarità molecolare e somme di cariche di Hirshfeld. Interfaccia grafica per selezione interattiva degli atomi e salva i risultati in file CSV (con estensione TXT) e XYZ.

- **Calcolo planarità molecolare** - determina il piano migliore per atomi selezionati e calcola la distanza media
- **Dimensioni molecolari** - calcola l'estensione sopra e sotto il piano lungo la normale
- **Cariche di Hirshfeld** - somma cariche Q-H, spin S-H e cariche CM5 Q-CM5 per gli atomi selezionati
- **Visualizzazione 3D interattiva** - interfaccia grafica per selezione atomi con matplotlib
- **Esportazione multipla** - file CSV con tutte le metriche e file XYZ con geometrie ottimizzate
- **Elaborazione file singolo** - selezione interattiva del file .log/.out da analizzare
- **File di output** - salvati in directory "Planarity_Analysis_Results" con nome file_input_planarity.txt

**Dipendenze**: numpy, matplotlib
**Istruzioni dettagliate** su come usare lo script sono riportate nella parte iniziale dello script stesso (file: molplane.py).


