#!/usr/bin/env python3
"""
Molecular Planarity Analyzer for Gaussian Output Files
======================================================

Description:
------------
This script analyzes molecular planarity and charge properties from Gaussian 16 
output files. It calculates best-fit planes, molecular dimensions, and sums 
Hirshfeld charges for selected atom groups. Perfect for studying planar 
conjugated systems and charge distribution in molecular cores.

Key Features:
-------------
- Interactive 3D molecular visualization and atom selection
- Best-fit plane calculation using SVD (Singular Value Decomposition)
- Molecular dimension analysis along plane normal direction
- Automated extraction and summation of Hirshfeld charges (Q-H, S-H, Q-CM5)
- Export results to CSV and XYZ formats
- User-friendly graphical interface

Workflow:
---------
1. Input: Gaussian 16 output file (.out or .log) from geometry optimization
2. 3D Visualization: Interactive molecular structure display
3. Atom Selection: Choose atoms for plane calculation (≥3 atoms)
4. Analysis: Automatic calculation of planarity, dimensions, and charge sums
5. Output: CSV file with metrics + XYZ file with optimized geometry

Requirements:
-------------
- Python 3.6+ (I am using Windows: so I have installed PyCharm and this program 
  deals with all the Python stuff, very convenient)
- NumPy (for mathematical calculations)
- Matplotlib (for 3D visualization)

Installation for Dummies:
-------------------------
0. If you don't have it, download and install PyCharm Community Edition: 
   it makes running Python scripts much easier.

1. Create and open the project in PyCharm.
   - Open PyCharm and create a **new project**. You can name it molplane.  
   - On the top left, you'll see something that looks like a folder.  
   - **Right-click** on it → choose **New → Python File**.  
   - Name it `molplane.py` and press **Enter**.  
   - Copy and paste the entire script into the text area that appears.  
   - Go to **File → Save All**.

2. Install required packages:
   - Open Terminal in PyCharm (small icon on the bottom left)
   - Type these commands one by one:
     pip install numpy
     pip install matplotlib
   - Wait for installation to complete
  Just so you know, this is necessary only the first time you are using this script on your laptop.

3. To run the script:
   - Click the green "Play" button above the text area
   - OR right-click in the script and select "Run 'molplane'"
   - Follow the on-screen instructions to select your Gaussian file

Usage:
------
1. Run the script
2. Select your Gaussian output file (.out or .log) using the file dialog
3. In the 3D visualization window:
   - Rotate: Click and drag to rotate the molecule
   - Select atoms: Click on atoms in the list (≥3 atoms required)
   - Selected atoms will be highlighted in the 3D view
4. Click "Confirm Selection" to proceed with analysis
5. View results in the pop-up window and check the output files


Input Requirements:
-------------------
- Gaussian 16 output file from geometry optimization calculation
- Must contain optimized coordinates (Geom=AllCheck section)
- Must include frequency calculation section in the route (even if not completed)
- For charge analysis: must include Hirshfeld charges section (from frequency calculation)

Note: The frequency calculation doesn't need to be completed - the script only 
requires that the frequency job was initiated and the optimized geometry was 
printed. Failed frequency calculations that still output the geometry are perfectly acceptable.


Output:
-------
- CSV file (with .txt extension): Contains all calculated metrics
  - Filename, AvgDistance, DimensionAbove, DimensionBelow
  - sum_Spin-Hirshfeld, sum_Charge-Hirshfeld, sum_Charge-CM5
- XYZ file: Optimized molecular geometry for further analysis
- Files are saved in "Planarity_Analysis_Results" folder


Troubleshooting:
----------------
- "No coordinates found": Check your file contains optimized geometry
- "No Hirshfeld charges": File missing frequency calculation section
- "Module not found": Run pip install for missing packages
- 3D view not displaying: Ensure matplotlib is properly installed
- Atom selection not working: Make sure to select ≥3 atoms
- "Frequency section missing": Ensure your Gaussian calculation includes 'freq' 
  in the route section, even if the frequency calculation didn't complete

  
Technical Details:
------------------
- Plane fitting: Singular Value Decomposition (SVD) on selected atom coordinates
- Planarity metric: Average absolute distance from best-fit plane
- Charge summation: Sum of Hirshfeld properties for selected atoms only
- Normalization: All distances in Angstroms (Å)


Author:
---------
Jenny G. Vitillo (University of Insubria) and DeepSeek 2024


License:
---------
MIT


Citation:
---------
J.G. Vitillo, molplane, https://github.com/jennygvitillo/molplane (2025).
Additionally, please cite relevant Gaussian and Hirshfeld methodology references.
"""


import os
import sys
import csv
import numpy as np
import tkinter as tk
from tkinter import filedialog, messagebox, ttk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
from mpl_toolkits.mplot3d import Axes3D

# Color map for elements (CPK colors)
ELEMENT_COLORS = {
    'H': 'white', 'C': 'gray', 'N': 'blue', 'O': 'red',
    'S': 'yellow', 'Si': 'orange', 'P': 'orange', 'F': 'green',
    'Cl': 'green', 'Br': 'darkred', 'I': 'darkviolet'
}


def check_dependencies():
    """Check for required Python packages"""
    try:
        import numpy as np
        import matplotlib
    except ImportError:
        sys.exit("Error: Required packages not found. Please install with: pip install numpy matplotlib")


class AtomSelector(tk.Toplevel):
    """GUI window for atom selection with 3D visualization"""

    def __init__(self, parent, elements, coordinates):
        super().__init__(parent)
        self.title("Select Atoms for Plane Calculation (≥3 atoms)")
        self.elements = elements
        self.coords = np.array(coordinates)
        self.selected_indices = []

        # Setup GUI layout
        self.listbox_frame = tk.Frame(self)
        self.visualization_frame = tk.Frame(self)

        # Listbox for atom selection
        self.listbox = tk.Listbox(self.listbox_frame, selectmode=tk.MULTIPLE,
                                  width=60, height=25, font=('Courier', 10))
        scrollbar = tk.Scrollbar(self.listbox_frame, orient="vertical")
        scrollbar.config(command=self.listbox.yview)
        self.listbox.config(yscrollcommand=scrollbar.set)

        # Populate listbox
        for i, (elem, coord) in enumerate(zip(elements, coordinates)):
            self.listbox.insert(tk.END,
                                f"Atom {i + 1:3d}: {elem:2s} | X: {coord[0]:10.6f} | Y: {coord[1]:10.6f} | Z: {coord[2]:10.6f}")

        # 3D Visualization setup
        self.fig = Figure(figsize=(6, 6), dpi=100)
        self.ax = self.fig.add_subplot(111, projection='3d')
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.visualization_frame)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        # Initial plot
        self.plot_structure()

        # Buttons
        btn_frame = tk.Frame(self)
        tk.Button(btn_frame, text="Confirm Selection", command=self.on_confirm,
                  bg="#4CAF50", fg="white").pack(side=tk.LEFT, padx=5, pady=5)
        tk.Button(btn_frame, text="Cancel", command=self.on_cancel,
                  bg="#f44336", fg="white").pack(side=tk.LEFT, padx=5, pady=5)

        # Layout
        self.listbox.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        self.listbox_frame.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        self.visualization_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)
        btn_frame.pack(fill=tk.X)

        # Bind selection event
        self.listbox.bind('<<ListboxSelect>>', self.on_selection_change)

        self.protocol("WM_DELETE_WINDOW", self.on_cancel)
        self.grab_set()

    def plot_structure(self, highlight_indices=None):
        """Plot the molecular structure with proper 3D aspect ratio"""
        self.ax.clear()

        # Plot atoms
        for i, (elem, coord) in enumerate(zip(self.elements, self.coords)):
            color = ELEMENT_COLORS.get(elem, 'purple')
            alpha = 0.3 if highlight_indices and i not in highlight_indices else 1.0
            size = 100 if highlight_indices and i in highlight_indices else 50
            self.ax.scatter(*coord, color=color, s=size, alpha=alpha, edgecolors='black')

        # Fix aspect ratio
        all_coords = np.vstack([self.coords, [0, 0, 0]])  # Add origin
        max_range = np.ptp(all_coords, axis=0).max() / 2

        mid_x = (all_coords[:, 0].max() + all_coords[:, 0].min()) * 0.5
        mid_y = (all_coords[:, 1].max() + all_coords[:, 1].min()) * 0.5
        mid_z = (all_coords[:, 2].max() + all_coords[:, 2].min()) * 0.5

        self.ax.set_xlim(mid_x - max_range, mid_x + max_range)
        self.ax.set_ylim(mid_y - max_range, mid_y + max_range)
        self.ax.set_zlim(mid_z - max_range, mid_z + max_range)

        self.ax.set_box_aspect([1, 1, 1])
        self.ax.set_xlabel('X (Å)')
        self.ax.set_ylabel('Y (Å)')
        self.ax.set_zlabel('Z (Å)')

        self.canvas.draw()

    def on_selection_change(self, event):
        """Update visualization when selection changes"""
        selected = self.listbox.curselection()
        self.selected_indices = list(selected)
        self.plot_structure(highlight_indices=self.selected_indices)

    def on_confirm(self):
        if len(self.selected_indices) < 3:
            messagebox.showerror("Selection Error", "Please select at least 3 atoms to define a plane")
            return
        self.destroy()

    def on_cancel(self):
        self.selected_indices = []
        self.destroy()


def parse_gaussian_file(filepath):
    """Extract optimized coordinates from Gaussian output file"""
    elements = []
    coordinates = []
    capture = False

    try:
        with open(filepath, 'r') as f:
            for line in f:
                line = line.strip()
                if "Geom=AllCheck" in line:
                    capture = True
                    continue
                if capture:
                    if "Recover connectivity data from disk." in line:
                        break
                    if line and line[0].isalpha():
                        parts = [p.strip() for p in line.split(',')]
                        if len(parts) >= 5:
                            try:
                                elements.append(parts[0])
                                coordinates.append([float(parts[2]), float(parts[3]), float(parts[4])])
                            except (ValueError, IndexError) as e:
                                print(f"Skipping malformed line: {line} (Error: {str(e)})")
                                continue
    except Exception as e:
        messagebox.showerror("File Error", f"Error reading {filepath}:\n{str(e)}")
        return [], []

    if not elements:
        messagebox.showerror("Parsing Error", "No valid coordinates found in the selected file")

    return elements, np.array(coordinates)


def parse_hirshfeld_charges(filepath, selected_indices):
    """Parse Hirshfeld charges and CM5 charges from Gaussian output file"""
    hirshfeld_data = {
        'Q-H': [],
        'S-H': [],
        'Q-CM5': []
    }
    capture = False
    line_count = 0

    try:
        with open(filepath, 'r') as f:
            for line in f:
                if "Hirshfeld charges, spin densities, dipoles, and CM5 charges" in line:
                    line_count += 1
                    if line_count == 2:  # Second occurrence
                        capture = True
                        next(f)  # Skip header line
                        continue

                if capture:
                    if line.strip() == "":  # Stop at empty line
                        break

                    # Parse data line (format: index element Q-H S-H Dx Dy Dz Q-CM5)
                    parts = line.split()
                    if len(parts) >= 7:
                        try:
                            atom_index = int(parts[0]) - 1  # Convert to 0-based index
                            if atom_index in selected_indices:
                                hirshfeld_data['Q-H'].append(float(parts[2]))
                                hirshfeld_data['S-H'].append(float(parts[3]))
                                hirshfeld_data['Q-CM5'].append(float(parts[7]))
                        except (ValueError, IndexError) as e:
                            print(f"Skipping malformed line: {line} (Error: {str(e)})")
                            continue

    except Exception as e:
        messagebox.showerror("File Error", f"Error reading Hirshfeld charges from {filepath}:\n{str(e)}")

    # Calculate sums
    sums = {
        'sum_Spin-Hirshfeld': sum(hirshfeld_data['S-H']),
        'sum_Charge-Hirshfeld': sum(hirshfeld_data['Q-H']),
        'sum_Charge-CM5': sum(hirshfeld_data['Q-CM5'])
    }

    return sums


def calculate_planarity(selected_coords):
    """Calculate best-fit plane and average distance"""
    centroid = np.mean(selected_coords, axis=0)
    centered = selected_coords - centroid
    _, _, vh = np.linalg.svd(centered)
    normal = vh[2]
    normal /= np.linalg.norm(normal)  # Normalize

    distances = np.abs(np.dot(centered, normal))
    avg_distance = np.mean(distances)
    return centroid, normal, avg_distance


def calculate_dimensions(centroid, normal, all_coords):
    """Calculate molecular dimensions along normal direction"""
    vectors = all_coords - centroid
    projections = np.dot(vectors, normal)
    return np.max(projections), np.abs(np.min(projections))


def save_results(output_dir, filename, elements, coords, avg_dist, above, below, hirshfeld_sums):
    """Save results to CSV file (with .txt extension for Excel compatibility)"""
    os.makedirs(output_dir, exist_ok=True)
    base_name = os.path.splitext(filename)[0]
    csv_path = os.path.join(output_dir, f"{base_name}_planarity.txt")  # .txt for Excel

    with open(csv_path, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)

        # Write headers
        writer.writerow([
            "Filename",
            "AvgDistance",
            "DimensionAbove",
            "DimensionBelow",
            "sum_Spin-Hirshfeld",
            "sum_Charge-Hirshfeld",
            "sum_Charge-CM5"
        ])

        # Write data row
        writer.writerow([
            filename,
            f"{avg_dist:.6f}",
            f"{above:.6f}",
            f"{below:.6f}",
            f"{hirshfeld_sums['sum_Spin-Hirshfeld']:.6f}",
            f"{hirshfeld_sums['sum_Charge-Hirshfeld']:.6f}",
            f"{hirshfeld_sums['sum_Charge-CM5']:.6f}"
        ])

    # Save XYZ file (unchanged)
    xyz_path = os.path.join(output_dir, f"{base_name}.xyz")
    with open(xyz_path, 'w') as xyzfile:
        xyzfile.write(f"{len(elements)}\n")
        xyzfile.write(f"Optimized structure from {filename}\n")
        for elem, coord in zip(elements, coords):
            xyzfile.write(f"{elem:2s} {coord[0]:12.8f} {coord[1]:12.8f} {coord[2]:12.8f}\n")

    return csv_path, xyz_path


def main():
    check_dependencies()

    root = tk.Tk()
    root.withdraw()

    # Ask user to select a Gaussian output file
    filepath = filedialog.askopenfilename(
        title="Select Gaussian Output File",
        filetypes=[("Gaussian files", "*.out *.log"), ("All files", "*.*")]
    )

    if not filepath:
        return

    filename = os.path.basename(filepath)
    output_dir = os.path.join(os.path.dirname(filepath), "Planarity_Analysis_Results")

    # Parse the file
    elements, coords = parse_gaussian_file(filepath)
    if not elements:
        return

    # Atom selection GUI
    selector = AtomSelector(root, elements, coords)
    root.wait_window(selector)

    if not selector.selected_indices:
        return

    # Calculate planarity and dimensions
    selected_coords = coords[selector.selected_indices]
    centroid, normal, avg_dist = calculate_planarity(selected_coords)
    above, below = calculate_dimensions(centroid, normal, coords)

    # Parse Hirshfeld charges for selected atoms
    hirshfeld_sums = parse_hirshfeld_charges(filepath, selector.selected_indices)

    # Save results (now including Hirshfeld sums)
    csv_path, xyz_path = save_results(output_dir, filename, elements, coords,
                                      avg_dist, above, below, hirshfeld_sums)

    # Update results message with CSV-style formatting
    result_msg = (
        "Filename,AvgDistance,DimensionAbove,DimensionBelow,sum_Spin-Hirshfeld,sum_Charge-Hirshfeld,sum_Charge-CM5\n"
        f"{filename},{avg_dist:.6f},{above:.6f},{below:.6f},"
        f"{hirshfeld_sums['sum_Spin-Hirshfeld']:.6f},"
        f"{hirshfeld_sums['sum_Charge-Hirshfeld']:.6f},"
        f"{hirshfeld_sums['sum_Charge-CM5']:.6f}\n\n"
        f"File saved to:\n{csv_path}\n{xyz_path}"
    )

    messagebox.showinfo("Analysis Results", result_msg)

if __name__ == "__main__":
    main()

