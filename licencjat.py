import os
import json
import urllib.request
import re
from pymol import cmd
from PyQt5 import QtWidgets, QtCore, QtGui
from matplotlib.colors import ListedColormap, BoundaryNorm
import numpy as np
from scipy.spatial import distance_matrix
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

def safe_int_resi(resi_str):
    """Extracts integer residue number from a string using regex."""
    match = re.search(r'(\d+)', resi_str.split()[-1])
    return int(match.group(1)) if match else 0

# --- COCOMAPS STYLE INTERACTION PLOT WINDOW ----
class InteractionPlotWindow(QtWidgets.QDialog):
    """Dialog window for displaying the 2D interaction heat map."""
    def __init__(self, extra_plot_data, data_dict, chain1_id, chain2_id):
        super().__init__()
        self.setWindowTitle(f"2D Interaction Map: {chain1_id} vs {chain2_id}")
        self.resize(900, 850)
        self.setWindowFlags(QtCore.Qt.WindowStaysOnTopHint)
        
        layout = QtWidgets.QVBoxLayout(self)
        self.canvas = InteractionMapCanvas(
            extra_plot_data['matrix'], 
            data_dict, 
            extra_plot_data['map1'], 
            extra_plot_data['map2'],
            extra_plot_data['labels1'], 
            extra_plot_data['labels2'],
            chain1_id,
            chain2_id
        )
        layout.addWidget(self.canvas)

class AdvancedInteractionWindow(QtWidgets.QDialog):
    """Window displaying detailed interaction tables (H-Bonds, Salt Bridges, etc.)."""
    def __init__(self, data_dict, bio_names=None):
        super().__init__()
        self.setWindowTitle("PDB Interaction Analysis (PLIP/HBAT Style)")
        self.resize(1000, 500)
        self.setWindowFlags(QtCore.Qt.WindowStaysOnTopHint)
        self.data_dict = data_dict
        
        layout = QtWidgets.QVBoxLayout(self)
        
        # Display molecule names from PDB metadata if available
        if bio_names:
            info_label = QtWidgets.QLabel(f"<b>PDB Database Info:</b><br>{bio_names}")
            info_label.setStyleSheet("background-color: #ecf0f1; padding: 10px; border-radius: 5px;")
            layout.addWidget(info_label)

        self.tabs = QtWidgets.QTabWidget()
        layout.addWidget(self.tabs)

        # Initialize tables for different interaction types
        headers_hbond = ["Residue A", "Chain A", "Atom A (Donor/Acc)", "Residue B", "Chain B", "Atom B (Donor/Acc)", "Distance [Å]", "Angle [°]"]
        self.tab_hbond = self.create_table(data_dict["hbonds"], headers_hbond)
        self.tabs.addTab(self.tab_hbond, "Hydrogen Bonds")

        headers_salt = ["Residue (+) / A", "Chain A", "Residue (-) / B", "Chain B", "Distance [Å]"]
        self.tab_salt = self.create_table(data_dict["salt"], headers_salt)
        self.tabs.addTab(self.tab_salt, "Salt Bridges")

        headers_hydro = ["Residue A", "Chain A", "Residue B", "Chain B", "Carbon Distance [Å]"]
        self.tab_hydro = self.create_table(data_dict["hydro"], headers_hydro)
        self.tabs.addTab(self.tab_hydro, "Hydrophobic Interactions")
        
        # Export functionality
        btn_export = QtWidgets.QPushButton("Export all tables to CSV")
        btn_export.setStyleSheet("background-color: #27ae60; color: white; font-weight: bold; padding: 10px; border-radius: 4px;")
        btn_export.clicked.connect(self.export_all_csv)
        layout.addWidget(btn_export)

    def create_table(self, data, headers):
        """Helper function to populate QTableWidget with analysis data."""
        table = QtWidgets.QTableWidget()
        table.setColumnCount(len(headers))
        table.setHorizontalHeaderLabels(headers)
        table.setAlternatingRowColors(True)
        table.setRowCount(len(data))
        for i, row in enumerate(data):
            for j, val in enumerate(row):
                item = QtWidgets.QTableWidgetItem(str(val))
                item.setTextAlignment(QtCore.Qt.AlignCenter)
                table.setItem(i, j, item)
        table.horizontalHeader().setSectionResizeMode(QtWidgets.QHeaderView.Stretch)
        return table

    def export_all_csv(self):
        """Saves interaction data to a CSV file."""
        path, _ = QtWidgets.QFileDialog.getSaveFileName(self, "Save Analysis", "", "CSV Files (*.csv)")
        if not path: return
        header = ["Type", "Res_A", "Chain_A", "Atom_A", "Res_B", "Chain_B", "Atom_B", "Distance", "Angle"]
        try:
            with open(path, 'w', encoding='utf-8') as f:
                f.write(",".join(header) + "\n")
                for cat, label in [("hbonds", "H-Bond"), ("salt", "Salt-Bridge"), ("hydro", "Hydrophobic")]:
                    for row in self.data_dict[cat]:
                        if cat == "hbonds":
                            f.write(f"{label}," + ",".join([str(x) for x in row]) + "\n")
                        else:
                            # Handling different column counts for non-hbond interactions
                            f.write(f"{label},{row[0]},{row[1]},N/A,{row[2]},{row[3]},N/A,{row[4]},N/A\n")
            QtWidgets.QMessageBox.information(self, "Success", "Data exported successfully!")
        except Exception as e:
            QtWidgets.QMessageBox.critical(self, "Error", f"Failed to export: {e}")

class InteractionMapCanvas(FigureCanvas):
    """Matplotlib canvas for rendering the interactive 2D contact map."""
    def __init__(self, dist_matrix, interactions, map1, map2, labels1, labels2, c1_id, c2_id):
        fig, (self.ax1, self.ax2) = plt.subplots(1,2, figsize=(9, 7), sharex=True, sharey=True)
        super().__init__(fig)
        self.fig = fig

        #  Define Distance Thresholds and Colors 
        bounds = [0, 3, 5, 8, 11] 
        colors = ['#e74c3c', '#f1c40f', '#3498db', '#ffffff']
        cmap = ListedColormap(colors)
        norm = BoundaryNorm(bounds, cmap.N)

        #  Draw Contact Matrix
        im = self.ax1.imshow(dist_matrix, cmap=cmap, norm=norm, origin='lower', aspect='equal')
        
        

        #  Labels and Titles
        self.ax1.set_ylabel(f"Chain {c1_id} ", fontsize=11, fontweight='bold')
        self.ax1.set_xlabel(f"Chain {c2_id} ", fontsize=11, fontweight='bold')
        self.ax1.set_title(f"Interactive Map: {c1_id} vs {c2_id}", pad=20, fontsize=13)
        
        # Legend Setup(ax1)
        dist_elements = [
            Patch(facecolor='#e74c3c', edgecolor='black', label='3 Å'),
            Patch(facecolor='#f1c40f', edgecolor='black', label='5 Å'),
            Patch(facecolor='#3498db', edgecolor='black', label='8 Å'),
            Patch(facecolor='#ffffff', edgecolor='black', label='11 Å')
        ]
        
        
        self.ax1.legend(handles=dist_elements, loc='upper left', bbox_to_anchor=(1.05, 1), 
                       title="Distance", title_fontproperties={'weight': 'bold'},
                       frameon=True, shadow=True)
        
        self.ax2.set_title(f"Interaction: {c1_id} vs {c2_id}", pad=20, fontsize=13)
        self.ax2.set_ylabel(f"Chain {c1_id} ", fontsize=11, fontweight='bold')
        self.ax2.set_xlabel(f"Chain {c2_id}", fontsize=11, fontweight='bold')
        
        self.ax2.set_xlim(-0.5, len(labels2) - 0.5)
        self.ax2.set_ylim(-0.5, len(labels1) - 0.5)
        self.ax2.set_aspect('equal')
        
        
        # Scale X-axis ticks (showing every ~20 residues for readability)
        step_x = max(1, len(labels2) // 20)
        xticks_pos = range(0, len(labels2), step_x)
        xticklabels = [re.search(r'\d+', labels2[i]).group() if re.search(r'\d+', labels2[i]) else "" for i in xticks_pos]
    
        
        # Scale Y-axis ticks
        step_y = max(1, len(labels1) // 20)
        yticks_pos = range(0, len(labels1), step_y)
        yticklabels = [re.search(r'\d+', labels1[i]).group() if re.search(r'\d+', labels1[i]) else "" for i in yticks_pos]
   

       
        # Grid Setup
        for ax in [self.ax1, self.ax2]:
         
         ax.set_xticks(np.arange(-0.5, len(labels2), 1), minor=True)
         ax.set_yticks(np.arange(-0.5, len(labels1), 1), minor=True)
         ax.grid(which='minor', color='#bdc3c7', linestyle='-', linewidth=0.4, alpha=0.5)
         ax.tick_params(which='minor', size=0) 

         ax.set_xticks(xticks_pos)
         ax.set_xticklabels(xticklabels, rotation=0, fontsize=9)
         ax.set_yticks(yticks_pos)
         ax.set_yticklabels(yticklabels, fontsize=9)
         self.ax2.tick_params(axis='y', labelleft=True)
        

        #  Initialize Interactive Tooltip
        self.annot = self.ax2.annotate("", xy=(0,0), xytext=(15, 15), textcoords="offset points",
                                      bbox=dict(boxstyle="round", fc="white", ec="#bdc3c7", alpha=0.95),
                                      arrowprops=dict(arrowstyle="->", color="black"))
        self.annot.set_visible(False)
        
        
        self.annot1 = self.ax1.annotate("", xy=(0,0), xytext=(15, 15), textcoords="offset points",
                                      bbox=dict(boxstyle="round", fc="white", ec="#bdc3c7", alpha=0.95),
                                      arrowprops=dict(arrowstyle="->", color="black"))
        self.annot1.set_visible(False)
        
        self.dist_matrix = dist_matrix
        self.labels1 = labels1
        self.labels2 = labels2
        
        self.scatter_plots = [] 

        #  Overlay Specific Interactions (H-Bonds, etc.)
        self.plot_interactions(interactions, map1, map2, labels1, labels2, c1_id, c2_id)

        #  Legend Setup(ax2)
        
        inter_elements = [
            Line2D([0], [0], marker='o', color='w', label='H-Bond', markerfacecolor='white', markeredgecolor='black', markersize=8),
            Line2D([0], [0], marker='s', color='w', label='Salt Bridge', markerfacecolor='magenta', markeredgecolor='black', markersize=8),
            Line2D([0], [0], marker='^', color='w', label='Hydrophobic', markerfacecolor='#2ecc71', markeredgecolor='black', markersize=8)
        ]
       
        
        self.ax2.legend(handles=inter_elements, loc='upper left', bbox_to_anchor=(1.05, 1), 
                       title="Interactions", title_fontproperties={'weight': 'bold'},
                       frameon=True, shadow=True)

        fig.tight_layout()
        self.fig.canvas.mpl_connect("motion_notify_event", self.hover)

    def plot_interactions(self, data, map1, map2, labels1, labels2, c1_id, c2_id):
        """Plots individual interaction markers over the distance matrix."""
        style = {
            "hbonds": {"c": "white", "marker": "o", "s": 80, "label": "H-Bond", "chain_idx": [1, 4]}, 
            "salt":   {"c": "magenta", "marker": "s", "s": 80, "label": "Salt Bridge", "chain_idx": [1, 3]},
            "hydro":  {"c": "#2ecc71", "marker": "^", "s": 80, "label": "Hydrophobic", "chain_idx": [1, 3]}
        }

        for itype, s_info in style.items():
            x_pts, y_pts, hover_texts = [], [], []
            
            for entry in data.get(itype, []):
                try:
                 
                    if itype == "hbonds":
                        res_a_id, chain_a, res_b_id, chain_b = entry[0], entry[1], entry[3], entry[4]
                    else:
                        res_a_id, chain_a, res_b_id, chain_b = entry[0], entry[1], entry[2], entry[3]

                  
                    r_a = re.search(r'(\d+[A-Z]?)', str(res_a_id)).group()
                    r_b = re.search(r'(\d+[A-Z]?)', str(res_b_id)).group()

                    y_idx, x_idx = None, None

                  
                    if chain_a == c1_id and r_a in map1 and chain_b == c2_id and r_b in map2:
                        y_idx, x_idx = map1[r_a], map2[r_b]
                    elif chain_a == c2_id and r_a in map2 and chain_b == c1_id and r_b in map1:
                        y_idx, x_idx = map1[r_b], map2[r_a]

                    if y_idx is not None and x_idx is not None:
                        y_pts.append(y_idx)
                        x_pts.append(x_idx)
                        hover_texts.append(
                            f"{labels1[y_idx]} ↔ {labels2[x_idx]}\nType: {s_info['label']}"
                        )
                except Exception as e:
                    continue

            if x_pts:
                sc = self.ax2.scatter(
                    x_pts, y_pts,
                    c=s_info["c"], marker=s_info["marker"], s=s_info["s"],
                    edgecolors="black", linewidths=0.8, zorder=10
                )
                self.scatter_plots.append((sc, hover_texts))

    def hover(self, event):
        """Handles mouse hover events for the interactive map."""
        vis = self.annot.get_visible()
        
        if event.inaxes ==self.ax1:
            x,y = int(round(event.xdata)), int(round(event.ydata))
            
            
            if 0 <= y < self.dist_matrix.shape[0] and 0 <= x < self.dist_matrix.shape[1]:
                dist_val = self.dist_matrix[y, x]
                
        
                self.annot1.xy = (event.xdata, event.ydata)
                text = f"{self.labels1[y]} ↔ {self.labels2[x]}\nDist: {dist_val:.2f} Å"
                self.annot1.set_text(text)
                
        
                self.annot1.get_bbox_patch().set_facecolor("#ecf0f1")
                
                self.annot1.set_visible(True)
                self.annot.set_visible(False) 
                self.fig.canvas.draw_idle()
                return
            
        elif event.inaxes == self.ax2:
            found=False
            
            for sc, hover_texts in self.scatter_plots:
                cont, ind = sc.contains(event)
                if cont:
                    point_idx = ind["ind"][0]
                    pos = sc.get_offsets()[point_idx]
                    
                    self.annot.xy = pos
                    self.annot.set_text(hover_texts[point_idx])
                    
                
                    self.annot.get_bbox_patch().set_facecolor(sc.get_facecolor()[0])
                    self.annot.get_bbox_patch().set_alpha(0.9)
                    
                    self.annot.set_visible(True)
                    found = True
                    break 
            
            if not found and self.annot.get_visible():
                self.annot.set_visible(False)
            
            self.fig.canvas.draw_idle()
            return
  
        if self.annot1.get_visible() or self.annot.get_visible():
                self.annot1.set_visible(False)
                self.annot.set_visible(False)
                self.fig.canvas.draw_idle()
class ProteinInterfaceAnalyzer(QtWidgets.QWidget):
    """Main Widget for PDB Interface Analysis UI and Logic."""
    def __init__(self):
        super().__init__()
        self.setWindowTitle("PDB Interface Analyzer")
        self.setWindowFlags(QtCore.Qt.WindowStaysOnTopHint) 
        self.resize(350, 400)
        
        self.loaded_obj_name = None
        self.valid_chain_1 = None
        self.valid_chain_2 = None
        self.pdb_id_for_db = None  
        self.init_ui()

    def init_ui(self):
        """Builds the PyQt5 user interface."""
        layout = QtWidgets.QVBoxLayout()
        layout.setSpacing(10)
        layout.setContentsMargins(15, 15, 15, 15)

        # Input Group: Chains and PDB source
        group_input = QtWidgets.QGroupBox("Structure Input")
        layout_input = QtWidgets.QVBoxLayout()
        
        layout_chains = QtWidgets.QHBoxLayout()
        self.input_chain_1 = QtWidgets.QLineEdit()
        self.input_chain_1.setPlaceholderText("Chain 1 (e.g. A)")
        self.input_chain_1.setAlignment(QtCore.Qt.AlignCenter)
        
        self.input_chain_2 = QtWidgets.QLineEdit()
        self.input_chain_2.setPlaceholderText("Chain 2 (e.g. B)")
        self.input_chain_2.setAlignment(QtCore.Qt.AlignCenter)
        
        layout_chains.addWidget(QtWidgets.QLabel("Chains:"))
        layout_chains.addWidget(self.input_chain_1)
        layout_chains.addWidget(QtWidgets.QLabel("-")) 
        layout_chains.addWidget(self.input_chain_2)
        layout_input.addLayout(layout_chains)

        layout_pdb = QtWidgets.QHBoxLayout()
        self.input_pdb_id = QtWidgets.QLineEdit()
        self.input_pdb_id.setPlaceholderText("PDB ID (e.g. 1kx5)")
        self.input_pdb_id.setMaxLength(4)
        
        btn_fetch = QtWidgets.QPushButton("Fetch")
        btn_fetch.setCursor(QtCore.Qt.PointingHandCursor)
        btn_fetch.clicked.connect(self.load_by_id)
        
        layout_pdb.addWidget(self.input_pdb_id)
        layout_pdb.addWidget(btn_fetch)
        layout_input.addLayout(layout_pdb)
        
        btn_browse = QtWidgets.QPushButton("Load file")
        btn_browse.setCursor(QtCore.Qt.PointingHandCursor)
        btn_browse.clicked.connect(self.browse_file)
        layout_input.addWidget(btn_browse)
        
        group_input.setLayout(layout_input)
        layout.addWidget(group_input)

        # Parameters Group: Distance Cutoff
        group_params = QtWidgets.QGroupBox("Analysis Parameters")
        layout_params = QtWidgets.QVBoxLayout()
        
        layout_dist = QtWidgets.QHBoxLayout()
        layout_dist.addWidget(QtWidgets.QLabel("Interface Distance (Å):"))
        self.spin_dist = QtWidgets.QDoubleSpinBox()
        self.spin_dist.setValue(5.0) 
        self.spin_dist.setRange(2.0, 15.0)
        self.spin_dist.setSingleStep(0.5)
        layout_dist.addWidget(self.spin_dist)
        layout_params.addLayout(layout_dist)

        self.btn_analyze = QtWidgets.QPushButton("RUN ANALYSIS")
        self.btn_analyze.setEnabled(False)
        self.btn_analyze.setMinimumHeight(45)
        self.btn_analyze.setStyleSheet("""
            QPushButton { background-color: #3498db; color: white; font-weight: bold; border-radius: 5px; font-size: 14px; }
            QPushButton:disabled { background-color: #bdc3c7; }
            QPushButton:hover { background-color: #2980b9; }
        """)
        self.btn_analyze.clicked.connect(self.run_analysis)
        layout_params.addWidget(self.btn_analyze)

        group_params.setLayout(layout_params)
        layout.addWidget(group_params)

        # Results Group: SASA and BSA
        self.group_results = QtWidgets.QGroupBox("Surface Analysis Results")
        layout_res = QtWidgets.QVBoxLayout()
        self.res_sasa = QtWidgets.QLabel("Total SASA: -")
        self.res_bsa = QtWidgets.QLabel("Buried Surface (BSA): -")
        
        res_style = "font-weight: bold; color: #2c3e50;"
        self.res_sasa.setStyleSheet(res_style)
        self.res_bsa.setStyleSheet(res_style)
        
        layout_res.addWidget(self.res_sasa)
        layout_res.addWidget(self.res_bsa)
        self.group_results.setLayout(layout_res)
        layout.addWidget(self.group_results)

        self.status_label = QtWidgets.QLabel("Ready")
        self.status_label.setAlignment(QtCore.Qt.AlignCenter)
        self.status_label.setStyleSheet("color: gray; font-size: 10px;")
        layout.addWidget(self.status_label)

        self.setLayout(layout)

    def set_status(self, message, error=False):
        """Updates the status bar message."""
        self.status_label.setText(message)
        color = "red" if error else "#27ae60" 
        self.status_label.setStyleSheet(f"color: {color}; font-weight: bold;")
        QtWidgets.QApplication.processEvents() 

    def fetch_biological_names(self, pdb_id, c1, c2):
        """Fetches molecule names from PDBe API for given chains."""
        if not pdb_id: return ""
        try:
            url = f"https://www.ebi.ac.uk/pdbe/api/pdb/entry/molecules/{pdb_id.lower()}"
            req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
            with urllib.request.urlopen(req, timeout=3) as response:
                data = json.loads(response.read().decode())
                name_c1, name_c2 = "Unknown", "Unknown"
                for mol in data[pdb_id.lower()]:
                    if mol.get('molecule_type') == 'water':
                       continue
                    if c1 in mol.get('in_chains', []): name_c1 = mol.get('molecule_name', [''])[0]
                    if c2 in mol.get('in_chains', []): name_c2 = mol.get('molecule_name', [''])[0]
                return f"Chain {c1}: {name_c1}<br>Chain {c2}: {name_c2}"
        except:
            return ""

    def cleanup_scene(self):
        """Deletes temporary PyMOL objects from the previous analysis."""
        objects_to_delete = [
            "Analysis_Group",
            "interface_*",
            "copy_chain_*",
            "temp_complex_for_sasa",
            "Interface_Surface",
            "dist_*",
            "HydrogenBonds",
            "SaltBridges",
            "Hydrophobic",
            "Interactions",
        ]
        for obj in objects_to_delete:
            cmd.delete(obj)
        
    def check_pdb_exists(self, pdb_id):
        """Checks if PDB ID is valid by pinging RCSB servers."""
        url = f"https://files.rcsb.org/header/{pdb_id}.pdb"
        try:
            with urllib.request.urlopen(url, timeout=5) as response:
                return response.status == 200
        except: 
            return False

    def validate_and_setup(self, obj_name):
        """Validates if selected chains exist in the loaded structure."""
        c1 = self.input_chain_1.text().strip().upper()
        c2 = self.input_chain_2.text().strip().upper()
        if not c1 or not c2:
            QtWidgets.QMessageBox.warning(self, "Input Error", "Please provide both chain identifiers.")
            cmd.delete(obj_name)
            return
        available_chains = cmd.get_chains(obj_name)
        if c1 not in available_chains or c2 not in available_chains:
            QtWidgets.QMessageBox.critical(self, "Structure Error", f"Chains {c1} or {c2} not found!")
            cmd.delete(obj_name)
            self.set_status("Error: Chains not found", error=True)
            return
        self.loaded_obj_name = obj_name
        self.valid_chain_1 = c1
        self.valid_chain_2 = c2
        self.btn_analyze.setEnabled(True)
        self.set_status(f"Loaded {obj_name}. Ready to analyze {c1}-{c2}")
        cmd.orient(obj_name)

    def load_by_id(self):
        """Downloads PDB structure by ID."""
        pdb_id = self.input_pdb_id.text().strip().lower()
        if not pdb_id: return
        self.set_status(f"Fetching {pdb_id}...")
        if self.check_pdb_exists(pdb_id):
            cmd.delete("analyzer_*") 
            cmd.fetch(pdb_id)
            self.pdb_id_for_db = pdb_id 
            self.validate_and_setup(pdb_id)
        else:
            self.set_status("Connection Error", error=True)

    def browse_file(self):
        """Loads a PDB/CIF file from local storage."""
        path, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Select PDB file", filter="PDB Files (*.pdb *.cif)")
        if path:
            name = os.path.splitext(os.path.basename(path))[0]
            cmd.delete("analyzer_*")
            cmd.load(path, name)
            self.pdb_id_for_db = None 
            self.validate_and_setup(name)

    # --- Shrake-Rupley VdW radii (in Angstroms) ---
    VDW_RADII = {
        'C': 1.70, 'N': 1.55, 'O': 1.52, 'S': 1.80,
        'P': 1.80, 'H': 1.20, 'F': 1.47, 'CL': 1.75,
        'BR': 1.85, 'I': 1.98, 'SE': 1.90,
    }
    PROBE_RADIUS = 1.4   # Water molecule probe radius (Å)
    N_SPHERE_POINTS = 960  # Number of sampling points on each atom sphere

    def _generate_sphere_points(self, n):
        
        pts = np.zeros((n, 3))
        golden = np.pi * (3.0 - np.sqrt(5.0))  # Golden angle ~2.399 rad
        for i in range(n):
            y = 1.0 - (i / (n - 1.0)) * 2.0    # y goes from +1 to -1
            r = np.sqrt(max(0.0, 1.0 - y * y))  # radius at this y
            theta = golden * i                   # azimuthal angle
            pts[i] = [np.cos(theta) * r, y, np.sin(theta) * r]
        return pts

    def _shrake_rupley(self, selection):
       
        # Collect atom coordinates and elements from PyMOL
        atom_list = []
        cmd.iterate_state(
            1, selection,
            "atom_list.append((x, y, z, elem.upper()))",
            space={'atom_list': atom_list}
        )
        if not atom_list:
            return 0.0

        coords = np.array([[a[0], a[1], a[2]] for a in atom_list], dtype=float)
        elems  = [a[3] for a in atom_list]

        # Resolve VdW radius for each atom (fall back to carbon if unknown)
        radii = np.array([
            self.VDW_RADII.get(e, self.VDW_RADII['C']) for e in elems
        ], dtype=float)

        probe       = self.PROBE_RADIUS
        ext_radii   = radii + probe          # extended radius = r_atom + r_probe
        sphere_pts  = self._generate_sphere_points(self.N_SPHERE_POINTS)
        n_pts       = self.N_SPHERE_POINTS
        total_sasa  = 0.0

        for i in range(len(coords)):
            # Scale unit-sphere points to atom i's extended sphere
            test_pts = sphere_pts * ext_radii[i] + coords[i]  # shape (N, 3)

            # Find all neighbours whose extended spheres could overlap
            dists_to_i   = np.linalg.norm(coords - coords[i], axis=1)
            neighbour_mask = (
                (dists_to_i < ext_radii[i] + ext_radii) &
                (dists_to_i > 0.0)                           # exclude self
            )
            neighbour_coords  = coords[neighbour_mask]
            neighbour_ext_r   = ext_radii[neighbour_mask]

            # Count accessible test points
            if len(neighbour_coords) == 0:
                accessible = n_pts
            else:
                # Vectorised burial check:
                # diff[k, j] = distance from test_point k to neighbour j
                diff = test_pts[:, np.newaxis, :] - neighbour_coords[np.newaxis, :, :]
                dist2 = np.sum(diff ** 2, axis=2)                    # (N, n_nb)
                buried_mask = np.any(dist2 < neighbour_ext_r ** 2, axis=1)
                accessible  = int(np.sum(~buried_mask))

            # SASA contribution of atom i
            fraction    = accessible / n_pts
            atom_sasa   = 4.0 * np.pi * ext_radii[i] ** 2 * fraction
            total_sasa += atom_sasa

        return total_sasa

    def calculate_sasa_stats(self, name_1, name_2):
        """
        Calculates SASA and Buried Surface Area (BSA) using the
        Shrake-Rupley algorithm 

        BSA formula:
            BSA = SASA(chain_A_isolated) + SASA(chain_B_isolated) − SASA(complex_AB)
        """
        area_1       = self._shrake_rupley(name_1)
        area_2       = self._shrake_rupley(name_2)
        # For the complex we pass a combined selection string
        complex_sel  = f"({name_1}) or ({name_2})"
        complex_area = self._shrake_rupley(complex_sel)
        bsa          = (area_1 + area_2) - complex_area
        return complex_area, bsa

    def get_contact_matrix(self, selection1, selection2):
       
        res1 = []
        res2 = []

        cmd.iterate(
            f"({selection1}) and polymer.protein",
            "res1.append((chain, resi, resn))",
            space={"res1": res1},
        )
        cmd.iterate(
            f"({selection2}) and polymer.protein",
            "res2.append((chain, resi, resn))",
            space={"res2": res2},
        )

      
        res1 = sorted(list(set(res1)), key=lambda x: safe_int_resi(x[1]))
        res2 = sorted(list(set(res2)), key=lambda x: safe_int_resi(x[1]))

        matrix = np.full((len(res1), len(res2)), np.nan, dtype=float)
        for i, (c1, r1, _n1) in enumerate(res1):
            sel1 = f"({selection1}) and chain {c1} and resi {r1}"
            coords1 = cmd.get_coords(sel1)
            if coords1 is None:
                continue

            for j, (c2, r2, _n2) in enumerate(res2):
                sel2 = f"({selection2}) and chain {c2} and resi {r2}"
                coords2 = cmd.get_coords(sel2)
                if coords2 is None:
                    continue

                d = distance_matrix(coords1, coords2)
                matrix[i, j] = float(np.min(d))

        return matrix


    def get_residue_mapping(self, selection):
        """Maps residue numbers to matrix indices for plotting."""
        residues = []
        cmd.iterate(
            f"({selection}) and polymer.protein",
            "residues.append((chain, resi, resn))",
            space={'residues': residues}
        )
        residues = sorted(list(set(residues)), key=lambda x: safe_int_resi(x[1]))
        mapping = {}
        labels = []
        for i, (_chain, resi_raw, resn) in enumerate(residues):
         
            if resi_raw not in mapping:
                mapping[resi_raw] = i
            labels.append(f"{resn}{resi_raw}")
        return mapping, labels
   

    def categorize_interactions(self, sel_1, sel_2, cutoff):
        """Identifies specific chemical interactions: Salt Bridges, H-Bonds, Hydrophobic."""
        results = {"hbonds": [], "salt": [], "hydro": []}
        cmd.h_add(f"({sel_1}) or ({sel_2})")
        atoms_data = {}
        cmd.iterate(f"({sel_1}) or ({sel_2})", 
                    "atoms_data[(model, index)] = {'resn': resn, 'resi': resi, 'chain': chain, 'name': name, 'elem': elem}",
                    space={'atoms_data': atoms_data})

        # 1. Salt Bridges Logic
       
        anionic_filter = "(resn ASP,GLU and name OD1+OD2+OE1+OE2)"
        cationic_filter = "(resn ARG,LYS,HIS and name NZ+NH1+NH2+ND1+NE2)"
        salt_pairs = cmd.find_pairs(
            f"({sel_1}) and {cationic_filter}",
            f"({sel_2}) and {anionic_filter}",
            cutoff=4.0,
        )
        salt_pairs += cmd.find_pairs(
            f"({sel_1}) and {anionic_filter}",
            f"({sel_2}) and {cationic_filter}",
            cutoff=4.0,
        )
        seen_salt = set()
        for i, (p1, p2) in enumerate(salt_pairs):
            a1, a2 = atoms_data.get(p1), atoms_data.get(p2)
            if not a1 or not a2: continue
            dist = cmd.get_distance(f"index {p1[1]} and model {p1[0]}", f"index {p2[1]} and model {p2[0]}")
            res_key = tuple(sorted([f"{a1['resn']}{a1['resi']}", f"{a2['resn']}{a2['resi']}"]))
            if res_key not in seen_salt:
                results["salt"].append([f"{a1['resn']} {a1['resi']}", a1['chain'], f"{a2['resn']} {a2['resi']}", a2['chain'], round(dist, 2)])
                seen_salt.add(res_key)
                cmd.distance(f"dist_salt_{i}", f"index {p1[1]} and model {p1[0]}", f"index {p2[1]} and model {p2[0]}")
                cmd.color("magenta", f"dist_salt_{i}") 

       
        # 2. Hydrogen Bonds Logic (including angle check)
        cmd.h_add(f"({sel_1}) or ({sel_2})")
        donor_filter = "donor"
        acceptor_filter = "acceptor"
        hb_pairs = cmd.find_pairs(
    f"({sel_1}) and {donor_filter}",
    f"({sel_2}) and {acceptor_filter}",
    cutoff=3.5
)

        hb_pairs1 = cmd.find_pairs(
    f"({sel_1}) and {acceptor_filter}",
    f"({sel_2}) and {donor_filter}",
    cutoff=3.5
)

        hb_pairs2=[(p2,p1)for p1, p2 in hb_pairs1]   
        
        all_hb_pairs = hb_pairs + hb_pairs2
        
        for i, (p1, p2) in enumerate(all_hb_pairs):
            
            a1, a2 = atoms_data.get(p1), atoms_data.get(p2)
            if not a1 or not a2:
                continue
            
            donor_sel, acc_sel = f"index {p1[1]} and model {p1[0]}", f"index {p2[1]} and model {p2[0]}"
         
           
            hydrogens = []
            cmd.iterate(f"({donor_sel}) extend 1 and elem H",
            "hydrogens.append((model,index))",
            space={'hydrogens': hydrogens})
            
            best_angle = 0
            best_H = None
            
            for h_mod, h_idx in hydrogens:
                H_sel = f"index {h_idx} and model {h_mod}"
                HA_dist = cmd.get_distance (H_sel, acc_sel)
                
               
                if HA_dist <= 2.5:
                    current_angle = cmd.get_angle(donor_sel,H_sel, acc_sel)
                    if current_angle > best_angle:
                        best_angle = current_angle
                        best_H = H_sel
                    
            if best_angle > 140.0:
                dist = cmd.get_distance(donor_sel, acc_sel)
                
                results["hbonds"].append([f"{a1['resn']}{a1['resi']}", a1['chain'], a1['name'],
                          f"{a2['resn']}{a2['resi']}", a2['chain'], a2['name'],
                          round(dist, 2), f"{round(best_angle, 1)}"])
                cmd.distance(f"dist_hbond_{i}", donor_sel, acc_sel)
                cmd.color("yellow", f"dist_hbond_{i}") 
            
    
                
        # 3. Hydrophobic Interactions (Carbon-Carbon contacts)
        hydro_res = "resn ALA+VAL+ILE+LEU+MET+PHE+TYR+TRP"
        hydro_pairs = cmd.find_pairs(f"({sel_1}) and {hydro_res} and elem C", f"({sel_2}) and {hydro_res} and elem C", cutoff=4.0)
        seen_hydro = set()
        for i, (p1, p2) in enumerate(hydro_pairs):
            a1, a2 = atoms_data.get(p1), atoms_data.get(p2)
            if not a1 or not a2: continue
            res_key = tuple(sorted([f"{a1['resn']}{a1['resi']}", f"{a2['resn']}{a2['resi']}"]))
            if res_key not in seen_hydro:
                dist = cmd.get_distance(f"index {p1[1]} and model {p1[0]}", f"index {p2[1]} and model {p2[0]}")
                results["hydro"].append([f"{a1['resn']} {a1['resi']}", a1['chain'], f"{a2['resn']} {a2['resi']}", a2['chain'], round(dist, 2)])
                seen_hydro.add(res_key)
                cmd.distance(f"dist_hydro_{i}", f"index {p1[1]} and model {p1[0]}", f"index {p2[1]} and model {p2[0]}")
                cmd.color("gray60", f"dist_hydro_{i}")
        
   
        
        cmd.group("HydrogenBonds", "dist_hbond_*")
        cmd.group("SaltBridges", "dist_salt_*")
        cmd.group("Hydrophobic", "dist_hydro_*")
        cmd.group("Interactions", "HydrogenBonds SaltBridges Hydrophobic")
        
        
        
        return results
    def run_analysis(self):
        """Main execution flow for structure profiling."""
        self.set_status("Advanced profiling in progress...")
        self.cleanup_scene()
        c1, c2, obj, dist_cutoff = self.valid_chain_1, self.valid_chain_2, self.loaded_obj_name, self.spin_dist.value()
        name_1, name_2 = f"copy_chain_{c1}", f"copy_chain_{c2}"
        sel_1_AB = f"interface_atoms_{c1}"
        sel_2_BA = f"interface_atoms_{c2}"
        
        try:
            cmd.disable(obj)
            cmd.create(name_1, f"{obj} and chain {c1}")
            cmd.create(name_2, f"{obj} and chain {c2}")
            
            cmd.hide("everything", obj)
            cmd.show("cartoon", f"{name_1} or {name_2}")
            cmd.color("palecyan", name_1)
            cmd.color("lightpink", name_2)
            
            # Select interface residues within cutoff
            cmd.select(sel_1_AB, f"byres ({name_1} within {dist_cutoff} of {name_2})")
            cmd.select(sel_2_BA, f"byres ({name_2} within {dist_cutoff} of {name_1})")
            
            cmd.show("sticks", f"interface_atoms_*")
            cmd.util.cnc(f"interface_atoms_*")
            
            if cmd.count_atoms(f"{sel_1_AB} or {sel_2_BA}") == 0:
                QtWidgets.QMessageBox.information(self, "Result", "No contact found within specified distance.")
                self.set_status("No contacts found.", error=True)
                return

            # Interface Visualization
            for sel, c_color in [(sel_1_AB, "blue"), (sel_2_BA, "magenta")]:
                cmd.show("sticks", sel)
                cmd.color("red", f"{sel} and elem O")
                cmd.color("blue", f"{sel} and elem N")
                cmd.color("yellow", f"{sel} and elem S")
                cmd.color(c_color, f"{sel} and elem C")

            # Interface Surface Generation
            surf_name = "Interface_Surface"
            cmd.create(surf_name, f"{sel_1_AB} or {sel_2_BA}")
            cmd.show("surface", surf_name)
            cmd.hide("sticks", surf_name)
            cmd.hide("cartoon", surf_name)
            cmd.set("transparency", 0.65, surf_name)
            cmd.color("white", surf_name)

            # Calculations and Map Generation
            interactions_dict = self.categorize_interactions(name_1, name_2, dist_cutoff)
            total_sasa, bsa_val = self.calculate_sasa_stats(name_1, name_2)
            
            self.res_sasa.setText(f"Total SASA: {total_sasa:.1f} Å²")
            self.res_bsa.setText(f"Buried Surface Area (BSA): {bsa_val:.1f} Å²")
            
            # Highlight high BSA values (typically > 500 Å² for stable complexes)
            if bsa_val > 500:
                self.res_bsa.setStyleSheet("font-weight: bold; color: #27ae60;")
            else:
                self.res_bsa.setStyleSheet("font-weight: bold; color: #e67e22;")

            m1, l1 = self.get_residue_mapping(name_1)
            m2, l2 = self.get_residue_mapping(name_2)
            d_matrix = self.get_contact_matrix(name_1, name_2)
            plot_data = {'matrix': d_matrix, 'map1': m1, 'map2': m2, 'labels1': l1, 'labels2': l2}
            bio_names = self.fetch_biological_names(self.pdb_id_for_db, c1, c2)
            
            # Show secondary windows
            self.table_window = AdvancedInteractionWindow(interactions_dict, bio_names)
            self.table_window.show()
            
            if plot_data['matrix'] is not None:
                self.plot_window = InteractionPlotWindow(plot_data, interactions_dict, c1, c2)
                self.plot_window.show()
                
            self.set_status("Analysis finished!", error=False)
            cmd.zoom(f"copy_chain_*", buffer=2.0)
            
        except Exception as e:
            self.set_status(f"Error: {str(e)}", error=True)

def __init_plugin__(app=None):
    from pymol.plugins import addmenuitemqt
    
    addmenuitemqt('Analizator Interfejsu PDB', lambda: setup_gui())


dialog = None

def setup_gui():
    global dialog
    if dialog is None:
        dialog = ProteinInterfaceAnalyzer()
    dialog.show()
    dialog.raise_()
    dialog.activateWindow()