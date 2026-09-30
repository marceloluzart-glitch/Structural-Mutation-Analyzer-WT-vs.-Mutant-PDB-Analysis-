"""
Structural Mutation Analyzer
Description: Compares Wild-Type and Mutant PDB structures (AlphaFold), computes 
global/local RMSD, residue distances, and visualizes using py3Dmol.
"""

from Bio.PDB import PDBParser
from Bio.PDB import Superimposer
import numpy as np
import matplotlib.pyplot as plt
import py3Dmol
import os

# =========================================================
# FILES
# =========================================================
WT_FILE = (
    "data/raw/WT_output/"
    "WT_unrelaxed_rank_001_alphafold2_ptm_model_1_seed_000.pdb"
)

MUT_FILE = (
    "data/raw/MUT_output/"
    "C1099S_unrelaxed_rank_001_alphafold2_ptm_model_1_seed_000.pdb"
)

# =========================================================
# PARSER
# =========================================================
parser = PDBParser(QUIET=True)

wt = parser.get_structure("WT", WT_FILE)
mut = parser.get_structure("MUT", MUT_FILE)

# =========================================================
# EXTRACT C-ALPHA ATOMS
# =========================================================
wt_atoms = []
mut_atoms = []

for wt_res, mut_res in zip(wt.get_residues(), mut.get_residues()):
    if "CA" in wt_res and "CA" in mut_res:
        wt_atoms.append(wt_res["CA"])
        mut_atoms.append(mut_res["CA"])

# =========================================================
# GLOBAL RMSD
# =========================================================
sup = Superimposer()
sup.set_atoms(wt_atoms, mut_atoms)

print("\n========================")
print("GLOBAL RMSD")
print("========================")
print(f"RMSD: {sup.rms:.3f} Å")

# =========================================================
# APPLY SUPERPOSITION
# =========================================================
sup.apply(mut.get_atoms())

# =========================================================
# HOTSPOT
# =========================================================
# Amino acid 1099 (region starts at 1000) -> 0-based Python index
hotspot = 99

wt_chain = list(wt[0].get_chains())[0]
mut_chain = list(mut[0].get_chains())[0]

wt_residues = list(wt_chain.get_residues())
mut_residues = list(mut_chain.get_residues())

wt_hot = wt_residues[hotspot]
mut_hot = mut_residues[hotspot]

print("\n========================")
print("HOTSPOT")
print("========================")
print("WT:", wt_hot)
print("MUT:", mut_hot)

# =========================================================
# LOCAL DISTANCE
# =========================================================
wt_ca = wt_hot["CA"].coord
mut_ca = mut_hot["CA"].coord
distance = np.linalg.norm(wt_ca - mut_ca)

print("\n========================")
print("HOTSPOT DISTANCE")
print("========================")
print(f"{distance:.3f} Å")

# =========================================================
# LOCAL RMSD
# =========================================================
window = 10
local_wt = []
local_mut = []

for i in range(hotspot - window, hotspot + window):
    try:
        local_wt.append(wt_residues[i]["CA"])
        local_mut.append(mut_residues[i]["CA"])
    except (KeyError, IndexError):
        pass

local_sup = Superimposer()
local_sup.set_atoms(local_wt, local_mut)

print("\n========================")
print("LOCAL RMSD")
print("========================")
print(f"{local_sup.rms:.3f} Å")

# =========================================================
# RESIDUE-WISE DISTANCES
# =========================================================
distances = []
positions = []

for i, (w, m) in enumerate(zip(wt_residues, mut_residues)):
    try:
        d = np.linalg.norm(w["CA"].coord - m["CA"].coord)
        distances.append(d)
        positions.append(i)
    except (KeyError, IndexError):
        pass

# =========================================================
# PLOT & SAVE RESULTS
# =========================================================
os.makedirs("data/results", exist_ok=True)

plt.figure(figsize=(14, 6))
plt.plot(positions, distances, linewidth=2, color='royalblue', label='C-alpha Displacement')
plt.axvline(hotspot, color='red', linestyle='--', label='C1099S Hotspot')
plt.xlabel("Amino Acid Position")
plt.ylabel("Structural Distance (Å)")
plt.title("WT vs C1099S Structural Deviation Profile")
plt.legend()
plt.grid(True, linestyle=':', alpha=0.6)

plot_path = "data/results/rmsd_profile_c1099s.png"
plt.savefig(plot_path, dpi=300, bbox_inches='tight')
print(f"\n[SUCCESS] Plot saved to {plot_path}")
plt.show()
