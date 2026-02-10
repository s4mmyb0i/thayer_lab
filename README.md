## **Project Overview**

This is a **molecular dynamics analysis pipeline** that studies protein-protein interactions and allosteric communication in protein complexes. The project analyzes four different protein constructs to understand how structural changes affect communication pathways within proteins.

### **The Four Protein Constructs:**
- **P**: Base protein construct (PDZ domain)
- **PL**: P + ligand binding
- **AP**: P + additional protein (Cdc42) 
- **APL**: P + additional protein + ligand (full complex)

The project uses **Markov chain analysis** and **random walk simulations** to study how information flows through these protein structures, particularly focusing on communication to a specific binding pocket (residues 50-55).

## **Workflow Breakdown**

### **1. Data Preparation** 
**Directory: `create_pdb/`**
- Converts molecular dynamics trajectory files (`.mdcrd`) and topology files (`.prmtop`) into PDB format
- Uses AMBER's `ambpdb` tool
- **Run first**: `./create_pdb/create_pdb.sh`

### **2. Dynamic Cross-Correlation Matrix (DCCM) Creation**
**Directory: `create_dccm/`**
- Analyzes molecular dynamics trajectories to compute correlation matrices between all atom pairs
- Uses AMBER's `cpptraj` tool
- Creates correlation matrices showing how atoms move together during simulation
- **Run second**: `./create_dccm/create_dccm.sh <mdcrd_dir> <prmtop_dir> <output_dir>`

### **3. Network Construction**
**Directory: `create_network/`**
- Converts DCCM matrices into weighted graphs
- Applies correlation thresholds to create sparse networks
- Builds transition matrices for Markov chain analysis
- Sets specific nodes as "absorbing" (binding pocket residues in code)
- **Run third**: `./create_network/run_create_network.sh <matrix_dir> <pdb_dir> <threshold> <graphml_dir> <adj_matrix_dir> <trans_matrix_dir> <trans_abs_matrix_dir>`

### **4. Random Walk Simulations**
**Directory: `simulate_walk/`**
- Performs random walks on the constructed networks
- Simulates how information/perturbations propagate through the protein
- Tracks paths from every atom to the binding pocket
- **Run fourth**: `./simulate_walk/run_simulate_walk.sh <matrix_dir> <max_steps> <n_walks> <output_dir>`

### **5. Absorption Probability Analysis**
**Directory: `linalg/absorption_probs/`**
- Computes theoretical absorption probabilities using Markov chain theory
- Calculates how likely each atom is to reach the binding pocket
- Provides analytical solution to complement random walk simulations
- **Run fifth**: `python3 -m linalg.absorption_probs.absorption_probs <matrix_file> <output_dir>`

### **6. Visualization and Analysis**
**Directory: `linalg/absorption_probs/`**
- Creates interactive HTML visualizations
- Generates heatmaps, bar charts, and scatter plots
- Analyzes entropy and influence patterns
- **Run sixth**: `python3 -m linalg.absorption_probs.visualize <csv_file> <output_dir>`

## **Recommended Execution Order**

```bash
# 1. Create PDB files from trajectories
cd create_pdb
./create_pdb.sh

# 2. Generate DCCM matrices
cd ../create_dccm
./create_dccm.sh ../source_files/mdcrd ../source_files/prmtop ../dccm

# 3. Build networks and transition matrices
cd ../create_network
./run_create_network.sh ../dccm ../pdb 0.8 ../matrices/graphml ../matrices/adj_matrix ../matrices/trans_matrix ../matrices/trans_abs_matrix

# 4. Run random walk simulations
cd ../simulate_walk
./run_simulate_walk.sh ../matrices/trans_abs_matrix 100 100 ../random_walk

# 5. Compute absorption probabilities
cd ../linalg/absorption_probs
python3 -m linalg.absorption_probs.absorption_probs ../matrices/trans_abs_matrix/P_trans_abs.npz .
python3 -m linalg.absorption_probs.absorption_probs ../matrices/trans_abs_matrix/PL_trans_abs.npz .
python3 -m linalg.absorption_probs.absorption_probs ../matrices/trans_abs_matrix/AP_trans_abs.npz .
python3 -m linalg.absorption_probs.absorption_probs ../matrices/trans_abs_matrix/APL_trans_abs.npz .

# 6. Generate visualizations
python3 -m linalg.absorption_probs.visualize P_trans_abs_probs.csv .
python3 -m linalg.absorption_probs.visualize PL_trans_abs_probs.csv .
python3 -m linalg.absorption_probs.visualize AP_trans_abs_probs.csv .
python3 -m linalg.absorption_probs.visualize APL_trans_abs_probs.csv .
```

## **Key Scientific Questions**

This pipeline addresses:
1. **How does ligand binding affect allosteric communication?** (P vs PL)
2. **How does protein complex formation affect communication?** (P vs AP)  
3. **What is the combined effect of both?** (P vs APL)
4. **Which atoms are most influential in transmitting information to the binding pocket?**
5. **How do communication pathways change between different protein states?**

The analysis focuses on residues 50-55 as the target binding pocket, using Markov chain theory to understand how perturbations anywhere in the protein can influence this critical region.