"""
Samvit Prem Singhal | spsinghal@wesleyan.edu

Takes in a .pdb file and returns a dictionary with keys as residue numbers and values as a list of atoms in that residue
"""

import mdtraj as md

def residue_atoms_dict(pdb_file):
    # Load the structure
    traj = md.load(pdb_file)

    # Initialize the dictionary
    residue_atoms_dict = {}

    # Iterate over residues in the topology
    for residue in traj.topology.residues:
        # Create a list of atom names in the current residue
        atom_indices = [(atom.index + 1) for atom in residue.atoms]
        # Populate dictionary with residue number as key and list of atom names as values
        residue_atoms_dict[residue.resSeq] = atom_indices

    return residue_atoms_dict

# # Usage example
# pdb_file = "../../../../from_leticia/P_1BFE_RIP/1bfe.pdb"  # Replace with your PDB file path
# residue_atoms = residue_atoms_dict(pdb_file)
# print(residue_atoms)




# import pytraj as pt

# # Load the PDB file
# pdb_file = "../../../../from_leticia/P_1BFE_RIP/1bfe.pdb"
# traj = pt.load(pdb_file)

# # Create a dictionary with residue numbers as keys and atom names as values
# residue_atoms = {}

# # Loop through each residue in the trajectory
# for residue in traj.topology.residues:
#     residue_number = residue.resid  # Get residue number
#     atoms = [atom.name for atom in residue.atoms]  # Get atom names in the residue
#     residue_atoms[residue_number] = atoms

# # Print the result
# print(residue_atoms)