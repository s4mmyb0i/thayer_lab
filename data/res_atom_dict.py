"""
Samvit Prem Singhal | spsinghal@wesleyan.edu

Takes in a .pdb file and returns a dictionary with keys as residue numbers and values as a list of atoms in that residue
"""

import mdtraj as md

def residue_atoms_dict(pdb_file:str) -> dict[int,list[int]]:
    """
    Input: .pdb file
    Output: dict(residue # : list[atom #s])
    """
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
