"""
Samvit Prem Singhal | spsinghal@wesleyan.edu

Takes in a .pdb file and returns a dictionary with keys as residue numbers and values as a list of atoms in that residue
"""

import mdtraj as md

def get_res_atom_map(pdb_file: str) -> dict[int, list[int]]:
    """Returns a mapping of residue numbers to the atom numbers from a PDB file.
    
    Args:
        pdb_file (str): Path to the PDB file
        
    Returns:
        dict (int, list[int]): A dictionary with residue number as the keys and a list of the atom numbers as values
    """
    # Load the structure
    traj = md.load(pdb_file)

    # Initialize the dictionary
    mapping = {}

    # Iterate over residues in the topology
    for residue in traj.topology.residues:
        # Create a list of atom names in the current residue
        atom_indices = [atom.index for atom in residue.atoms]
        # Populate dictionary with residue number as key and list of atom names as values
        mapping[residue.resSeq] = atom_indices

    return mapping


# pdb = "/Users/sam/Documents/Wesleyan/Thayer_Lab/2025/summer/bharat/pdb/P.pdb"

# print(len(get_res_atom_map(pdb)))

# for res_num, atom_list in get_res_atom_map(pdb).items():
#     print(f"Residue {res_num}: Atoms {atom_list}")


"""
Do we need scaling factor on A side since it's bigger than L side?
Estimate the surface area of both surfaces and get scalar value from it
Scalar value can be applied post priori. When making any comparison
"""