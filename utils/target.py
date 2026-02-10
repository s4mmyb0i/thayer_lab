"""
Given a PDB file and a set of target residue IDs, returns a mapping from atom indices to residue IDs and a set of atom indices corresponding to the target residues.
"""

from utils.walk_utils import get_res_atom_map

PDB_FILE = "/Users/sam/Documents/Wesleyan/Thayer_Lab/2025/summer/bharat/pdb/P.pdb" # all 4 constructs have the same structure, so can just use any one of them

# TARGET_RES = {32, 38, 39, 40, 41, 42, 43, 44, 45, 98, 101, 102, 105}
# TARGET_RES = {3, 4, 6, 9, 10, 11, 13, 14, 17, 19, 15, 18, 26, 78} # allosteric binding pocket
# TARGET_RES = {38, 42, 66, 102, 105} # ligand binding pocket




# TARGET_NODES = {
#     a for res in TARGET_RES
#     for a in res2atoms.get(res, [])
# }

def get_target_nodes(target_residues: set[int], pdb_file: str = PDB_FILE) -> tuple[dict[int, int], set[int]]:
    """
    Given a PDB file and a set of target residue IDs, returns a mapping from atom indices
    to residue IDs and a set of atom indices corresponding to the target residues.

    Args:
        pdb_file (str): Path to the PDB file.
        target_residues (set[int]): Set of residue IDs representing the binding pocket.

    Returns:
        tuple:
            - atom2res (dict[int, int]): Mapping from atom index to residue ID.
            - target_nodes (set[int]): Set of atom indices that belong to the target residues.
    """
    res2atoms = get_res_atom_map(pdb_file)

    atom2res = {
        atom: res for res,
        atoms in res2atoms.items() for atom in atoms
    }
    assert 0 in atom2res, "Atom 0 not found in atom2res – index mismatch?"

    target_nodes = {
        atom for res in target_residues
        for atom in res2atoms.get(res, [])
    }
    if not target_nodes:
        print("No target atoms found for the provided TARGET_RES.")

    return atom2res, target_nodes