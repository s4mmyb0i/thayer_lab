"""
Utilities for parsing paths, working with coordinates, and residue mappings.
"""

import pandas as pd
from typing import Dict, List, Tuple, Set
from walk_utils import get_res_atom_map
from constants import START_RES, END_RES


def parse_path(path_str: str) -> List[int]:
    """
    Parse a comma-separated path string into a list of atom IDs.
    
    Args:
        path_str (str): Comma-separated string of atom IDs
        
    Returns:
        List[int]: List of atom IDs as integers
    """
    if pd.isna(path_str) or path_str == '':
        return []
    return [int(x.strip()) for x in str(path_str).split(',') if x.strip()]


def create_atom_coord_map(pdb_file: str) -> Dict[int, Tuple[float, float, float]]:
    """
    Create a mapping from atom_id to (x, y, z) coordinates.
    Handles both 0-based and 1-based atom IDs.
    
    Args:
        pdb_file (str): Path to the PDB file
        
    Returns:
        Dict[int, Tuple[float, float, float]]: Mapping from atom_id to coordinates
    """
    df = extract_xyz_dataframe(pdb_file)
    coord_map = {}
    # Use vectorized operations instead of iterrows for better performance
    atom_ids_1based = df['atom_id'].astype(int).values
    atom_ids_0based = atom_ids_1based - 1
    coords = list(zip(df['x'].values, df['y'].values, df['z'].values))
    
    # Map both 0-based and 1-based to same coordinates
    for atom_id_1based, atom_id_0based, coord in zip(atom_ids_1based, atom_ids_0based, coords):
        coord_map[atom_id_1based] = coord
        coord_map[atom_id_0based] = coord
    return coord_map


def get_residue_atom_sets(pdb_file: str) -> Tuple[Set[int], Set[int]]:
    """
    Get sets of atom IDs for START_RES and END_RES residues.
    
    Args:
        pdb_file (str): Path to the PDB file
        
    Returns:
        Tuple[Set[int], Set[int]]: Tuple of (start_atoms, end_atoms) sets
    """
    res_atom_map = get_res_atom_map(pdb_file)
    
    start_atoms: Set[int] = set()
    end_atoms: Set[int] = set()
    
    for res_id in START_RES:
        if res_id in res_atom_map:
            # res_atom_map returns 0-based indices, but we need to handle both
            for atom_idx in res_atom_map[res_id]:
                start_atoms.add(atom_idx)
                start_atoms.add(atom_idx + 1)  # Also add 1-based for compatibility
    
    for res_id in END_RES:
        if res_id in res_atom_map:
            for atom_idx in res_atom_map[res_id]:
                end_atoms.add(atom_idx)
                end_atoms.add(atom_idx + 1)  # Also add 1-based for compatibility
    
    return start_atoms, end_atoms

import mdtraj as md
from typing import List, Tuple
import pandas as pd

def extract_xyz(pdb_file: str) -> List[Tuple[int, int, float, float, float]]:
    """
    Extract xyz coordinates from a PDB file.
    
    Args:
        pdb_file (str): Path to the PDB file
        
    Returns:
        List[Tuple[int, int, float, float, float]]: List of tuples containing
            (atom_id, res_id, x, y, z) for each atom
    """
    # Load the PDB file
    traj = md.load(pdb_file)
    
    # Extract coordinates and metadata
    coordinates = []
    
    # Get the first frame (PDB files typically have one frame)
    positions = traj.xyz[0]
    
    # Iterate over all atoms
    for atom in traj.topology.atoms:
        atom_id = atom.index  # 0-based index
        res_id = atom.residue.resSeq  # Residue sequence number
        x, y, z = positions[atom_id]
        coordinates.append((1 + atom_id, res_id, 10 * x, 10 * y, 10 * z))
        # coordinates.append((atom_id, res_id, float(x), float(y), float(z)))
    
    return coordinates


def extract_xyz_dataframe(pdb_file: str):
    """
    Extract xyz coordinates from a PDB file and return as a pandas DataFrame.
    
    Args:
        pdb_file (str): Path to the PDB file
        
    Returns:
        pd.DataFrame: DataFrame with columns ['atom_id', 'res_id', 'x', 'y', 'z']
    """
    
    coordinates = extract_xyz(pdb_file)
    # Create DataFrame from list of tuples
    df = pd.DataFrame(
        coordinates,
        columns=['atom_id', 'res_id', 'x', 'y', 'z']  # type: ignore
    )
    
    return df