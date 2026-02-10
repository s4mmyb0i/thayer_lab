import mdtraj as md
import numpy as np
import networkx as nx
from scipy.sparse import csr_matrix

def atom_list(pdb_file: str) -> list:
    """Generates a list of atoms from a PDB file.
    
    Args:
        pdb_file (str): Path to the PDB file
        
    Returns:
        list: List of atoms in the PDB file
    """
    pdb = md.load(pdb_file)
    atoms = pdb.topology.atoms
    return [res for res in atoms]

def res_list(pdb_file: str) -> list:
    """Generates a set of the residues from a PDB file.
    
    Args:
        pdb_file (str): Path to the PDB file
        
    Returns:
        list: List of residues in the PDB file
    """
    pdb = md.load(pdb_file)
    residues = pdb.topology.residues
    return [res for res in residues]

def find_num_atoms(pdb_file: str) -> int:
    """Finds the number of atoms in a PDB file.
    
    Args:
        pdb_file (str): Path to the PDB file
    
    Returns:
        int: Number of atoms in the PDB file
    """
    return len(atom_list(pdb_file))

def find_num_residues(pdb_file: str) -> int:
    """Finds the number of residues in a PDB file.
    
    Args:
        pdb_file (str): Path to the PDB file
    
    Returns:
        int: Number of residues in the PDB file
    """
    return len(res_list(pdb_file))

def get_reachable_transients(P: csr_matrix, absorbing: set[int]) -> np.ndarray:
    """
    Get the transient nodes in the given transition matrix that are reachable from an absorbing node.

    Args:
        P (csr_matrix): The transition matrix.
        absorbing (set[int]): The indices of the absorbing nodes.

    Returns:
        np.ndarray: The indices of the transient nodes that are reachable from an absorbing node.
    """
    G = nx.from_scipy_sparse_array(P, create_using=nx.DiGraph)
    n = P.shape[0] # type: ignore
    transient = [i for i in range(n) if i not in absorbing]
    reachable = []
    for t in transient:
        if any(r in absorbing for r in nx.descendants(G, t)):
            reachable.append(t)
    return np.array(reachable)

def check_transition_matrix(P: csr_matrix, absorbing: set[int], verbose: bool = True) -> bool:
    """
    Sanity checks for a transition matrix.
    Checks:
    - Square matrix
    - Rows sum to 1
    - Absorbing nodes have only self-loop of 1
    - Every transient state can reach at least one absorbing node
    - No dangling non-absorbing nodes

    Args:
        P: Transition matrix
        absorbing: Set of absorbing node indices
        verbose: Whether to print out the results of each test

    Returns:
        bool: Whether all tests passed
    """
    try:
        if P.shape[0] != P.shape[1]: # type: ignore
            if verbose: print("❌ Matrix is not square.")
            return False

        row_sums = P.sum(axis=1).A1
        if not np.allclose(row_sums, 1.0, atol=1e-8):
            if verbose:
                print(f"❌ Some rows are not stochastic. Min={row_sums.min():.4f}, Max={row_sums.max():.4f}")
            return False

        if not absorbing:
            if verbose: print("❌ No absorbing nodes provided.")
            return False

        for node in absorbing:
            outgoing = P[node].nonzero()[1] # type: ignore
            if len(outgoing) != 1 or outgoing[0] != node or not np.isclose(P[node, node], 1.0): # type: ignore
                if verbose: print(f"❌ Node {node} is not a proper absorbing state.")
                return False

        G = nx.from_scipy_sparse_array(P, create_using=nx.DiGraph)
        transient_nodes = [i for i in range(P.shape[0]) if i not in absorbing] # type: ignore

        # for t in transient_nodes:
        #     reachable = nx.descendants(G, t)
        #     if not any(r in absorbing for r in reachable):
        #         if verbose: print(f"❌ Transient node {t} cannot reach any absorbing state.")
        #         return False

        dead_ends = np.where(row_sums == 0)[0]
        non_absorbing_dead = [i for i in dead_ends if i not in absorbing]
        if non_absorbing_dead:
            if verbose: print(f"❌ Dangling non-absorbing nodes: {non_absorbing_dead}")
            return False

        if verbose: print("✅ Transition matrix passed all sanity checks.")
        return True

    except Exception as e:
        if verbose: print(f"❌ Exception during check: {e}")
        return False



# print("1bfe.pdb")
# print(res_list("/Users/sam/Documents/Wesleyan/Thayer_Lab/2025/summer/1bfe.pdb"))

# print("1be9.pdb")
# print(res_list("/Users/sam/Documents/Wesleyan/Thayer_Lab/2025/summer/1be9.pdb"))

# print("pdz3lc_3.pdb")
# print(res_list("/Users/sam/Documents/Wesleyan/Thayer_Lab/2025/summer/pdz3lc_3.pdb"))

# print("pdz3cdc42.pdb")
# print(res_list("/Users/sam/Documents/Wesleyan/Thayer_Lab/2025/summer/pdz3cdc42.pdb"))


# print("P")
# print(res_list("/Users/sam/Documents/Wesleyan/Thayer_Lab/2025/summer/bharat/pdb/P.pdb"))

# print("P")
# print(find_num_residues("/Users/sam/Documents/Wesleyan/Thayer_Lab/2025/summer/bharat/pdb/P.pdb"))

# print("1bfe.pdb")
# print(res_list("/Users/sam/Documents/Wesleyan/Thayer_Lab/2025/summer/1bfe.pdb"))

# print("PL")
# print(res_list("/Users/sam/Documents/Wesleyan/Thayer_Lab/2025/summer/bharat/pdb/PL.pdb"))

# print("AP")
# print(res_list("/Users/sam/Documents/Wesleyan/Thayer_Lab/2025/summer/bharat/pdb/AP.pdb"))

# print("APL")
# print(res_list("/Users/sam/Documents/Wesleyan/Thayer_Lab/2025/summer/bharat/pdb/APL.pdb"))


# P = "[ILE1, VAL2, ILE3, SER4, MET5, PRO6, GLN7, ASP8, PHE9, ARG10, PRO11, VAL12, SER13, SER14, ILE15, ILE16, ASP17, VAL18, ASP19, ILE20, LEU21, PRO22, GLU23, THR24, HIS25, ARG26, ARG27, VAL28, ARG29, LEU30, CYS31, LYS32, TYR33, GLY34, THR35, GLU36, LYS37, PRO38, LEU39, GLY40, PHE41, TYR42, ILE43, ARG44, ASP45, GLY46, SER47, SER48, VAL49, ARG50, VAL51, THR52, PRO53, HIS54, GLY55, LEU56, GLU57, LYS58, VAL59, PRO60, GLY61, ILE62, PHE63, ILE64, SER65, ARG66, LEU67, VAL68, PRO69, GLY70, GLY71, LEU72, ALA73, GLN74, SER75, THR76, GLY77, LEU78, LEU79, ALA80, VAL81, ASN82, ASP83, GLU84, VAL85, LEU86, GLU87, VAL88, ASN89, GLY90, ILE91, GLU92, VAL93, SER94, GLY95, LYS96, SER97, LEU98, ASP99, GLN100, VAL101, THR102, ASP103, MET104, MET105, ILE106, ALA107, ASN108, SER109, ARG110, ASN111, LEU112, ILE113, ILE114, THR115, VAL116, ARG117, PRO118, ALA119, ASN120, GLN121, ARG122, ASN123]"

# PL = "[ILE1, VAL2, ILE3, SER4, MET5, PRO6, GLN7, ASP8, PHE9, ARG10, PRO11, VAL12, SER13, SER14, ILE15, ILE16, ASP17, VAL18, ASP19, ILE20, LEU21, PRO22, GLU23, THR24, HIS25, ARG26, ARG27, VAL28, ARG29, LEU30, CYS31, LYS32, TYR33, GLY34, THR35, GLU36, LYS37, PRO38, LEU39, GLY40, PHE41, TYR42, ILE43, ARG44, ASP45, GLY46, SER47, SER48, VAL49, ARG50, VAL51, THR52, PRO53, HIS54, GLY55, LEU56, GLU57, LYS58, VAL59, PRO60, GLY61, ILE62, PHE63, ILE64, SER65, ARG66, LEU67, VAL68, PRO69, GLY70, GLY71, LEU72, ALA73, GLN74, SER75, THR76, GLY77, LEU78, LEU79, ALA80, VAL81, ASN82, ASP83, GLU84, VAL85, LEU86, GLU87, VAL88, ASN89, GLY90, ILE91, GLU92, VAL93, SER94, GLY95, LYS96, SER97, LEU98, ASP99, GLN100, VAL101, THR102, ASP103, MET104, MET105, ILE106, ALA107, ASN108, SER109, ARG110, ASN111, LEU112, ILE113, ILE114, THR115, VAL116, ARG117, PRO118, ALA119, ASN120, GLN121, ARG122, ASN123]"

# AP = "[ILE1, VAL2, ILE3, SER4, MET5, PRO6, GLN7, ASP8, PHE9, ARG10, PRO11, VAL12, SER13, SER14, ILE15, ILE16, ASP17, VAL18, ASP19, ILE20, LEU21, PRO22, GLU23, THR24, HIS25, ARG26, ARG27, VAL28, ARG29, LEU30, CYS31, LYS32, TYR33, GLY34, THR35, GLU36, LYS37, PRO38, LEU39, GLY40, PHE41, TYR42, ILE43, ARG44, ASP45, GLY46, SER47, SER48, VAL49, ARG50, VAL51, THR52, PRO53, HIS54, GLY55, LEU56, GLU57, LYS58, VAL59, PRO60, GLY61, ILE62, PHE63, ILE64, SER65, ARG66, LEU67, VAL68, PRO69, GLY70, GLY71, LEU72, ALA73, GLN74, SER75, THR76, GLY77, LEU78, LEU79, ALA80, VAL81, ASN82, ASP83, GLU84, VAL85, LEU86, GLU87, VAL88, ASN89, GLY90, ILE91, GLU92, VAL93, SER94, GLY95, LYS96, SER97, LEU98, ASP99, GLN100, VAL101, THR102, ASP103, MET104, MET105, ILE106, ALA107, ASN108, SER109, ARG110, ASN111, LEU112, ILE113, ILE114, THR115, VAL116, ARG117, PRO118, ALA119, ASN120, GLN121, ARG122, ASN123]"

# APL = "[ILE1, VAL2, ILE3, SER4, MET5, PRO6, GLN7, ASP8, PHE9, ARG10, PRO11, VAL12, SER13, SER14, ILE15, ILE16, ASP17, VAL18, ASP19, ILE20, LEU21, PRO22, GLU23, THR24, HIS25, ARG26, ARG27, VAL28, ARG29, LEU30, CYS31, LYS32, TYR33, GLY34, THR35, GLU36, LYS37, PRO38, LEU39, GLY40, PHE41, TYR42, ILE43, ARG44, ASP45, GLY46, SER47, SER48, VAL49, ARG50, VAL51, THR52, PRO53, HIS54, GLY55, LEU56, GLU57, LYS58, VAL59, PRO60, GLY61, ILE62, PHE63, ILE64, SER65, ARG66, LEU67, VAL68, PRO69, GLY70, GLY71, LEU72, ALA73, GLN74, SER75, THR76, GLY77, LEU78, LEU79, ALA80, VAL81, ASN82, ASP83, GLU84, VAL85, LEU86, GLU87, VAL88, ASN89, GLY90, ILE91, GLU92, VAL93, SER94, GLY95, LYS96, SER97, LEU98, ASP99, GLN100, VAL101, THR102, ASP103, MET104, MET105, ILE106, ALA107, ASN108, SER109, ARG110, ASN111, LEU112, ILE113, ILE114, THR115, VAL116, ARG117, PRO118, ALA119, ASN120, GLN121, ARG122, ASN123]"

# print(P == PL)
# print(P == AP)
# print(P == APL)