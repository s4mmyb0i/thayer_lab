"""
Samvit Prem Singhal | spsinghal@wesleyan.edu

Returns the res sequence in one letter format
"""

from Bio.PDB import PDBParser, Polypeptide

# Load the PDB file
pdb_file = "your_file.pdb"
parser = PDBParser(QUIET=True)
structure = parser.get_structure("protein", pdb_file)

# Iterate over all chains
for model in structure:
    for chain in model:
        seq = ""
        for residue in chain:
            if Polypeptide.is_aa(residue, standard=True):
                seq += Polypeptide.three_to_one(residue.get_resname())
        print(f"Chain {chain.id}: {seq}")