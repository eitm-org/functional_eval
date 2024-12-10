#!/usr/bin/env python

import os, sys
from Bio import PDB, SeqIO


def save_fasta(fasta_str, filename="test_sequence.fasta"):
    
    
    # Open the file in write mode and save the string
    with open(filename, 'w') as file:
        file.write(fasta_str)



def get_pdb_name(fpath): 
    return os.path.splitext(os.path.basename(fpath))[0]

def pdb_to_fasta(pdb_file, *args): 

    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)

    fasta_string = ""
    
    # Iterate through the structure to extract chains and sequences
    for model in structure:
        for chain in model:
            # Extract the sequence from the chain, exclude water molecules if generated
            seq = ','.join([res.get_resname() for res in chain if res.get_resname() != 'HOH'])            
            seq_1letter = ''.join(PDB.Polypeptide.index_to_one(PDB.Polypeptide.three_to_index(c)) for c in seq.split(","))
            
            fasta_string += f">{get_pdb_name(pdb_file)}_{chain.get_id()}\n{seq_1letter}\n"

    return fasta_string
            






if __name__ == "__main__": 
    pdb_file = "/home/ntangella/test_data/results_genie2_motif_20241016234531/generation/8hgo_0.pdb"
    fasta_seq = pdb_to_fasta(pdb_file)   
    save_fasta(fasta_seq)
    

