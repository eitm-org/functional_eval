from Bio import PDB

def pdb_to_sequence(pdb_file, *args):
    """
    take in a pdb file and output a string of one letter AA codes
    irrespective of chain for use in immunogenicity predictions
    """
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("seq_to_struct", pdb_file)
    amino_acids = []

    # Iterate over the residues in the structure
    for model in structure:
        for chain in model:
            for residue in chain:
                if PDB.is_aa(residue):
                    amino_acids.append(residue.get_resname())

    #there is, in fact a three to one conversion, but is in a release that is not on conda-forge
    one_letter_codes = [PDB.Polypeptide.index_to_one(PDB.Polypeptide.three_to_index(res)) for res in amino_acids]

    return "".join(one_letter_codes)

if __name__ == "__main__": 
    # Example usage
    pdb_file = ""  
    amino_acid_sequence = pdb_to_sequence(pdb_file)
    print(amino_acid_sequence)