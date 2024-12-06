import glob
import os, sys
from optparse import OptionParser
import pandas as pd
from collections import namedtuple

"""
local imports
"""
from utils.pdb_to_fasta import pdb_to_fasta as p2f
from utils.pdb_to_sequence import pdb_to_sequence as p2s
from utils.mhcii_netmhciipan_4_2.mhcii_netmhciipan_4_2_el_percentile_data import percentile_manager as mhcii_netmhciipan_42_el_percentile_manager
from utils.netmhciipan_4_2_executable.netmhciipan_4_2_executable import single_prediction as single_prediction_netmhciipan42

"""
imports via env
"""
# from immunogenicity.utils import pdb_to_fasta as p2f
# from immunogenicity.utils import pdb_to_sequence as p2s
# from immunogenicity.utils.mhcii_netmhciipan_4_2.mhcii_netmhciipan_4_2_el_percentile_data import percentile_manager as mhcii_netmhciipan_42_el_percentile_manager
# from immunogenicity.utils.netmhciipan_4_2_executable.netmhciipan_4_2_executable import single_prediction as single_prediction_netmhciipan42



def eval_immunogenicity(results_dir, window_size=15, allele_name="DRB1*01:01"): 
    files =  glob.glob(os.path.join(results_dir, 'designability_eval', 'designs', '*.pdb'))
    
    for PDBFile in files:
        design_name = os.path.splitext(os.path.basename(PDBFile))[0]
        mhc2_pred = predict_mhc_2_binding(PDBFile, window_size, allele_name)
        mhc2_pred.to_csv(os.path.join(results_dir, f"functional_eval/mhc2_binding_{design_name}.csv"))
                                                 
def do_netmhciipan_42_el_prediction(sequence_list, allele_length_pairs, coreseq_len=9):

        
        allele_length_pairs = tuple(allele_length_pairs)

        _, binding_length = allele_length_pairs
            
        if binding_length > 30 or binding_length < 11:
            msg = 'Only a binding_length between 11-30 is supported, found {}' \
                .format(binding_length)
            raise Exception(msg)
        if coreseq_len != 9:
            msg = 'Only a core sequence length of 9 is supported, found {}' \
                .format(coreseq_len)
            raise Exception(msg)

        """
        returns 2d array with information in the col_names variable
        """
        single_prediction_result = single_prediction_netmhciipan42(sequence_list, allele_length_pairs, el=True)
        

        mhcii_col_names = [
            "pos", "mhc_allele", "peptide", "of", "core", "core_rel", 
            "identity", "score_el", "prcentile_rank_el", "exp_bin", "bind_level"
        ]

        df = pd.DataFrame(single_prediction_result, columns=mhcii_col_names)
        df.set_index('pos', inplace=True)
        
        return df
   
def predict_mhc_2_binding(pdb_file, window_size=15, allele_name="DRB1*01:01"): 
    
    seq = p2s(pdb_file) 
    parsed_allele=allele_name
    pair = tuple((parsed_allele, window_size))
    return do_netmhciipan_42_el_prediction(seq, pair)

if __name__ == "__main__":  
    pdb_file = ""
    df = predict_mhc_2_binding(pdb_file)
    print(df.head())

    


    
    