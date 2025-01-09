import glob
import os, sys
from optparse import OptionParser
import pandas as pd
import numpy as np
import pathlib
from collections import namedtuple
import shutil


"""
local imports
"""
from utils.pdb_to_sequence import pdb_to_sequence as p2s
from utils.bp3 import bepipred3 as bp3

"""
env imports
"""
# from immunogenicity.utils.pdb_to_sequence import pdb_to_sequence as p2s
# from immunogenicity.utils.bp3 import bepipred3 as bp3

def eval_b_cell_epitope(results_dir, window_size=10): 
    files =  glob.glob(os.path.join(results_dir, 'designability_eval', 'designs', '*.pdb'))
    save_dir = os.path.join(results_dir, 'functional_eval')
    
    for PDBFile in files:
        predict_b_cell_epitope(PDBFile, save_dir, window_size=window_size)
        


def extract_pdb_fname(pdb_path): 
    retval, _ = os.path.splitext(os.path.basename(pdb_path))
    return retval

def predict_b_cell_epitope(pdb_fpath, save_dir, window_size=10):

    design_name = extract_pdb_fname(pdb_fpath)
    aa_seq = p2s(pdb_fpath)
    esm_temp_dir = os.path.join(save_dir, f"esm_encodings_{design_name}")
    
    #generates esm embeddings in dir which gets deleted later...might be better way to do this
    antigens = bp3.Antigens(aa_seq,pathlib.Path(esm_temp_dir), design_name)
    #prediction object
    bp3_predict = bp3.BP3EnsemblePredict(antigens, window_size)
    #do prediction
    bp3_predict.run_bp3_ensemble()
    #log results
    bp3_predict.create_csvfile(pathlib.Path(save_dir), design_name)
    ##cleanup, remove esm embeddings dir
    shutil.rmtree(esm_temp_dir)

    return 


if __name__ == "__main__": 
    # test_pdb = "/home/xchen/projects/salt/results_rfdiffusion_denovo_20241018225424/folding/rf_design_0/unrelaxed_model_1_pred_0.pdb"
    
    # aa_seq = p2s(test_pdb)
    # print(aa_seq)
    # predict_b_cell_epitope(aa_seq, extract_pdb_fname(test_pdb))

    #batch test 
    dir_name = "/home/ntangella/test_data/results_genie2_motif_20241016234531/"
    eval_b_cell_epitope(dir_name)

