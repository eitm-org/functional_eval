import pandas as pd
import sys, os
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio import SeqIO

def make_row(dict, dir, PDBFile, pH = []):
    sequence = ""

    with open(dir + 'folding/' + PDBFile, 'r') as pdb_file:
        for record in SeqIO.parse(pdb_file, 'pdb-atom'):
            sequence += str(record.seq)

    analysis = ProteinAnalysis(sequence)

    params = {'count_amino_acids_'+a: v for a, v in analysis.count_amino_acids().items()}
    params.update({'amino_acids_percent_'+a: v for a, v in analysis.get_amino_acids_percent().items()})
    params['molecular_weight'] = "%0.2f" % analysis.molecular_weight()
    params['aromaticity'] = "%0.2f" % analysis.aromaticity()
    params['instability_index'] = "%0.2f" % analysis.instability_index()
    params['flexibility'] = analysis.flexibility()
    params['gravy'] = analysis.gravy(scale="KyteDoolitle")
    params['isoelectric_point'] = "%0.2f" % analysis.isoelectric_point()
    params['extinction_coeff'] = analysis.molar_extinction_coefficient()[0]
    params['extinction_reduced'] = analysis.molar_extinction_coefficient()[1]

    if len(pH) > 0:
        for p in pH: 
            params['charge_at_pH_'+str(p)] = analysis.charge_at_pH(p)

    dict[PDBFile] = params
    return dict

def eval(dir, pH_list):

    files =  os.listdir(dir + "folding/")
    props_dict = {}
    for PDBFile in files:
        props_dict = make_row(props_dict, dir, PDBFile, pH_list)

    pd.DataFrame.from_dict(props_dict, orient = 'index').to_csv(dir + "functional_eval/protparam.csv")