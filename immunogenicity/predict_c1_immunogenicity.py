import os, sys
from optparse import OptionParser
import pandas as pd

if __name__ == '__main__':
    from immunogenicity.utils.pdb_to_sequence import *
else: 
    from immunogenicity.utils import pdb_to_sequence    

import argparse
from Bio import PDB

def sliding_window_string(s, k):
    """
    for a given string representing a linear sequence of amino acids 
    this will return an array of strings with a newline space using a 
    sliding window of size k and wil throw a runtime error if  s < k
    """
    
    if(len(s) < k): 
        raise RuntimeError("window size larger than length of string")
    
    return "\n".join(s[i:i+k] for i in range(len(s) - k + 1))
    

def predict_c1_immunogenicity(pdb_fpath, c_allele="HLA-A0101", c_window_size=9, custom_mask=None): 
    """
    this wrapper around the model to predict CLass I immunogenicity (http://tools.iedb.org/immunogenicity/)
    takes in an amino acid, constructs a set of strings of length c_window_size and predicts binding to the
    c_allele choice and returns a dataframe over all peptides along with the allele and length information. 

    Specifically, this model uses features of the amino acid sequence to predict the immunogenicity of the 
    pMHC COMPLEX and not MHC I binding alone. 
    """
    aa_seq = pdb_to_sequence(pdb_fpath)
    parsed_seq = sliding_window_string(aa_seq, c_window_size).split()

    pred = Prediction()

    prediction_options = {
        "allele": str(c_allele), #predict binding to human leukocyte antigen
        "custom_mask": custom_mask
    }

    return pred.predict(pred.prep_EIT(parsed_seq,prediction_options,None))

class Prediction():
    '''
    Created on 08.10.2015
    @author: Dorjee Gyaltsen
    '''
    
    def isint(self, x):
        try:
            a = float(x)
            b = int(a)
        except ValueError: return False
        else: return a == b
    
    def prep_EIT(self, seq, options, args): 

        allele_dict = {"H-2-Db":"2,5,9","H-2-Dd":"2,3,5","H-2-Kb":"2,3,9","H-2-Kd":"2,5,9","H-2-Kk":"2,8,9","H-2-Ld":"2,5,9","HLA-A0101":"2,3,9","HLA-A0201":"1,2,9","HLA-A0202":"1,2,9","HLA-A0203":"1,2,9","HLA-A0206":"1,2,9","HLA-A0211":"1,2,9","HLA-A0301":"1,2,9","HLA-A1101":"1,2,9","HLA-A2301":"2,7,9","HLA-A2402":"2,7,9","HLA-A2601":"1,2,9","HLA-A2902":"2,7,9","HLA-A3001":"1,3,9","HLA-A3002":"2,7,9","HLA-A3101":"1,2,9","HLA-A3201":"1,2,9","HLA-A3301":"1,2,9","HLA-A6801":"1,2,9","HLA-A6802":"1,2,9","HLA-A6901":"1,2,9","HLA-B0702":"1,2,9","HLA-B0801":"2,5,9","HLA-B1501":"1,2,9","HLA-B1502":"1,2,9","HLA-B1801":"1,2,9","HLA-B2705":"2,3,9","HLA-B3501":"1,2,9","HLA-B3901":"1,2,9","HLA-B4001":"1,2,9","HLA-B4002":"1,2,9","HLA-B4402":"2,3,9","HLA-B4403":"2,3,9","HLA-B4501":"1,2,9","HLA-B4601":"1,2,9","HLA-B5101":"1,2,9","HLA-B5301":"1,2,9","HLA-B5401":"1,2,9","HLA-B5701":"1,2,9","HLA-B5801":"1,2,9"}

        for peptide in seq:
            for amino_acid in peptide.strip():
                if not amino_acid.upper() in "ACDEFGHIKLMNPQRSTVWY":
                    print("Sequence: '%s' contains an invalid character: '%c' at position %d." %(peptide, amino_acid, peptide.find(amino_acid)))
                    sys.exit(1)

        custom_mask = options["custom_mask"]
        
        if custom_mask:
            custom_mask_list = list(map(int, custom_mask.split(",")))
            if sum(n < 0 for n in custom_mask_list) > 0:
                print("custom-mask should be greater then zero.")
                sys.exit(1)

            max_length = max(custom_mask_list)
            if not all([len(sequence) >= max_length for sequence in seq]):
                print("custom length '{}' cannot be greater then the peptide length.".format(max_length))
                sys.exit(1)

        if custom_mask:
            if not all(self.isint(num) for num in custom_mask.split(",")):
                print("custom-mask should be a single number or comma-separated list of numbers.")
                sys.exit(1)
        
        allele = options["allele"]
        allele = allele.replace("*","").replace(":","") if allele else None
        
        # Check if both custom_mask and allele options given
        if custom_mask and allele:
            print("* Allele {} has default value {}.\n* When both \'custom_mask\' and \'allele\' options are used, latter takes precedence over former.\n".format(allele, allele_dict[allele], '--custom_mask'))
        
        # Check if allele is included in the available alleles
        if allele in allele_dict:
            custom_mask = allele_dict[allele]
        
        # Check if allele option is used and is in the available alleles
        if allele:
            if allele not in allele_dict:
                print("Allele {} is not available.".format(allele))
                sys.exit(1)

        return seq, custom_mask, options["allele"]
    
    def predict(self, cleaned_data):
        '''Returns the prediction result.'''
        
        immunoscale = {"A":0.127, "C":-0.175, "D":0.072, "E":0.325, "F":0.380, "G":0.110, "H":0.105, "I":0.432, "K":-0.700, "L":-0.036, "M":-0.570, "N":-0.021, "P":-0.036, "Q":-0.376, "R":0.168, "S":-0.537, "T":0.126, "V":0.134, "W":0.719, "Y":-0.012}
        immunoweight = [0.00, 0.00, 0.10, 0.31, 0.30, 0.29, 0.26, 0.18, 0.00]
        
        result_list = []
        
        
        sequence_text, custom_mask, allele = cleaned_data
        
        for pep in sequence_text:
            peptide = pep.upper()
            peplen = len(peptide)
            
            cterm = peplen - 1
            score = 0
            count = 0
            
            if not custom_mask:
                mask_num  = [0, 1, cterm]
                mask_out = [1, 2, "cterm"]
            elif custom_mask:
                try: 
                    mask_str = custom_mask.split(",")
                    mask_num = list(map(int, mask_str))
                    mask_num = list(map(lambda x: x - 1, mask_num))
                    mask_out = list(map(lambda x: x + 1, mask_num))
                except IOError as e:
                    print("I/O error({0}): {1}".format(e.errno, e.strerror))
            else:
                self.mask_num = []
                self.mask_out = [1,2, "cterm"]
                    
            if peplen > 9:
                pepweight = immunoweight[:5] + ((peplen - 9) * [0.30]) + immunoweight[5:]
            else:
                pepweight = immunoweight
                
            try:
                for pos in peptide:
                    if pos not in immunoscale.keys():
                        print(pos)
                        raise KeyError()
                    elif count not in mask_num:
                        score += pepweight[count] * immunoscale[pos]
                        count += 1
                    else:
                        count += 1
                result_list.append([peptide, len(peptide), round(score, 5)])
            
            except IOError as e:
                print("I/O error({0}): {1}".format(e.errno, e.strerror))
#                     shutil.rmtree(atemp_dir)  
#                     raise ("Error: Please make sure you are entering in correct amino acids.")
            except:
                print("Unexpected error:", sys.exc_info()[0])
                raise

        # Sort by the last column value (score)       
        result_list.sort(key=lambda tup: tup[-1], reverse=True)
        
        # Column headers for the result
        header_list= ['peptide','length','score']
        result_list.insert(0, header_list)
        
        
        # if allele:
        #     print("allele: {}".format(allele))
        # print("masking: {0}".format('custom' if custom_mask else 'default'))
        # print("masked variables: {0}\n".format(mask_out))

        result_list = [(pep, l, immuno_score, str(allele)) for (pep, l, immuno_score) in result_list[1::]] #also output allele of binder

        out_df = pd.DataFrame(result_list, 
                              columns=["peptide", "length", "score", "allele"]
        )
        return out_df

    def create_csv(self, mask_choice, mask_out, data):
        import csv 
        import tempfile
        
        tmpdir = './output'
        
        # Create a temporary file inside the tmp/ directory
        tmpfile = tempfile.NamedTemporaryFile(prefix="immunogenicity_", suffix=".csv", dir=tmpdir, delete=False)
        
        with open(tmpfile.name, 'wb') as result:
            writer = csv.writer(result, delimiter=',')
            data.insert(0, ['masking: ', '{0}'.format(mask_choice)])
            data.insert(1, ['masked variables: ', '{0}'.format(mask_out)])
            for score in data:
                writer.writerow(score)
        tmpfile.close()
        return tmpfile.name
    
    

"""
========================================================================= 
| Alleles available for Class-I Immunogenicity:                         | 
|-----------------------------------------------------------------------| 
| H-2-Db    | H-2-Dd    | H-2-Kb    | H-2-Kd    | H-2-Kk    | H-2-Ld    | 
| HLA-A0101 | HLA-A0201 | HLA-A0202 | HLA-A0203 | HLA-A0206 | HLA-A0211 | 
| HLA-A0301 | HLA-A1101 | HLA-A2301 | HLA-A2402 | HLA-A2601 | HLA-A2902 | 
| HLA-A3001 | HLA-A3002 | HLA-A3101 | HLA-A3201 | HLA-A3301 | HLA-A6801 | 
| HLA-A6802 | HLA-A6901 | HLA-B0702 | HLA-B0801 | HLA-B1501 | HLA-B1502 | 
| HLA-B1801 | HLA-B2705 | HLA-B3501 | HLA-B3901 | HLA-B4001 | HLA-B4002 | 
| HLA-B4402 | HLA-B4403 | HLA-B4501 | HLA-B4601 | HLA-B5101 | HLA-B5301 | 
| HLA-B5401 | HLA-B5701 | HLA-B5801                                     |
-------------------------------------------------------------------------
"""

"""
args    = ['example/test.txt']
options = {'allele': None, 'allele_list': False, 'custom_mask': '2,3,9'}
"""
if __name__ == '__main__':
    test_path = "/home/xchen/projects/salt/results_rfdiffusion_denovo_20241018225424/generation/rf_design_0.pdb"
    a = predict_c1_immunogenicity(test_path,"HLA-A0201",9)
    print(a.head())






















'''
original validation function that was modified for EIT workflow
    def validate(self, options, args):
        
        allele_dict = {"H-2-Db":"2,5,9","H-2-Dd":"2,3,5","H-2-Kb":"2,3,9","H-2-Kd":"2,5,9","H-2-Kk":"2,8,9","H-2-Ld":"2,5,9","HLA-A0101":"2,3,9","HLA-A0201":"1,2,9","HLA-A0202":"1,2,9","HLA-A0203":"1,2,9","HLA-A0206":"1,2,9","HLA-A0211":"1,2,9","HLA-A0301":"1,2,9","HLA-A1101":"1,2,9","HLA-A2301":"2,7,9","HLA-A2402":"2,7,9","HLA-A2601":"1,2,9","HLA-A2902":"2,7,9","HLA-A3001":"1,3,9","HLA-A3002":"2,7,9","HLA-A3101":"1,2,9","HLA-A3201":"1,2,9","HLA-A3301":"1,2,9","HLA-A6801":"1,2,9","HLA-A6802":"1,2,9","HLA-A6901":"1,2,9","HLA-B0702":"1,2,9","HLA-B0801":"2,5,9","HLA-B1501":"1,2,9","HLA-B1502":"1,2,9","HLA-B1801":"1,2,9","HLA-B2705":"2,3,9","HLA-B3501":"1,2,9","HLA-B3901":"1,2,9","HLA-B4001":"1,2,9","HLA-B4002":"1,2,9","HLA-B4402":"2,3,9","HLA-B4403":"2,3,9","HLA-B4501":"1,2,9","HLA-B4601":"1,2,9","HLA-B5101":"1,2,9","HLA-B5301":"1,2,9","HLA-B5401":"1,2,9","HLA-B5701":"1,2,9","HLA-B5801":"1,2,9"}

        sequence_text = open(args[0], "r").read().split()
        for peptide in sequence_text:
            for amino_acid in peptide.strip():
                if not amino_acid.upper() in "ACDEFGHIKLMNPQRSTVWY":
                    print("Sequence: '%s' contains an invalid character: '%c' at position %d." %(peptide, amino_acid, peptide.find(amino_acid)))
                    sys.exit(1)

        custom_mask = options.custom_mask
        
        if custom_mask:
            custom_mask_list = list(map(int, custom_mask.split(",")))
            if sum(n < 0 for n in custom_mask_list) > 0:
                print("custom-mask should be greater then zero.")
                sys.exit(1)

            max_length = max(custom_mask_list)
            if not all([len(sequence) >= max_length for sequence in sequence_text]):
                print("custom length '{}' cannot be greater then the peptide length.".format(max_length))
                sys.exit(1)

        if custom_mask:
            if not all(self.isint(num) for num in custom_mask.split(",")):
                print("custom-mask should be a single number or comma-separated list of numbers.")
                sys.exit(1)
        
        allele = options.allele
        allele = allele.replace("*","").replace(":","") if allele else None
        
        # Check if both custom_mask and allele options given
        if custom_mask and allele:
            print("* Allele {} has default value {}.\n* When both \'custom_mask\' and \'allele\' options are used, latter takes precedence over former.\n".format(allele, allele_dict[allele], '--custom_mask'))
        
        # Check if allele is included in the available alleles
        if allele in allele_dict:
            custom_mask = allele_dict[allele]
        
        # Check if allele option is used and is in the available alleles
        if allele:
            if allele not in allele_dict:
                print("Allele {} is not available.".format(allele))
                sys.exit(1)
                
        return sequence_text, custom_mask, allele
'''

