import os
import sys
import re
import shutil
from pathlib import Path
from haddock.gear.config import load, save
from docking.utils import passive_from_active
from docking.utils import active_passive_to_ambig
from docking.utils import pdb_reres

def convert_motif_to_dict(input_str, name):
    input_ranges = input_str.split('/')
    structures = []

    total_min_length = 0
    total_max_length = 0

    # Match chain letters
    pattern = re.compile(r'([A-Z])?(\d+)-(\d+)')

    for r in input_ranges:
        match = pattern.match(r)
        if match:
            chain, start, end = match.groups()
            start = int(start)
            end = int(end)
            if chain:
                # Add motif with chain info
                structures.append({
                    "type": "motif",
                    "chain": chain,
                    "start_index": start,
                    "end_index": end,
                    "group": chain
                })
            else:
                # Add scaffold without chain info
                structures.append({
                    "type": "scaffold",
                    "min_length": start,
                    "max_length": end
                })
            # For now, default max and min to the actual total values
            total_min_length += start
            total_max_length += end

    # Build the final dictionary
    result = {
        "name": name,
        "structures": structures,
        "min_total_length": total_min_length,
        "max_total_length": total_max_length
    }

    return result

def generate_tbl(ambig_fname, molecules, chains, candidate_specs, target_specs):

    mol1 = molecules[0]
    mol2 = molecules[1]

    if candidate_specs != None:
        if len(candidate_specs) == 1: # motif dict
            mol1_dict = convert_motif_to_dict(candidate_specs[0], Path(molecules[0]).stem)
            mol1_active_list = []
            for struct in mol1_dict['structures']:
                if struct['type'] == 'motif':
                    mol1_active_list += list(range(struct['start_index'], struct['end_index'] + 1))
        else: # user specified residues 
            mol1_active_list = candidate_specs
        
        mol1_passive_list = passive_from_active.run(pdb_file = mol1, active_list=mol1_active_list, chain_id = chains[0]) 
        
    else:
        raise Exception('Specify a list of active residues for non-generated protein or both proteins')  

    if target_specs != None:
            mol2_active_list = target_specs
            mol2_passive_list = passive_from_active.run(pdb_file = mol2, active_list=mol2_active_list, chain_id = chains[1])
    else:
        raise Exception('Specify a list of active residues for non-generated protein or both proteins')
    
    ambig_fname = active_passive_to_ambig.active_passive_to_ambig(ambig_fname, mol1_active_list, mol1_passive_list, mol2_active_list, mol2_passive_list, chains[0][0], chains[1][0])
    return ambig_fname

# code adapted from pdb-tools (haddock suite)
def select_chain(out_dir, fpath, chain_set):
    residues = []
    records = ('ATOM', 'HETATM', 'TER', 'ANISOU')

    with open(fpath, 'r') as pdb_file:
        residues = pdb_file.readlines()

    molname = os.path.basename(fpath)
    with open(Path(out_dir, molname), 'w') as pdb_file:
        for line in residues:
            if line.startswith(records):
                if line[21] not in chain_set:
                    continue
                pdb_file.write(line)
    
    return Path(out_dir, molname)

def build_config(out_dir, config_file, molecules, chains, candidate_specs, target_specs):

    try:
        os.makedirs(Path(out_dir, "docking"))
    except OSError:
        pass

    out_dir = Path(out_dir, "docking")
    run_name = "run-" + Path(molecules[0]).stem + "-" + Path(molecules[1]).stem

    filepath = Path(out_dir, run_name + "-files")
    
    if os.path.exists(filepath):
        shutil.rmtree(filepath)

    if os.path.exists(Path(out_dir, run_name)):
        shutil.rmtree(Path(out_dir, run_name))

    os.makedirs(filepath)

    molecules[0] = pdb_reres.main(filepath, molecules[0], 1)
    molecules[1] = pdb_reres.main(filepath, molecules[1], 1)

    if len(chains) == 2:

        if chains[0] != ["ALL"]:
            molecules[0] = select_chain(filepath, molecules[0], chains[0])
            
        if chains[1] != ["ALL"]:
            molecules[1] = select_chain(filepath, molecules[1], chains[1])

    else: 
        raise Exception("Please provide only one list of chains to keep per molecule.")

    config_path = config_file
    config_dict = load(config_path)['loaded_cleaned_input']
    config_path = Path(filepath, os.path.basename(config_file))
    
    config_dict['run_dir'] = Path(out_dir, run_name)
    config_dict['molecules'] = molecules

    ambig_fname = Path(filepath, "ambig_fname.tbl")
    generate_tbl(ambig_fname, molecules, chains, candidate_specs, target_specs)

    config_dict['rigidbody.1']['ambig_fname'] = ambig_fname   # it0
    config_dict['flexref.1']['ambig_fname'] = ambig_fname   # it1
    config_dict['emref.1']['ambig_fname'] = ambig_fname   # itw    

    save(config_dict, config_path, False)

    return config_path