import os
import sys
import re
from pathlib import Path
from haddock.gear.config import load, save
from docking.utils import passive_from_active
from docking.utils import active_passive_to_ambig

## REPLACE WITH KENNETH'S CODE
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

def generate_tbl(ambig_fname, molecules, motif_specs, residue_specs):

    mol1 = molecules[0]
    mol2 = molecules[1]

    if residue_specs != None: # determining active / passive residues using motif specs
        if len(residue_specs) == 1:
            mol2_active_list = residue_specs[0]
            mol2_passive_list = passive_from_active.run(pdb_file = mol2, active_list=mol2_active_list)
            
        elif len(residue_specs) == 2:
            mol1_active_list = residue_specs[0]
            mol2_active_list = residue_specs[1]

            mol1_passive_list = passive_from_active.run(pdb_file = mol1, active_list=mol1_active_list)
            mol2_passive_list = passive_from_active.run(pdb_file = mol2, active_list=mol2_active_list)
        else:
            raise Exception('Specify a list of active residues for non-generated protein or both proteins')
            

    if motif_specs != None: # user gives list with active / passive residues

        mol1_dict = convert_motif_to_dict(motif_specs, Path(molecules[0]).stem)
        mol1_active_list = []

        for struct in mol1_dict['structures']:
            if struct['type'] == 'motif':
                mol1_active_list += list(range(struct['start_index'], struct['end_index'] + 1))
        mol1_passive_list = passive_from_active.run(pdb_file = mol1, active_list=mol1_active_list)
                
    ambig_fname = active_passive_to_ambig.active_passive_to_ambig(ambig_fname, mol1_active_list, mol1_passive_list, mol2_active_list, mol2_passive_list)
    return ambig_fname

# code adapted from pdb-tools (haddock suite)
def select_chain(out_dir, fpath, chain_set):
    residues = []
    records = ('ATOM', 'HETATM', 'TER', 'ANISOU')

    with open(fpath, 'r') as pdb_file:
        residues = pdb_file.readlines()

    file = os.path.basename(fpath)
    with open(Path(out_dir, "docking", file), 'w') as pdb_file:
        for line in residues:
            if line.startswith(records):
                if line[21] in chain_set:
                    pdb_file.write(line)


def build_config(out_dir, config_file, molecules, chains, motif_specs, residue_specs):

    if len(chains) == 2:
        if chains[0] != ["ALL"]:
            select_chain(out_dir, molecules[0], chains[0])
        elif chains[1] != ["ALL"]:
            select_chain(out_dir, molecules[1], chains[1])
    else: 
        raise Exception("Please provide only one list of chains to keep per molecule.")

    config_path = "docking_inputs/config.cfg"
    config_dict = load(config_path)['loaded_cleaned_input']
    config_path = Path(out_dir, "docking", config_file)

    run_name = "run-" + Path(molecules[0]).stem + "-" + Path(molecules[1]).stem
    config_dict['run_dir'] = Path(out_dir, "docking", run_name)
    config_dict['molecules'] = molecules

    ambig_fname = Path(out_dir, "docking", "ambig_fname")
    generate_tbl(ambig_fname, molecules, motif_specs, residue_specs)

    config_dict['rigidbody.1']['ambig_fname'] = ambig_fname   # it0
    config_dict['flexref.1']['ambig_fname'] = ambig_fname   # it1
    config_dict['emref.1']['ambig_fname'] = ambig_fname   # itw    

    save(config_dict, config_path, False)

    return config_path