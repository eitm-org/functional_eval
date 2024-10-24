import toml, tomlkit
from pathlib import Path
from haddock3.src.haddock.gear.config import load, save
from utils import passive_from_active
from utils import active_passive_to_ambig

def convert_motif_to_dict(motif_specs):
    return {"name": "8HGO", 
            "structures": [{"type": "scaffold", "min_length": 0, "max_length": 72}, 
                           {"type": "motif", "chain": "C", "start_index": 24, "end_index": 51, "group": "A"}, 
                           {"type": "scaffold", "min_length": 0, "max_length": 72}], 
            "min_total_length": 100, "max_total_length": 100} 

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
        mol1_dict = convert_motif_to_dict(motif_specs)
        mol1_active_list = []

        for struct in mol1_dict['structures']:
            if struct['type'] == 'motif':
                mol1_active_list += list(range(struct['start_index'], struct['end_index'] + 1))
        mol1_passive_list = passive_from_active.run(pdb_file = mol1, active_list=mol1_active_list)
                
    ambig_fname = active_passive_to_ambig.active_passive_to_ambig(ambig_fname, mol1_active_list, mol1_passive_list, mol2_active_list, mol2_passive_list)
    return ambig_fname

def build_config(out_dir, config_file, ambig_fname, molecules, motif_specs, residue_specs):

    inp = """\
    run_dir: 
    mode: local
    ncores: 40
    molecules: []

    topoaa:
        autohis: true

    rigidbody:
        ambig_fname: 
        sampling: 20
        tolerance: 20

    caprieval.1: {}

    flexref:
        ambig_fname: 
        tolerance: 20

    caprieval.2: {}

    emref:
        ambig_fname: 
        tolerance: 20

    clustfcc:
        min_population: 1

    caprieval.3: {}
    """

    # yaml = YAML()
    # config = yaml.load(inp)

    config_path = Path(config_file)
    config_dict = load(config_path)['loaded_cleaned_input']
    
    config_dict['run_dir'] = out_dir
    config_dict['molecules'] = molecules

    generate_tbl(ambig_fname, molecules, motif_specs, residue_specs)

    config_dict['rigidbody.1']['ambig_fname'] = ambig_fname   # it0
    config_dict['flexref.1']['ambig_fname'] = ambig_fname   # it1
    config_dict['emref.1']['ambig_fname'] = ambig_fname   # itw    

    save(config_dict, config_path, False)

    return config_file