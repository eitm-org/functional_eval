''' 
Haddock3 Protein-Protein Docking
Example: E2A (PDB ID 1F3G) + HPr (PDB ID 1HDN)
         Native complex (PDB ID 1GGR)

Software needs: PyMOL, HADDOCK2.4

it0: randomization of orientations & rigid-body minimization
- interacting partners treated as rigid bodies (all geometrical parameters are frozen)
- partners separated in space and rotated randomly
- rigid body energy minimization step: allowed to rotate + optimize interaction

it1: semi-flexible simulated annealing in torsion angle space
- introduces flexibility to interacting partners
    - interface region, then backbone and side-chains

user can decide: refinement in cartesian space w explicit solvent (water)
- supports water (TIP3P) and DMSO environments (membrane mimic)

'''

'''
pseudocode steps

1. use pdb-tools to clean pdbs
2. remove water molecules (maybe this happens in step 1)
3. add these to the list to be set as molecules in config
4. run haddock-tools reduce to extract protonation state (reduce must be in path)
5. take output to define alternative protonation states for histidines
6. this needs to be set for topoaa (maybe just set autohis = true?????)

7. rigidbody - sampling module
    - rigid body energy minimization with CNS
    - caprieval (this calculates metrics - good for checkpointing)
    - tolerance + ambig_fname + sampling
8. flexref - model refinement module
    - semi-flexible refinement 
    - tolerance + ambig_fname
    - caprieval
9. emref - refinement by energy minimization
    - tolerance + ambig_fname
10. mdref 0 
11. clustfcc -> min_population = 1
'''



## preparation steps
# pdb-tools -> pdb_validate.py 
# remove any irrelevant water and other small molecules
#   do leave relevant co-factors
# HADDOCK server can set passive residues for us

### configuration 
# run_dir (output directory)
# execution mode -> mode + ncores
# molecules = list of molecules to be docked

## parameters
# default HADDOCK treats all histidines as doubly protonated + positively charged
# if structure contains Histidines, need to check what protonation state should be
    # use PROPKA or MolProbity (Haddock webserver)
    # haddock-tools contains script that runs reduce

## topoaa : all atom topology
# autohis
