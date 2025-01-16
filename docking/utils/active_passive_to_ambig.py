#!/usr/bin/env python

"""
Python script to convert a list of active and passive residues into 
ambiguous interaction restraints for HADDOCK
"""


def active_passive_to_ambig(ambig_fname, active1, passive1, active2, passive2, segid1, segid2):
    """Convert active and passive residues to Ambiguous Interaction Restraints

    Parameters
    ----------
    active1 : list
        List of active residue numbers of the first segid

    passive1 : list
        List of passive residue numbers of the first segid

    passive2 : list
        List of passive residue numbers of the second segid

    active2 : list
        List of active residue numbers of the second segid
    
    active2 : list
        List of passive residue numbers of the second segid

    segid1 : string
        Segid to use for the first model

    segid2 : string
        Segid to use for the second model

    """
    passive1 = [int(x) for x in passive1.split()]
    passive2 = [int(x) for x in passive2.split()]
    
    all1 = active1 + passive1
    all2 = active2 + passive2
    
    f = open(ambig_fname, 'w')

    f.write('! HADDOCK AIR restraints\n')
    f.write('! HADDOCK AIR restraints for 1st selection\n')
    f.write('!\n')
    
    for resi1 in active1:
        f.write('assign (resid {:d} and segid {:s})\n'.format(resi1, segid1))
        f.write('(\n')
        c = 0
        for resi2 in all2:
            f.write('       (resid {:d} and segid {:s})\n'.format(resi2, segid2))
            c += 1
            if c != len(all2):
                f.write('        or\n')

        f.write(') 2.0 2.0 0.0\n')
    
    f.write('\n! HADDOCK AIR restraints for 2nd selection\n')
    f.write('!\n')
            
    for resi2 in active2:
        f.write('assign (resid {:d} and segid {:s})\n'.format(resi2, segid2))
        f.write('(\n')
        c = 0
        for resi1 in all1:
            f.write('       (resid {:d} and segid {:s})\n'.format(resi1, segid1))
            c += 1
            if c != len(all1):
                f.write('        or\n')

        f.write(') 2.0 2.0 0.0\n')

    f.close()

def main():
    import sys
    if len(sys.argv) != 3:
        print('\nUsage:\n     python active-passive_to_ambig.py <active-passive-file1> <active-passive-file2>\n\n' +
              'where <active-passive-file> is a file consisting of two space-delimited lines with\n' +
              'the first line active residues numbers and the second line passive residue numbers\n')
        sys.exit()

    active1, passive1 = [[int(x) for x in line.split()] for line in open(sys.argv[1])]
    active2, passive2 = [[int(x) for x in line.split()] for line in open(sys.argv[2])]
    active_passive_to_ambig(active1, passive1, active2, passive2)


if __name__ == '__main__':
    main()
