'''
Calculate the RMSD between the crystal structures of two CD98hc proteins or between loop regions within these structures.
'''
import Bio.PDB
from Bio.PDB import PDBParser
from Bio.PDB import Vector

pdb1 = '2dh2-hcd98hc.pdb'
AF3_1 = 'AF3-HCD98HC.pdb'
pdb2 = '6i9q-mcd98hc.pdb'
AF3_2 = 'AF3-MCD98HC.pdb'
region1 = [263, 275, 280, 291]
region2 = [263, 275, 280, 291]


region3 = [264, 275, 280, 291]
region4 = [264, 275, 280, 291]

def RMSD(structure_1, structure_2):
	
    parser = PDBParser()
    STR1 = parser.get_structure('0', structure_1)[0]['A']
    STR2 = parser.get_structure('1', structure_2)[0]['A']

    STR1 = [res for res in STR1]
    #print(fixed)
    STR2 = [res for res in STR2]

	
    #fixed  = [atom for aa in STR1 for atom in aa]
    #moving = [atom for aa in STR2 for atom in aa]

    fixed  = [atom for aa in STR1 for atom in aa if atom.get_id() in ['CA', 'C', 'N']]
    #print(fixed)
    moving = [atom for aa in STR2 for atom in aa if atom.get_id() in ['CA', 'C', 'N']]
    #print(moving)
	#print(fixed)
	#print(moving)
    lengths = [len(fixed), len(moving)]
    print(lengths)
    smallest = min(lengths)
    sup = Bio.PDB.Superimposer()
    sup.set_atoms(fixed[:smallest], moving[:smallest])
    sup.apply(moving)
    RMSD = round(sup.rms, 4)
    print(RMSD)

def RMSD_loop(structure_1, structure_2, region1, region2):
	
    parser = PDBParser()
    STR1 = parser.get_structure('0', structure_1)[0]['A']
    STR2 = parser.get_structure('1', structure_2)[0]['A']

    STR1 = [res for res in STR1]
    #print(fixed)
    STR2 = [res for res in STR2]

    STR1 = STR1[region1[0]: region1[1]] + STR1[region1[2]: region1[3]]
    print(STR1)

    STR2 = STR2[region2[0]: region2[1]] + STR2[region2[2]: region2[3]]
    print(STR2)
	
    #fixed  = [atom for aa in STR1 for atom in aa]
    #moving = [atom for aa in STR2 for atom in aa]

    fixed  = [atom for aa in STR1 for atom in aa if atom.get_id() in ['CA', 'C', 'N']]
    #print(fixed)
    moving = [atom for aa in STR2 for atom in aa if atom.get_id() in ['CA', 'C', 'N']]
    #print(moving)
	#print(fixed)
	#print(moving)
    lengths = [len(fixed), len(moving)]
    print(lengths)
    smallest = min(lengths)
    sup = Bio.PDB.Superimposer()
    sup.set_atoms(fixed[:smallest], moving[:smallest])
    sup.apply(moving)
    RMSD = round(sup.rms, 4)
    print(RMSD)

RMSD(pdb1, AF3_1)
RMSD(pdb2, AF3_2)

RMSD_loop(pdb1, AF3_1, region1, region2)
RMSD_loop(pdb2, AF3_2, region3, region4)
