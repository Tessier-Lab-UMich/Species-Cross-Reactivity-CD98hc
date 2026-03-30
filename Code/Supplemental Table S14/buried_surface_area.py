'''
This code is used to identify amino acids within a specified distance range in an antibody–antigen complex structure (defined as epitope and paratope) and to calculate the buried surface area resulting from antibody–antigen binding.
'''

from afpdb.afpdb import Protein,RS,RL,ATS
import numpy as np
import pandas as pd




def get_epitope_paratope(pdb_file, antibody_chains, antigen_chain, dis_threshold, file_path):
    # load the ab-ag complex structure using PDB code name
    #p=Protein("test5cil.pdb")

    p=Protein(file_path + pdb_file + '.pdb')
    #p=Protein(pdb_file)

    #print(file_path + pdb_file + '.pdb')
    # show key statistics summary of the structure
    #p.summary()
    #print("Old P chain residue numbering:", p.rs("P").name(), "\n")
    #q, old_num=p.renumber("RESTART")
    #rint("New P chain residue numbering:", q.rs("P").name(), "\n")
    #q.summary()

    # We can further remove the insertion code by
    # renumber() returns a tuple, (Protein obj, old numbers)
    #q=p.renumber("NOCODE")[0].renumber("RESTART")[0]
    #q.summary()

    #print("Sequence for AlphaFold modeling, with the 20 missing residues replaced by Glycine:")
    #print(">5cil\n"+p.seq(gap="G")+"\n")

    # identify H,L chain residues within 4A to antigen P chain
    binder, target, df_dist=p.rs_around(antigen_chain, dist=dis_threshold)

    # show the distance of binder residues to antigen P chain
    #df_dist[:5]
    #print(binder)
    #print(target)
    #print(df_dist)
    df_dist = pd.DataFrame(df_dist)
    df_dist.to_excel("epitope_paratope_%s.xlsx" % pdb_file)
    return binder, target

    
def get_buried_surface_area(pdb_file, antibody_chains, antigen_chain, file_path):
      # create a new PDB file only containing the antigen and binder residues
    #p=p.extract(binder | "C")
    # save the new structure into a local PDB file
    #p.save("binding.pdb")

    # display the PDB struture, default is show ribbon and color by chains.
    #p.show(show_sidechains=True)
    #all_sasa = p.sasa()
    #antibody_sasa = p.sasa("F:J")
    #antigen_sasa = p.sasa("C")
    #print(file_path + pdb_file + '.pdb')
    p=Protein(file_path + pdb_file + '.pdb')
    #p=Protein(pdb_file)


    all_sasa_list=p.sasa().SASA.values

    chains = (antibody_chains[0] + ":" + antibody_chains[1]) if len(antibody_chains) == 2 else antibody_chains[0]

    antibody_sasa_list=p.sasa(chains).SASA.values
    antigen_sasa_list=p.sasa(antigen_chain).SASA.values

    #p.b_factors(sasa/sasa.max())
    #p.show(style="sphere", color="b")

    complex_sasa = all_sasa_list.sum()
    antibodyHL_sasa = antibody_sasa_list.sum()
    antigenCD_sasa = antigen_sasa_list.sum()
    Buried_surface_area = (antibodyHL_sasa + antigenCD_sasa - complex_sasa)/2.
    return Buried_surface_area


    #df_all = pd.DataFrame(all_sasa)
    #df_antibody = pd.DataFrame(antibody_sasa)
    #df_antigen = pd.DataFrame(antigen_sasa)


    #df_all.to_excel("7df1_sasa.xlsx")
    #df_antibody.to_excel("7df1_FJ_sasa.xlsx")
    #df_antigen.to_excel("7df1_C_sasa.xlsx") 

complex_crystal_list = ['7DF1-1', '8G0M-1', '6S8V-1', '6JMR-1']
antibody_chain_list = [['F', 'J'], ['B'], ['A'], ['E', 'F']]
antigen_chain_list = ['C', 'A', 'B', 'B']

Paratopes = []
Epitopes = []
Buried_surface_areas = []

for i in range(len(complex_crystal_list)):
    pdb_file = complex_crystal_list[i]
    antibody_chains = antibody_chain_list[i]
    antigen_chain = antigen_chain_list[i]
    paratope, epitope = get_epitope_paratope(pdb_file, antibody_chains, antigen_chain, 5., '')
    #print('.....................')
    #rint(paratope)
    #rint(paratope.__str__())
    #print('........................')
    Paratopes.append(paratope.__str__())
    Epitopes.append(epitope.__str__())
    buried_surface_area = get_buried_surface_area(pdb_file, antibody_chains, antigen_chain, '')
    Buried_surface_areas.append(buried_surface_area)

print(Paratopes)
print(Epitopes)
print(Buried_surface_areas)