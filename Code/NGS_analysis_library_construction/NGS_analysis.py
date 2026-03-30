'''
This code was developed to analyze next-generation sequencing data from antibody libraries generated during the directed evolution workflow. 
It processes sequencing-derived variant information to calculate multiple sequence- and enrichment-based metrics, including variant frequency, enrichment ratio, position-specific scoring matrices (PSSM), 
and position-specific enrichment ratio matrices (PSERM). These metrics are then used to score antibody variants and to predict properties related to species cross-reactivity. 
In addition, the code was used to guide the rational design of a redesigned HCDR3 library, which was subsequently employed for the selection and enrichment of antibody variants 
with improved cross-species reactivity
library: 
input: MACS3-enriched library (input library) sorted against mouse CD98hc
f1mp: FACS1-enriched library sorted against mouse CD98hc
f1hp: FACS1-enriched library sorted against human CD98hc
f2mp: FACS2-enriched library sorted against mouse CD98hc
f2hp: FACS2-enriched library sorted against human CD98hc
df1mhp: FACS1-enriched library sorted against mouse and human CD98hc
df2mhp: FACS2-enriched library sorted against mouse and human CD98hc
'''

#Add new directory to path to import ngs and pserm
import sys
import os

path = os.getcwd()
path_with_ngs_and_pserm = os.path.dirname(path)
sys.path.append(path_with_ngs_and_pserm)
import math

#Imports 
import tqdm
import glob
import matplotlib
import numpy as np
import pandas as pd
import seaborn as sns
from ngs import NGS_round_data
from pserm import ngs_analysis, generate_clone_set
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from matplotlib.colors import TwoSlopeNorm

#specify font properties for better export into adobe illustrator

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

#Add myriad pro font
fpath = "/home/adduser/PSERM_main/Library/Fonts/Myriad Pro Regular.ttf"
aapath = '/home/adduser/PSERM_main/Library/Fonts/cour.ttf'
data_path = '/mnt/d/NGS/S1F4-NGS-DATA/NGS2/'

prop = fm.FontProperties(fname=fpath, size = 30)
tickprop = fm.FontProperties(fname = fpath, size = 25)
aaprop = fm.FontProperties(fname = aapath, size = 25)


#template DNA sequence in HCDR3 library

s1f4_wt = 'YYDILGSRPI'

#possible mutants at each position without WT resdiues
s1f4_muts_no_wt = {
0: ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W'],
1: ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W'],
2: ['A', 'C', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'],
3: ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'],
4: ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'],
5: ['A', 'C', 'D', 'E', 'F', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'],
6: ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'T', 'V', 'W', 'Y'],
7: ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'S', 'T', 'V', 'W', 'Y'],
8: ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'],
9: ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'],
}

# Generate an NGS_round_data object by importing sequence:frequency data for each sample,
# which will be used for downstream enrichment ratio calculations.
s1f4_r = NGS_round_data(

    Round = 4, 
    sequence_type = 'mutations', 
    samples =[
       'input', 'f1mp', 'f1mn', 'f1hp', 'f1hn', 'f2mp', 'f2hp', 'df1mhp', 'df2mhp', 'df2mp', 'df2hp', 'df2mhn', 'ovap', 'ovan'],
    sample_of_interest = 'input', 
    path = os.path.join(data_path, 'cross-data'), 
    wild_type = s1f4_wt, 
    mutations_dict = s1f4_muts_no_wt
)


# Generate the union of all observed clones across the selected samples.
s1f4_clone_set = generate_clone_set(s1f4_r, ['input', 'f1mp', 'f1mn', 'f1hp', 'f1hn', 'f2mp', 'f2hp', 'df1mhp', 'df2mhp', 'df2mp', 'df2hp', 'df2mhn', 'ovap', 'ovan'])



s1f4_clone_set_trimmed = list(set(s1f4_clone_set))

# Initialize the NGS analysis object using the trimmed clone set.
s1f4_data = ngs_analysis([s1f4_r], ['input', 'f1mp', 'f1mn', 'f1hp', 'f1hn', 'f2mp', 'f2hp', 'df1mhp', 'df2mhp', 'df2mp', 'df2hp', 'df2mhn', 'ovap', 'ovan'], clone_set = s1f4_clone_set_trimmed)

# Generate the count/frequency matrix D for all clones across samples.
s1f4_data.generate_D()

# Generate and export the PSSM/PSERM matrix for each sample.
for sample in tqdm.tqdm(s1f4_data.samples):
    s1f4_data.generate_PSSM(sample)
    s1f4_data.PSSM[sample].to_excel(os.path.join(data_path, 'cross-data/R4/PSSM', f'{sample}.xlsx'))


for sample in ['f1mp', 'f1mn', 'f1hp', 'f1hn', 'df1mhp', 'f2mp', 'f2hp', 'df2mhp', 'df2mp', 'df2hp', 'df2mhn']:
   s1f4_data.generate_PSERM(In_sample = 'input', Out_sample = sample) 
   s1f4_data.PSERM[sample].to_excel(os.path.join(data_path, 'cross-data/R4/PSERM', f'{sample}.xlsx'))

# Score all clones using the PSSM/PSERM model.
s1f4_data.score_all_clones_mp(method = 'PSSM')
s1f4_data.score_all_clones_mp(method = 'PSERM')


# Visualize the PSERM matrix as a heatmap.

def draw_PSERM_Matrix(s1f4_data, sample): 
    fig, axs = plt.subplots(1, 1, figsize = (8.0, 6))

    aas = sorted(s1f4_data.AA_order)

    mask = np.ones_like(s1f4_data.PSERM[sample])
    for j in s1f4_data.library.keys():
        for aa in s1f4_data.library[j]:
            i = aas.index(aa)
            mask[i, j] = 0

    rows_to_drop = []
    for i in range(mask.shape[0]):
        if mask[i, :].mean() == 1:
            rows_to_drop.append(aas[i])

    PSERM = s1f4_data.PSERM[sample].drop(rows_to_drop, axis = 'index').sort_index()

    updated_mask = np.zeros_like(PSERM)
    row_ind = 0
    for i in range(mask.shape[0]):
        if mask[i, :].mean() != 1:
            updated_mask[row_ind, :] = mask[i, :]
            row_ind += 1

    sns.heatmap(PSERM, robust = True, cmap = 'bwr', ax = axs, center = 0, annot = False)

    #YYDILGSRPI   
    xlabels = ['H96', 'H97', 'H98', 'H99', 'H100', 'H100B', 'H100C', 'H100D', 'H100E', 'H102']

    axs.set_xticks([i + 0.5 for i in range(len(xlabels))])
    axs.set_xticklabels(xlabels, rotation = 90)
    
    colorbar = axs.collections[0].colorbar
    colorbar.set_label(f'amino acid enrichment ratio', fontproperties = tickprop)
    for tick in colorbar.ax.get_yticklabels():
        tick.set_fontproperties(tickprop)

    axs.set_facecolor('lightgrey')
    for tick in axs.get_yticklabels():
        tick.set_fontproperties(aaprop)
        tick.set_rotation(0)
    for tick in axs.get_xticklabels():
        tick.set_fontproperties(tickprop)

    axs.set_title(f'{sample} PSERM Score', fontproperties = prop)


    plt.tight_layout()
    plt.savefig(os.path.join(data_path, 'cross-data/R4/new_data/', f'{sample}_PSERM_Score.pdf'), bbox_inches = 'tight')

    plt.show()



samples = ['df1mhp', 'df2mhp']
# Generate PSERM heatmap visualizations for the selected samples (Figure3c and S4b)
#for sample in samples:
#    draw_PSERM_Matrix(s1f4_data, sample)

def generate_for_data(s1f4_data, file):
    s1f4_binding = pd.read_csv(file, index_col=0)

    sample_dict = {
        'input': 'in the MACS3-enriched library (input library) sorted against mouse CD98hc',
        'f1mp': 'in the FACS1-enriched library sorted against mouse CD98hc',
        'f1hp': 'in the FACS1-enriched library sorted against human CD98hc',
        'df1mhp': 'in the FACS1-enriched library sorted against mouse and human CD98hc',
        'f2mp': 'in the FACS2-enriched library sorted against mouse CD98hc',
        'f2hp': 'in the FACS2-enriched library sorted against human CD98hc',
        'df2mhp': 'in the FACS2-enriched library sorted against mouse and human CD98hc'
    }

    for idx, seq in s1f4_binding['HCDR3_seqs'].items():
        if pd.isna(seq):
            continue
        seq = str(seq)

        # Convert the original HCDR3 sequence into the trimmed sequence format
        if len(seq) > 11:
            seq1 = seq[:5] + seq[6:10] + seq[-1]
        elif len(seq) == 11:
            seq1 = seq[:9] + seq[-1]
        else:
            seq1 = seq

        # Retrieve and store frequency values for each library
        for sample in ['input', 'f1mp', 'f1hp', 'f2mp', 'f2hp', 'df1mhp', 'df2mhp']:
            try:
                fo = s1f4_data.D.loc[seq1, sample]
                s1f4_binding.loc[idx, f'Frequency {sample_dict[sample]}'] = fo
            except:
                pass

        # Calculate and store enrichment ratio, PSSM score, and PSERM score
        for sample in [['input', 'f1mp'], ['input', 'f1hp'], ['input', 'df1mhp'], ['input', 'df2mhp'], ['input', 'f2mp'], ['input', 'f2hp']]:
            try:
                fo = s1f4_data.D.loc[seq1, sample[1]]
                fi = s1f4_data.D.loc[seq1, sample[0]]

                # Calculate log2 enrichment ratio if both frequencies are non-zero
                if fo != 0 and fi != 0:
                    s1f4_binding.loc[idx, f'Enrichment ratio {sample_dict[sample[1]]}'] = np.log2(fo / fi)
                else:
                    s1f4_binding.loc[idx, f'Enrichment ratio {sample_dict[sample[1]]}'] = 'N/A'
            except:
                pass

            try:
                s1f4_binding.loc[idx, f'PSSM score {sample_dict[sample[1]]}'] = s1f4_data.scores.loc[seq1, f'{sample[1]} PSSM Score']
            except:
                s1f4_binding.loc[idx, f'PSSM score {sample_dict[sample[1]]}'] = 'N/A'

            try:
                s1f4_binding.loc[idx, f'PSERM score {sample_dict[sample[1]]}'] = s1f4_data.scores.loc[seq1, f'{sample[1]} PSERM Score']
            except:
                s1f4_binding.loc[idx, f'PSERM score {sample_dict[sample[1]]}'] = 'N/A'

    # Replace all remaining NaN values with 0
    s1f4_binding = s1f4_binding.fillna(0)

    # Save the updated DataFrame back to the CSV file
    s1f4_binding.to_csv(file)

    return s1f4_binding

# generate frequency, enrichment ratio, PSSM, PSERM data for each variants in csv file.
s1f4_binding = generate_for_data(s1f4_data, os.path.join('/mnt/d/Research_umich/Species-cross-reactivity-paper/Final submission/code/', 'Variants_FRE_ER_PSSM_PSERM.csv'))


# Redesign the HCDR3 library.
# Select sequences with df1mhp PSERM score >= 1.686468629 for amino acid distribution analysis.
Seq_high_PSERM = s1f4_data.scores[s1f4_data.scores.loc[:, 'df1mhp PSERM Score'] >= 1.686468629].index


# Define the list of datasets/conditions used to identify sequences observed in the screening data.
filterlist = ['input', 'f1mp', 'f1hp', 'f2mp', 'f2hp', 'df1mhp', 'df2mhp']
#filterlist = ['input', 'f1mp', 'f1mn', 'f1hp', 'f1hn', 'f2mp', 'f2hp', 'df1mhp', 'df2mhp', 'df2mp', 'df2hp', 'df2mhn', 'ovap', 'ovan']


s1f4_data_seq = set()

for fi in filterlist:

    each_seq = s1f4_data.D[s1f4_data.D.loc[:, fi] != 0.0].index
    print(len(each_seq))

    s1f4_data_seq = s1f4_data_seq.union(set(each_seq))

Seq_high_PSERM_1 = list(set(Seq_high_PSERM))
# Retain only high-PSERM sequences that are also present in the selected screening datasets/conditions.
Seq_high_PSERM_2= list(set(Seq_high_PSERM_1).intersection(set(s1f4_data_seq)))
# Calculate the amino acid frequency distribution at each position.
def get_aa_pos(seqs):
    aa_pos = [{}, {}, {}, {}, {}, {}, {}, {}, {}, {}]
    for seq in seqs:
        for i in range(len(seq)):
            if seq[i] not in aa_pos[i].keys():
                aa_pos[i][seq[i]] = 1
            else:
                aa_pos[i][seq[i]] += 1
    for i in range(len(aa_pos)):
        aa_pos[i] =  dict(sorted(aa_pos[i].items(), key=lambda x: x[1], reverse=True))
    for dic in aa_pos:
        sum = 0
        for i in dic.values():
            sum += i
        for key in dic.keys():
            dic[key] = round(dic[key] / sum * 1.0, 4)

    return aa_pos



# Select amino acids at each position until the cumulative frequency exceeds the specified threshold.
def get_aa_pos_prob(lib_con, prob, all_num):
    pos = 10
    d_dict = {}
    for i in range(pos):
        d_dict[i] = []
    
    for i in range(len(lib_con)):
        sum = 0.0
        for key, val in lib_con[i].items():
           if sum > prob:
               break 
           else:
               sum += val/all_num
               d_dict[i].append(key)
    return d_dict

# Compute amino acid frequency distributions across the selected high-PSERM sequences.
lib_con = get_aa_pos(Seq_high_PSERM_2)


#redesign_lib = get_aa_pos_prob(lib_con, 0.90, len(Seq_high_PSERM_2))


print(lib_con)


