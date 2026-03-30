#Create theoretical clones and calculate PSERM for each of the clones for the non-optimized H3 library
#Import packages

import numpy as np
import pandas as pd
import matplotlib as plt
import random

#Load in the PSERM data
#for FACS1 data
PSERM_table = pd.read_excel('df1mhp.xlsx')
# for FACS2 data
#PSERM_table = pd.read_excel('df2mhp.xlsx')

PSERM_table = PSERM_table.sort_values('Amino Acid')
Amino_Acids = ['A','C','D','E','F','G','H','I','K', 'L', 'M', 'N', 'P', 'Q', 'R','S','T','V','W','Y']

PSERM_sums = []
for i in range(0,10000):

#Generate a random non-optimized HCDR3 amino acid sequence and save it
    HCDR3_seq = []
    for x in np.arange(0,10):
        HCDR3_seq.append(np.random.choice(Amino_Acids))
    HCDR3_seq_df = pd.DataFrame({'Amino Acid':HCDR3_seq})
    HCDR3_PSERM_blank = HCDR3_seq_df.assign(Position = ['H96','H97','H98','H99','H100','H100B','H100C','H100D','H100E','H102'])

    #Obtain PSERM values for every amino acid in the generated sequence and save it

    HCDR3_PSERM_blank_merged = HCDR3_PSERM_blank.merge(PSERM_table)
    HCDR3_PSERMs_only = HCDR3_PSERM_blank_merged.iloc[:,2:]
    diag_vals = np.diag(HCDR3_PSERMs_only)
    HCDR3_PSERM_complete = HCDR3_PSERM_blank.assign(PSERM = diag_vals)
    sums = HCDR3_PSERM_complete['PSERM'].sum()
    average_PSERM = sums/10
    PSERM_sums.append(average_PSERM)

PSERM_sums_df = pd.DataFrame({'PSERM average':PSERM_sums})
PSERM_sums_df.to_csv('HCDR3_PSERMs_average_unopt_S1_10000.csv')
#PSERM_sums_df.to_csv('HCDR3_PSERMs_average_unopt_S2_10000.csv')




