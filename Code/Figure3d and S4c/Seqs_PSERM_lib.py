#Create theoretical clones and calculate PSERM for each of the clones for the optimized H3 library
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

#Create amino acid tables for each position
AA_H96 = ['F','Y']
AA_H97 = ['W']
AA_H98 = ['D','E']
AA_H99 = ['A','D','E','G','H','I','L','M','N','Q','P','R','S','T','V']
AA_H100 = ['A','F','I','L','V']
AA_H100B = ['A','G','P','S','V']
AA_H100C = ['C','L','F','S','W']
AA_H100D = ['P']
AA_H100E = ['R','S','T', 'P']
AA_H102 = ['A','D','E','G','H','I','L','M','N','Q','P','R','S','T','V']

PSERM_sums = []
for i in range(0,10000):

#Generate a random optimized HCDR3 amino acid sequence and save it
    HCDR3_seq = [0,0,0,0,0,0,0,0,0,0]
    HCDR3_seq[0] = np.random.choice(AA_H96)
    HCDR3_seq[1] = np.random.choice(AA_H97)
    HCDR3_seq[2] = np.random.choice(AA_H98)
    HCDR3_seq[3] = np.random.choice(AA_H99)
    HCDR3_seq[4] = np.random.choice(AA_H100)
    HCDR3_seq[5] = np.random.choice(AA_H100B)
    HCDR3_seq[6] = np.random.choice(AA_H100C)
    HCDR3_seq[7] = np.random.choice(AA_H100D)
    HCDR3_seq[8] = np.random.choice(AA_H100E)
    HCDR3_seq[9] = np.random.choice(AA_H102)
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
PSERM_sums_df.to_csv('HCDR3_PSERMs_average_opt_S2_10000.csv')
#PSERM_sums_df.to_csv('HCDR3_PSERMs_average_opt_S2_10000.csv')


