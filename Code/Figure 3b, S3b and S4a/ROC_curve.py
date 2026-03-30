'''
This code is used to generate ROC curves to evaluate the predictive performance of four metrics—frequency, enrichment ratio, PSSM, and PSERM—across different antibody libraries, 
assessing how effectively these metrics can predict whether an antibody variant exhibits species cross-reactivity.
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
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from matplotlib.colors import TwoSlopeNorm
from sklearn.metrics import roc_auc_score, roc_curve, RocCurveDisplay, auc
from sklearn.datasets import make_classification
from matplotlib import pyplot as plt

#generate ROC curves for dual-selection 1 and 2 (Figure3b, S3a and S4a)
#here the file 'Variants_FRE_ER_PSSM_PSERM - Copy' combined the data from raw data of Figure3b, S3a and S4a
s1f4_binding = pd.read_csv(os.path.join('/mnt/d/Research_umich/Species-cross-reactivity-paper/Final submission/code', 'Variants_FRE_ER_PSSM_PSERM - Copy.csv'), index_col = 0)


#get the scores of variants based on Dual selection 1 and 2
la = list(s1f4_binding['Cross-reactive label (1 = yes; 0 = no)'])
fre1 = np.array(list(s1f4_binding['Log frequency in the FACS1-enriched library sorted against mouse and human CD98hc']))
fre2 = np.array(list(s1f4_binding['Log frequency in the FACS2-enriched library sorted against mouse and human CD98hc']))
e1 = list(s1f4_binding['Enrichment ratio in the FACS1-enriched library sorted against mouse and human CD98hc'])
e2 = list(s1f4_binding['Enrichment ratio in the FACS2-enriched library sorted against mouse and human CD98hc'])
pssm1 = np.array(list(s1f4_binding['PSSM score in the FACS1-enriched library sorted against mosue and human CD98hc']))
pssm2 = np.array(list(s1f4_binding['PSSM score in the FACS2-enriched library sorted against mosue and human CD98hc']))
pserm1 = np.array(list(s1f4_binding['PSERM score in the FACS1-enriched library sorted against mosue and human CD98hc']))
pserm2 = np.array(list(s1f4_binding['PSERM score in the FACS2-enriched library sorted against mosue and human CD98hc']))

psermm1 = np.array(list(s1f4_binding['PSERM score in the FACS1-enriched library sorted against mouse CD98hc']))
psermh1 = np.array(list(s1f4_binding['PSERM score in the FACS1-enriched library sorted against human CD98hc']))



labels = np.array(la)

label1 = []
er1 = []
label2 = []
er2 = []
for i in range(len(e1)):
    
    if not math.isnan(e1[i]):
        label1.append(la[i])
        er1.append(e1[i])

for i in range(len(e2)):
    
    if not math.isnan(e2[i]):
        label2.append(la[i])
        er2.append(e2[i])

label1 = np.array(label1)
label2 = np.array(label2)
er1 = np.array(er1)
er2 = np.array(er2)


def roc_display(labels, preds):
    fprs = dict()
    tprs = dict()
    roc_aucs = dict()
    for i in range(len(preds)):
        fpr, tpr, thresholds = roc_curve(labels[i], preds[i])
        roc_auc = auc(fpr, tpr)
        fprs[i] = fpr
        tprs[i] = tpr
        roc_aucs[i] = roc_auc
    plt.figure()
    lw = 2
    colors = ['black', 'red', 'blue', 'green', 'yellow', 'pink', 'lightblue', 'orange']
    #names = ['fre1', 'fre2', 'er1', 'er2', 'pssm1', 'pssm2', 'pserm1', 'pserm2']
    names = ['fre', 'er', 'pssm', 'pserm']

    for i in range(len(preds)):
        plt.plot(fprs[i], tprs[i], color=colors[i],
             lw=lw, label='ROC curve of %s (AUC = %0.2f)' % (names[i], roc_aucs[i]))
    plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Receiver Operating Characteristic example')
    bwith = 2.0
    ax = plt.gca() 
    ax.spines['bottom'].set_linewidth(bwith)
    ax.spines['left'].set_linewidth(bwith)
    ax.spines['top'].set_linewidth(bwith)
    ax.spines['right'].set_linewidth(bwith)

    plt.legend(loc="lower right")
    plt.show()
#label_all1 = [labels, labels, label1, label2, labels, labels, labels, labels]  
#preds = [fre1, fre2, er1, er2, pssm1, pssm2, pserm1, pserm2]
label_all1 = [labels, label1, labels, labels]  
preds1 = [fre1, er1, pssm1, pserm1]
label_all2 = [labels, label2, labels, labels]  
preds2 = [fre2, er2, pssm2, pserm2]
label_all3 = [labels, labels, labels]  
preds3 = [pserm1, psermm1, psermh1]

#figure3b
roc_display(label_all1, preds1)
#figures3a
roc_display(label_all2, preds2)
#figures4a
roc_display(label_all3, preds3)


