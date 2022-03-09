# -*- coding: utf-8 -*-
"""
Created on Tue Mar  8 10:52:07 2022

@author: lawashburn
"""

import pandas as pd
import csv


working_directory = input('Input path to working directory, where all results will be stored: ')
#ADvCtrl_path = r"C:\Users\lawashburn\Documents\LiP-MS\Out_update\ADvCtrl_nofilter.csv"
#ADvMCI_path = r"C:\Users\lawashburn\Documents\LiP-MS\Out_update\ADvMCI_nofilter.csv"
#MCIvCtrl_path = r"C:\Users\lawashburn\Documents\LiP-MS\Out_update\MCIvCtrl_nofilter.csv"
#working_directory = r"C:\Users\lawashburn\Documents\LiP-MS\up1.5_down.66"
ADvCtrl_path = input('Enter path to ADvCtrl .csv file updated with p values: ')
ADvMCI_path = input('Enter path to ADvMCI .csv file updated with p values: ')
MCIvCtrl_path = input('Enter path to MCIvCtrl .csv file updated with p values: ')
p_cutoff = input('Enter p-value cutoff (0.05 is recommended): ')
up_thresh = input('Enter upregulated peptide threshold (recommended is 2): ')
down_thresh = input('Enter downregulated peptide threshold (recommended is -0.5): ')
p_cutoff = float(p_cutoff)
up_thresh = float(up_thresh)
down_thresh = float(down_thresh)

p_cutoff = 0.05
up_thresh = 1.5
down_thresh = 2/3

ADvCtrl = pd.read_csv(ADvCtrl_path)
ADvMCI = pd.read_csv(ADvMCI_path)
MCIvCtrl = pd.read_csv(MCIvCtrl_path)

ADvCtrl_mask = ADvCtrl['p value']<= p_cutoff
filtered_ADvCtrl = ADvCtrl[ADvCtrl_mask]


ADvMCI_mask = ADvMCI['p value']<= p_cutoff
filtered_ADvMCI = ADvMCI[ADvMCI_mask]


MCIvCtrl_mask = MCIvCtrl['p value']<= p_cutoff
filtered_MCIvCtrl = MCIvCtrl[MCIvCtrl_mask]

ADvCtrlUp_mask = filtered_ADvCtrl['AD vs Ctrl fold change']>=up_thresh
ADvCtrlUp = filtered_ADvCtrl[ADvCtrlUp_mask]

ADvCtrlDown_mask = filtered_ADvCtrl['AD vs Ctrl fold change']<=down_thresh
ADvCtrlDown = filtered_ADvCtrl[ADvCtrlDown_mask]

ADvMCIUp_mask = filtered_ADvMCI['AD vs MCI fold change']>=up_thresh
ADvMCIUp = filtered_ADvMCI[ADvMCIUp_mask]

ADvMCIDown_mask = filtered_ADvMCI['AD vs MCI fold change']<=down_thresh
ADvMCIDown = filtered_ADvMCI[ADvMCIDown_mask]

MCIvCtrlUp_mask = filtered_MCIvCtrl['MCI vs Ctrl fold change']>=up_thresh
MCIvCtrlUp = filtered_MCIvCtrl[MCIvCtrlUp_mask]

MCIvCtrlDown_mask = filtered_MCIvCtrl['MCI vs Ctrl fold change']<=down_thresh
MCIvCtrlDown = filtered_MCIvCtrl[MCIvCtrlDown_mask]

ADvCtrlUp_out = working_directory + '\\ADvCtrl_upregulated.csv'
ADvMCIUp_out = working_directory + '\\ADvMCI_upregulated.csv'
MCIvCtrlUp_out = working_directory + '\\MCIvCtrl_Upregulated.csv'

ADvCtrlDown_out = working_directory + '\\ADvCtrl_downregulated.csv'
ADvMCIDown_out = working_directory + '\\ADvMCI_downregulated.csv'
MCIvCtrlDown_out = working_directory + '\\MCIvCtrl_downregulated.csv'

with open(ADvCtrlUp_out,'w',newline='') as filed:
    writerd = csv.writer(filed)
    ADvCtrlUp.to_csv(filed,index=False)
    
with open(ADvMCIUp_out,'w',newline='') as filed:
    writerd = csv.writer(filed)
    ADvMCIUp.to_csv(filed,index=False)

with open(MCIvCtrlUp_out,'w',newline='') as filed:
    writerd = csv.writer(filed)
    MCIvCtrlUp.to_csv(filed,index=False)
    
with open(ADvCtrlDown_out,'w',newline='') as filed:
    writerd = csv.writer(filed)
    ADvCtrlDown.to_csv(filed,index=False)
    
with open(ADvMCIDown_out,'w',newline='') as filed:
    writerd = csv.writer(filed)
    ADvMCIDown.to_csv(filed,index=False)

with open(MCIvCtrlDown_out,'w',newline='') as filed:
    writerd = csv.writer(filed)
    MCIvCtrlDown.to_csv(filed,index=False)
    
MCIvCtrlUpSeq = MCIvCtrlUp['Sequence'].values.tolist()
MCIvCtrlDownSeq = MCIvCtrlDown['Sequence'].values.tolist()
#MCIvCtrlDownSeq = ['TYDFHG','KAAL','SVPDHL','DDD']
MCIvCtrlDownHit = []
MCIvsCtrlUpHit = []

for a in MCIvCtrlDownSeq:
    for b in MCIvCtrlUpSeq:
        if b in a:
            MCIvCtrlDownHit.append(a)
            MCIvsCtrlUpHit.append(b)
            
for a in MCIvCtrlUpSeq:
    for b in MCIvCtrlDownSeq:
        if b in a:
            MCIvCtrlDownHit.append(b)
            MCIvsCtrlUpHit.append(a)
            
MCIvCtrl_Hits = pd.DataFrame()
MCIvCtrl_Hits['Upregulated'] = MCIvsCtrlUpHit
MCIvCtrl_Hits['Downregulated'] = MCIvCtrlDownHit
MCIvCtrl_Hits['Sequence'] = MCIvCtrl_Hits['Upregulated']
MCIvCtrl_Hits_merge = pd.merge(MCIvCtrl_Hits, MCIvCtrlUp, on=['Sequence'], how='inner')

MCIvCtrl_out = pd.DataFrame()
MCIvCtrl_out['Upregulated'] = MCIvCtrl_Hits_merge['Upregulated']
MCIvCtrl_out['Downregulated'] = MCIvCtrl_Hits_merge['Downregulated']
MCIvCtrl_out['Protein IDs'] = MCIvCtrl_Hits_merge['Protein IDs']
MCIvCtrl_out['Protein Names'] = MCIvCtrl_Hits_merge['Protein Names']
MCIvCtrl_out['Gene Names'] = MCIvCtrl_Hits_merge['Gene Names']
MCIvCtrl_out['Fasta headers'] = MCIvCtrl_Hits_merge['Fasta headers']

print(MCIvCtrl_out)
MCIvCtrl_hits_out = working_directory + '\\MCIvCtrl_hits.csv'
with open(MCIvCtrl_hits_out,'w',newline='') as filed:
    writerd = csv.writer(filed)
    MCIvCtrl_out.to_csv(filed,index=False)
           
ADvCtrlUpSeq = ADvCtrlUp['Sequence'].values.tolist()
ADvCtrlDownSeq = ADvCtrlDown['Sequence'].values.tolist()

ADvCtrlDownHit = []
ADvsCtrlUpHit = []

for a in ADvCtrlDownSeq:
    for b in ADvCtrlUpSeq:
        if b in a:
            ADvCtrlDownHit.append(a)
            ADvsCtrlUpHit.append(b)
            
for a in ADvCtrlUpSeq:
    for b in ADvCtrlDownSeq:
        if b in a:
            ADvCtrlDownHit.append(b)
            ADvsCtrlUpHit.append(a)

ADvCtrl_Hits = pd.DataFrame()
ADvCtrl_Hits['Upregulated'] = ADvsCtrlUpHit
ADvCtrl_Hits['Downregulated'] = ADvCtrlDownHit

ADvCtrl_Hits['Sequence'] = ADvCtrl_Hits['Upregulated']
ADvCtrl_Hits_merge = pd.merge(ADvCtrl_Hits, ADvCtrlUp, on=['Sequence'], how='inner')

ADvCtrl_out = pd.DataFrame()
ADvCtrl_out['Upregulated'] = ADvCtrl_Hits_merge['Upregulated']
ADvCtrl_out['Downregulated'] = ADvCtrl_Hits_merge['Downregulated']
ADvCtrl_out['Protein IDs'] = ADvCtrl_Hits_merge['Protein IDs']
ADvCtrl_out['Protein Names'] = ADvCtrl_Hits_merge['Protein Names']
ADvCtrl_out['Gene Names'] = ADvCtrl_Hits_merge['Gene Names']
ADvCtrl_out['Fasta headers'] = ADvCtrl_Hits_merge['Fasta headers']

ADvCtrl_hits_out = working_directory + '\\ADvCtrl_hits.csv'
with open(ADvCtrl_hits_out,'w',newline='') as filed:
    writerd = csv.writer(filed)
    ADvCtrl_out.to_csv(filed,index=False)
            
ADvMCIUpSeq = ADvMCIUp['Sequence'].values.tolist()
ADvMCIDownSeq = ADvMCIDown['Sequence'].values.tolist()

ADvMCIDownHit = []
ADvMCIUpHit = []

for a in ADvMCIDownSeq:
    for b in ADvMCIUpSeq:
        if b in a:
            ADvMCIDownHit.append(a)
            ADvMCIUpHit.append(b)
            
for a in ADvMCIUpSeq:
    for b in ADvMCIDownSeq:
        if b in a:
            ADvMCIDownHit.append(b)
            ADvMCIUpHit.append(a)
            
ADvMCI_Hits = pd.DataFrame()
ADvMCI_Hits['Upregulated'] = ADvMCIUpHit
ADvMCI_Hits['Downregulated'] = ADvMCIDownHit

ADvMCI_Hits['Sequence'] = ADvMCI_Hits['Upregulated']
ADvMCI_Hits_merge = pd.merge(ADvMCI_Hits, ADvMCIUp, on=['Sequence'], how='inner')

ADvMCI_out = pd.DataFrame()
ADvMCI_out['Upregulated'] = ADvMCI_Hits_merge['Upregulated']
ADvMCI_out['Downregulated'] = ADvMCI_Hits_merge['Downregulated']
ADvMCI_out['Protein IDs'] = ADvMCI_Hits_merge['Protein IDs']
ADvMCI_out['Protein Names'] = ADvMCI_Hits_merge['Protein Names']
ADvMCI_out['Gene Names'] = ADvMCI_Hits_merge['Gene Names']
ADvMCI_out['Fasta headers'] = ADvMCI_Hits_merge['Fasta headers']

ADvMCI_Hits_out = working_directory + '\\ADvMCI_hits.csv'
with open(ADvMCI_Hits_out,'w',newline='') as filed:
    writerd = csv.writer(filed)
    ADvMCI_out.to_csv(filed,index=False)