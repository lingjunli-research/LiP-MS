# -*- coding: utf-8 -*-
"""
Created on Mon Mar  7 15:48:03 2022

@author: lawashburn
"""

import pandas as pd
import numpy as np
import csv
from scipy import stats
import math

working_directory = input('Input path to working directory, where all results will be stored: ')
protein_path = input('Input path to Proteingroup.csv file: ')
peptides_path = input('Input path to Peptidegroup.csv file: ')
protein_path = r"C:\Users\lawashburn\Documents\LiP-MS\15_Proteingroup.csv"
peptides_path = r"C:\Users\lawashburn\Documents\LiP-MS\15_peptidesgroup.csv"
working_directory = r"C:\Users\lawashburn\Documents\LiP-MS\Out_update"

protein = pd.read_csv(protein_path)
peptide = pd.read_csv(peptides_path)
protein = protein.rename(columns={'LFQ intensity MCI_003':'MCI1','LFQ intensity MCI_154':'MCI2','LFQ intensity MCI_242':'MCI3',
                                  'LFQ intensity MCI_219':'MCI4','LFQ intensity MCI_068':'MCI5','LFQ intensity CTRL_019':'Ctrl1',
                                  'LFQ intensity CTRL_275':'Ctrl2','LFQ intensity CTRL_375':'Ctrl3','LFQ intensity CTRL_351':'Ctrl4',
                                  'LFQ intensity CTRL_111':'Ctrl5','LFQ intensity AD_214':'AD1','LFQ intensity AD_082':'AD2',
                                  'LFQ intensity AD_075':'AD3','LFQ intensity AD_078':'AD4','LFQ intensity AD_059':'AD5'})

peptide = peptide.rename(columns={'Protein names':'Protein name'})
AD = pd.DataFrame()
AD['AD1'] = protein['AD1']
AD['AD2'] = protein['AD2']
AD['AD3'] = protein['AD3']
AD['AD4'] = protein['AD4']
AD['AD5'] = protein['AD5']
AD = AD.iloc[1:,:]
AD['mean'] = AD.astype(float).replace(0, np.nan).mean(axis=1, skipna=True)


MCI = pd.DataFrame()
MCI['MCI1'] = protein['MCI1']
MCI['MCI2'] = protein['MCI2']
MCI['MCI3'] = protein['MCI3']
MCI['MCI4'] = protein['MCI4']
MCI['MCI5'] = protein['MCI5']
MCI = MCI.iloc[1:,:]
MCI['mean'] = MCI.astype(float).replace(0, np.nan).mean(axis=1, skipna=True)

Ctrl = pd.DataFrame()
Ctrl['Ctrl1'] = protein['Ctrl1']
Ctrl['Ctrl2'] = protein['Ctrl2']
Ctrl['Ctrl3'] = protein['Ctrl3']
Ctrl['Ctrl4'] = protein['Ctrl4']
Ctrl = Ctrl.iloc[1:,:]
Ctrl['mean'] = Ctrl.astype(float).replace(0, np.nan).mean(axis=1, skipna=True)

protein_calc = pd.DataFrame()
protein_calc['Protein ID'] = protein['Protein IDs']
protein_calc['Protein name'] = protein['Protein names']
protein_calc['Gene name'] = protein['Gene names']
protein_calc['Fasta header'] = protein['Fasta headers']
protein_calc['AD avg'] = AD['mean']
protein_calc['MCI avg'] = MCI['mean']
protein_calc['Ctrl avg'] = Ctrl['mean']
protein_calc['AD/Ctrl fold change'] = protein_calc['AD avg'] / protein_calc['Ctrl avg']
protein_calc['AD/MCI fold change'] = protein_calc['AD avg'] / protein_calc['MCI avg']
protein_calc['MCI/Ctrl fold change'] = protein_calc['MCI avg'] / protein_calc['Ctrl avg']
print(protein_calc)
df_common = pd.merge(protein_calc, peptide, on=['Protein name'], how='inner')

MCIvCtrl_out2 = working_directory + '\\common.csv'
with open(MCIvCtrl_out2,'w',newline='') as filed:
    writerd = csv.writer(filed)
    df_common.to_csv(filed,index=False)


df_common.dropna(subset = ["Protein name"], inplace=True)

df_common = df_common.rename(columns={'LFQ intensity CTRL_351':'Pep Ctrl 1','LFQ intensity CTRL_375':'Pep Ctrl 2','LFQ intensity CTRL_275':'Pep Ctrl 3',
                                  'LFQ intensity CTRL_019':'Pep Ctrl 4','LFQ intensity CTRL_111':'Pep Ctrl 5','LFQ intensity AD_059':'Pep AD 1',
                                  'LFQ intensity AD_214':'Pep AD 2','LFQ intensity AD_075':'Pep AD 3','LFQ intensity AD_082':'Pep AD 4',
                                  'LFQ intensity AD_078':'Pep AD 5','LFQ intensity MCI_154':'Pep MCI 1','LFQ intensity MCI_068':'Pep MCI 2',
                                  'LFQ intensity MCI_219':'Pep MCI 3','LFQ intensity MCI_242':'Pep MCI 4','LFQ intensity MCI_003':'Pep MCI 5'})


df_common = df_common.drop(df_common[df_common['Pep Ctrl 1'] == 'CTRL'].index)
df_common = df_common.drop(df_common[df_common['Pep AD 1'] == 'AD'].index)
df_common = df_common.drop(df_common[df_common['Pep MCI 1'] == 'MCI'].index)

ADvCtrl = pd.DataFrame()
ADvCtrl['AD 1 vs Ctrl normalized'] = (df_common['Pep AD 1'].astype(float)) / (df_common['AD/Ctrl fold change'].astype(float))
ADvCtrl['AD 2 vs Ctrl normalized'] = (df_common['Pep AD 2'].astype(float)) / (df_common['AD/Ctrl fold change'].astype(float))
ADvCtrl['AD 3 vs Ctrl normalized'] = (df_common['Pep AD 3'].astype(float)) / (df_common['AD/Ctrl fold change'].astype(float))
ADvCtrl['AD 4 vs Ctrl normalized'] = (df_common['Pep AD 4'].astype(float)) / (df_common['AD/Ctrl fold change'].astype(float))
ADvCtrl['AD 5 vs Ctrl normalized'] = (df_common['Pep AD 5'].astype(float)) / (df_common['AD/Ctrl fold change'].astype(float))
ADvCtrl['AD vs Ctrl normalized average'] = ADvCtrl.astype(float).replace(0, np.nan).mean(axis=1, skipna=True)
Ctrl = pd.DataFrame()
Ctrl['Pep Ctrl 1'] = df_common['Pep Ctrl 1'].astype(float)
Ctrl['Pep Ctrl 2'] = df_common['Pep Ctrl 2'].astype(float)
Ctrl['Pep Ctrl 3'] = df_common['Pep Ctrl 3'].astype(float)
Ctrl['Pep Ctrl 4'] = df_common['Pep Ctrl 4'].astype(float)
Ctrl['Pep Ctrl 5'] = df_common['Pep Ctrl 5'].astype(float)
Ctrl['Ctrl mean'] = Ctrl.astype(float).replace(0, np.nan).mean(axis=1, skipna=True)
ADvCtrl['Pep Ctrl 1'] = Ctrl['Pep Ctrl 1']
ADvCtrl['Pep Ctrl 2'] = Ctrl['Pep Ctrl 2']
ADvCtrl['Pep Ctrl 3'] = Ctrl['Pep Ctrl 3']
ADvCtrl['Pep Ctrl 4'] = Ctrl['Pep Ctrl 4']
ADvCtrl['Pep Ctrl 5'] = Ctrl['Pep Ctrl 5']
ADvCtrl['Ctrl mean'] = Ctrl['Ctrl mean']

ADvMCI = pd.DataFrame()
ADvMCI['AD 1 vs MCI normalized'] = (df_common['Pep AD 1'].astype(float)) / (df_common['AD/MCI fold change'].astype(float))
ADvMCI['AD 2 vs MCI normalized'] = (df_common['Pep AD 2'].astype(float)) / (df_common['AD/MCI fold change'].astype(float))
ADvMCI['AD 3 vs MCI normalized'] = (df_common['Pep AD 3'].astype(float)) / (df_common['AD/MCI fold change'].astype(float))
ADvMCI['AD 4 vs MCI normalized'] = (df_common['Pep AD 4'].astype(float)) / (df_common['AD/MCI fold change'].astype(float))
ADvMCI['AD 5 vs MCI normalized'] = (df_common['Pep AD 5'].astype(float)) / (df_common['AD/MCI fold change'].astype(float))
ADvMCI['AD vs MCI normalized average'] = ADvMCI.astype(float).replace(0, np.nan).mean(axis=1, skipna=True)
MCI = pd.DataFrame()
MCI['Pep MCI 1'] = df_common['Pep MCI 1'].astype(float)
MCI['Pep MCI 2'] = df_common['Pep MCI 2'].astype(float)
MCI['Pep MCI 3'] = df_common['Pep MCI 3'].astype(float)
MCI['Pep MCI 4'] = df_common['Pep MCI 4'].astype(float)
MCI['Pep MCI 5'] = df_common['Pep MCI 5'].astype(float)
MCI['MCI mean'] = MCI.astype(float).replace(0, np.nan).mean(axis=1, skipna=True)
ADvMCI['Pep MCI 1'] = MCI['Pep MCI 1']
ADvMCI['Pep MCI 2'] = MCI['Pep MCI 2']
ADvMCI['Pep MCI 3'] = MCI['Pep MCI 3']
ADvMCI['Pep MCI 4'] = MCI['Pep MCI 4']
ADvMCI['Pep MCI 5'] = MCI['Pep MCI 5']
ADvMCI['MCI mean'] = MCI['MCI mean']

MCIvCtrl = pd.DataFrame()
MCIvCtrl['MCI 1 vs Ctrl normalized'] = (df_common['Pep MCI 1'].astype(float)) / (df_common['MCI/Ctrl fold change'].astype(float))
MCIvCtrl['MCI 2 vs Ctrl normalized'] = (df_common['Pep MCI 2'].astype(float)) / (df_common['MCI/Ctrl fold change'].astype(float))
MCIvCtrl['MCI 3 vs Ctrl normalized'] = (df_common['Pep MCI 3'].astype(float)) / (df_common['MCI/Ctrl fold change'].astype(float))
MCIvCtrl['MCI 4 vs Ctrl normalized'] = (df_common['Pep MCI 4'].astype(float)) / (df_common['MCI/Ctrl fold change'].astype(float))
MCIvCtrl['MCI 5 vs Ctrl normalized'] = (df_common['Pep MCI 5'].astype(float)) / (df_common['MCI/Ctrl fold change'].astype(float))
MCIvCtrl['MCI vs Ctrl normalized average'] = MCIvCtrl.astype(float).replace(0, np.nan).mean(axis=1, skipna=True)
Ctrl = pd.DataFrame()
Ctrl['Pep Ctrl 1'] = df_common['Pep Ctrl 1'].astype(float)
Ctrl['Pep Ctrl 2'] = df_common['Pep Ctrl 2'].astype(float)
Ctrl['Pep Ctrl 3'] = df_common['Pep Ctrl 3'].astype(float)
Ctrl['Pep Ctrl 4'] = df_common['Pep Ctrl 4'].astype(float)
Ctrl['Pep Ctrl 5'] = df_common['Pep Ctrl 5'].astype(float)
Ctrl['Ctrl mean'] = Ctrl.astype(float).replace(0, np.nan).mean(axis=1, skipna=True)
MCIvCtrl['Pep Ctrl 1'] = Ctrl['Pep Ctrl 1']
MCIvCtrl['Pep Ctrl 2'] = Ctrl['Pep Ctrl 2']
MCIvCtrl['Pep Ctrl 3'] = Ctrl['Pep Ctrl 3']
MCIvCtrl['Pep Ctrl 4'] = Ctrl['Pep Ctrl 4']
MCIvCtrl['Pep Ctrl 5'] = Ctrl['Pep Ctrl 5']
MCIvCtrl['Ctrl mean'] = Ctrl['Ctrl mean']

Ctrl = pd.DataFrame()
Ctrl['Pep Ctrl 1'] = df_common['Pep Ctrl 1'].astype(float)
Ctrl['Pep Ctrl 2'] = df_common['Pep Ctrl 2'].astype(float)
Ctrl['Pep Ctrl 3'] = df_common['Pep Ctrl 3'].astype(float)
Ctrl['Pep Ctrl 4'] = df_common['Pep Ctrl 4'].astype(float)
Ctrl['Pep Ctrl 5'] = df_common['Pep Ctrl 5'].astype(float)
Ctrl['Ctrl average'] = Ctrl.astype(float).replace(0, np.nan).mean(axis=1, skipna=True)

ADvCtrl['Protein IDs'] = df_common['Protein ID']
ADvCtrl['Protein Names'] = df_common['Protein name']
ADvCtrl['Gene Names'] = df_common['Gene names']
ADvCtrl['Fasta headers'] = df_common['Fasta header']
ADvCtrl['Sequence'] = df_common['Sequence']

ADvMCI['Protein IDs'] = df_common['Protein ID']
ADvMCI['Protein Names'] = df_common['Protein name']
ADvMCI['Gene Names'] = df_common['Gene names']
ADvMCI['Fasta headers'] = df_common['Fasta header']
ADvMCI['Sequence'] = df_common['Sequence']

MCIvCtrl['Protein IDs'] = df_common['Protein ID']
MCIvCtrl['Protein Names'] = df_common['Protein name']
MCIvCtrl['Gene Names'] = df_common['Gene names']
MCIvCtrl['Fasta headers'] = df_common['Fasta header']
MCIvCtrl['Sequence'] = df_common['Sequence']

ADvCtrl['AD vs Ctrl fold change'] = ADvCtrl['AD vs Ctrl normalized average'] / Ctrl['Ctrl average']
ADvMCI['AD vs MCI fold change'] = ADvMCI['AD vs MCI normalized average'] / Ctrl['Ctrl average']
MCIvCtrl['MCI vs Ctrl fold change'] = MCIvCtrl['MCI vs Ctrl normalized average'] / Ctrl['Ctrl average']

#ADvCtrl.dropna(subset = ['AD vs Ctrl fold change'], inplace=True)
ADvMCI.dropna(subset = ['AD vs MCI fold change'], inplace=True)
MCIvCtrl.dropna(subset = ['MCI vs Ctrl fold change'], inplace=True)

ADvCtrl_out = working_directory + '\\ADvCtrl_nofilter.csv' 
ADvMCI_out = working_directory + '\\ADvMCI_nofilter.csv' 
MCIvCtrl_out = working_directory + '\\MCIvCtrl_nofilter.csv' 

with open(ADvCtrl_out,'w',newline='') as filed:
    writerd = csv.writer(filed)
    ADvCtrl.to_csv(filed,index=False)
    
with open(ADvMCI_out,'w',newline='') as filed:
    writerd = csv.writer(filed)
    ADvMCI.to_csv(filed,index=False)

with open(MCIvCtrl_out,'w',newline='') as filed:
    writerd = csv.writer(filed)
    MCIvCtrl.to_csv(filed,index=False)

    
print('Unfiltered data has been exported to working directory')