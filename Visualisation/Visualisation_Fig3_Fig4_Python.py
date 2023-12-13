#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  4 11:34:14 2023

@author: riannehaartsen
"""

# Import libraries
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Number of clean trials 1
Pipeline_Ntrls = pd.read_csv('Pipeline_InclRates.csv', sep = ',', header=0) 

Roche_Ntrls = pd.read_csv('Roche1_InclRates.csv', sep = ',', header=0) 
MADE_Ntrls = pd.read_csv('MADE_InclRates.csv', sep = ',', header=0) 
BOND_Ntrls = pd.read_csv('BOND_InclRates.csv', sep = ',', header=0) 
HAPPE_Ntrls = pd.read_csv('HAPPE_InclRates.csv', sep = ',', header=0) 


# Numbers of participants included
PL_names = ('Roche', 'MADE', 'BOND', 'HAPPE')
pd.value_counts(Roche_Ntrls['Ntrls_tot'] == 0)
pd.value_counts((Roche_Ntrls['Ntrls_tot'] > 0) & (Roche_Ntrls['Ntrls_tot'] < 20))
pd.value_counts((Roche_Ntrls['Ntrls_tot'] >= 20) & (Roche_Ntrls['Ntrls_tot'] < 90))
pd.value_counts(Roche_Ntrls['Ntrls_tot'] >= 90)

pd.value_counts(MADE_Ntrls['Ntrls_tot'] == 0)
pd.value_counts((MADE_Ntrls['Ntrls_tot'] > 0) & (MADE_Ntrls['Ntrls_tot'] < 20))
pd.value_counts((MADE_Ntrls['Ntrls_tot'] >= 20) & (MADE_Ntrls['Ntrls_tot'] < 90))
pd.value_counts(MADE_Ntrls['Ntrls_tot'] >= 90)

pd.value_counts(BOND_Ntrls['Ntrls_tot'] == 0)
pd.value_counts((BOND_Ntrls['Ntrls_tot'] > 0) & (BOND_Ntrls['Ntrls_tot'] < 20))
pd.value_counts((BOND_Ntrls['Ntrls_tot'] >= 20) & (BOND_Ntrls['Ntrls_tot'] < 90))
pd.value_counts(BOND_Ntrls['Ntrls_tot'] >= 90)

pd.value_counts(HAPPE_Ntrls['Neps_tot'] == 0)
pd.value_counts((HAPPE_Ntrls['Neps_tot'] > 0) & (HAPPE_Ntrls['Neps_tot'] < 20))
pd.value_counts((HAPPE_Ntrls['Neps_tot'] >= 20) & (HAPPE_Ntrls['Neps_tot'] < 90))
pd.value_counts(HAPPE_Ntrls['Neps_tot'] >= 90)


Data_counts = {'Threshold for power and connectivity (90 trls)': np.array([129, 124, 126, 122]),
    'Threshold for power (20 trls)': np.array([0, 3, 1, 1]),
    'Not enough data': np.array([0, 2, 0, 0]),
    'No data': np.array([2, 2, 4, 8])}

width = 0.6

fig, ax = plt.subplots()
bottom = np.zeros(4)

for Incl, Incl_count in Data_counts.items():
    p = ax.bar(PL_names, Incl_count, width, label=Incl, bottom=bottom)
    bottom += Incl_count
    
ax.set_title('Data thresholds for inclusion: all conditions')
ax.legend(loc='best')

plt.show()  
plt.tight_layout()  



# Numbers of participants included: condition differences
PL_names = ('Roche', 'MADE', 'BOND', 'HAPPE')
pd.value_counts((Roche_Ntrls['Ntrls_soc'] == 0) & (Roche_Ntrls['Ntrls_toy'] == 0))
pd.value_counts((Roche_Ntrls['Ntrls_soc'] >= 20) & (Roche_Ntrls['Ntrls_soc'] < 90) &
                (Roche_Ntrls['Ntrls_toy'] >= 20) & (Roche_Ntrls['Ntrls_toy'] < 90))
pd.value_counts((Roche_Ntrls['Ntrls_soc'] >= 90) & (Roche_Ntrls['Ntrls_toy'] >= 90))

pd.value_counts((MADE_Ntrls['Ntrls_soc'] == 0) & (MADE_Ntrls['Ntrls_toy'] == 0))
pd.value_counts((MADE_Ntrls['Ntrls_soc'] >= 20) & (MADE_Ntrls['Ntrls_soc'] < 90) &
                (MADE_Ntrls['Ntrls_toy'] >= 20) & (MADE_Ntrls['Ntrls_toy'] < 90))
pd.value_counts((MADE_Ntrls['Ntrls_soc'] >= 90) & (MADE_Ntrls['Ntrls_toy'] >= 90))

pd.value_counts((BOND_Ntrls['Ntrls_soc'] == 0) & (BOND_Ntrls['Ntrls_toy'] == 0))
pd.value_counts((BOND_Ntrls['Ntrls_soc'] >= 20) & (BOND_Ntrls['Ntrls_soc'] < 90) &
                (BOND_Ntrls['Ntrls_toy'] >= 20) & (BOND_Ntrls['Ntrls_toy'] < 90))
pd.value_counts((BOND_Ntrls['Ntrls_soc'] >= 90) & (BOND_Ntrls['Ntrls_toy'] >= 90))

pd.value_counts((HAPPE_Ntrls['Neps_soc'] == 0) & (HAPPE_Ntrls['Neps_toy'] == 0))
pd.value_counts((HAPPE_Ntrls['Neps_soc'] >= 20) & (HAPPE_Ntrls['Neps_soc'] < 90) &
                (HAPPE_Ntrls['Neps_toy'] >= 20) & (HAPPE_Ntrls['Neps_toy'] < 90))
pd.value_counts((HAPPE_Ntrls['Neps_soc'] >= 90) & (HAPPE_Ntrls['Neps_toy'] >= 90))

# 131-128-0-2
# 131-114-7-2
# 131-121-3-4
# 131-121-0-8

Data_counts = {'Threshold for power and connectivity (90 trls)': np.array([128, 114, 121, 121]),
    'Threshold for power (20 trls)': np.array([0, 7, 3, 0]),
    'Not enough data in both conditions': np.array([1, 8, 3, 2]),
    'No data for either condition': np.array([2, 2, 4, 8])}

width = 0.6

fig, ax = plt.subplots()
bottom = np.zeros(4)

for Incl, Incl_count in Data_counts.items():
    p = ax.bar(PL_names, Incl_count, width, label=Incl, bottom=bottom)
    bottom += Incl_count
    
ax.set_title('Data thresholds for inclusion: both conditions')
ax.legend(loc='best')

plt.show()  
plt.tight_layout()  

# Correlation matrices: comparison between pipelines
# Power - all trials
Pow_ICCs_r = pd.read_csv('Power_ICCs_alltrls_rvals.csv', sep = ',', header=None) 
r_values = Pow_ICCs_r.to_numpy()

mask = np.triu(np.ones_like(r_values, dtype=bool)) # this masks the top part of the heatmap
cmap2 = sns.color_palette('Spectral', as_cmap=True)
Measures = ['Roche \u03B4', 'Roche \u03B8', 'Roche \u03B1', 'Roche \u03B2',
             'MADE \u03B4', 'MADE \u03B8', 'MADE \u03B1', 'MADE \u03B2',
             'BOND \u03B4', 'BOND \u03B8', 'BOND \u03B1', 'BOND \u03B2',
             'HAPPE \u03B4', 'HAPPE \u03B8', 'HAPPE \u03B1', 'HAPPE \u03B2']

# Set the theme for the plot and plot it (make sure to execute all the following lines together)
fig, ax = plt.subplots()
sns.set_theme(style="white", font='DejaVu Sans', font_scale=1.2, rc={'figure.figsize':(10, 7)})
sns.heatmap(np.round(r_values, decimals=2), # make sure to round the correlations to fit the values within the heatmap boxes
            vmin=-1.0, vmax=1.0, mask = mask, square = True, center = 0, # vmin and vmax set the minimum and maximum values for the scale
            cmap=cmap2, annot=True, linewidths=.5,
            cbar_kws={"shrink": .8},
            annot_kws={'fontsize': 10},
            xticklabels=Measures,
            yticklabels=Measures)
ax.set_title('ICCs for power across all trials')
plt.tight_layout()

 
# Connectivity - all trials
FC_ICCs_r = pd.read_csv('FC_ICCs_alltrls_rvals.csv', sep = ',', header=None) 
r_values = FC_ICCs_r.to_numpy()
fig, ax = plt.subplots()
sns.set_theme(style="white", font='DejaVu Sans', font_scale=1.2, rc={'figure.figsize':(10, 7)})
sns.heatmap(np.round(r_values, decimals=2), # make sure to round the correlations to fit the values within the heatmap boxes
            vmin=-1.0, vmax=1.0, mask = mask, square = True, center = 0, # vmin and vmax set the minimum and maximum values for the scale
            cmap=cmap2, annot=True, linewidths=.5,
            cbar_kws={"shrink": .8},
            annot_kws={'fontsize': 10},
            xticklabels=Measures,
            yticklabels=Measures)
ax.set_title('ICCs for connectivity across all trials')
plt.tight_layout()

# Power - condition differences
Pow_ICCs_r = pd.read_csv('Pow_ICCs_diffs_rvals.csv', sep = ',', header=None) 
r_values = Pow_ICCs_r.to_numpy()
fig, ax = plt.subplots()
sns.set_theme(style="white", font='DejaVu Sans', font_scale=1.2, rc={'figure.figsize':(10, 7)})
sns.heatmap(np.round(r_values, decimals=2), # make sure to round the correlations to fit the values within the heatmap boxes
            vmin=-1.0, vmax=1.0, mask = mask, square = True, center = 0, # vmin and vmax set the minimum and maximum values for the scale
            cmap=cmap2, annot=True, linewidths=.5,
            cbar_kws={"shrink": .8},
            annot_kws={'fontsize': 10},
            xticklabels=Measures,
            yticklabels=Measures)
ax.set_title('ICCs for power differences between conditions')
plt.tight_layout()

# Connectivity - condition differences
FC_ICCs_r = pd.read_csv('FC_ICCs_diff_rvals.csv', sep = ',', header=None) 
r_values = FC_ICCs_r.to_numpy()
fig, ax = plt.subplots()
sns.set_theme(style="white", font='DejaVu Sans', font_scale=1.2, rc={'figure.figsize':(10, 7)})
sns.heatmap(np.round(r_values, decimals=2), # make sure to round the correlations to fit the values within the heatmap boxes
            vmin=-1.0, vmax=1.0, mask = mask, square = True, center = 0, # vmin and vmax set the minimum and maximum values for the scale
            cmap=cmap2, annot=True, linewidths=.5,
            cbar_kws={"shrink": .8},
            annot_kws={'fontsize': 10},
            xticklabels=Measures,
            yticklabels=Measures)
ax.set_title('ICCs for connectivity differences between conditions')
plt.tight_layout()



# Correlation matrices: split-half reliability
# Power - all trials
Pow_corrs_r = pd.read_csv('SplitHalf_all_rvals.csv', sep = ',', header=None) 
r_values = Pow_corrs_r.to_numpy()

cmap2 = sns.color_palette('Spectral', as_cmap=True)
Pipelines = ['Roche', 'MADE', 'BOND', 'HAPPE']
FreqBand = ['\u03B4', '\u03B8', '\u03B1', '\u03B2',]

# Set the theme for the plot and plot it (make sure to execute all the following lines together)
fig, ax = plt.subplots()
sns.set_theme(style="white", font='DejaVu Sans', font_scale=1.2, rc={'figure.figsize':(10, 7)})
sns.heatmap(np.round(r_values, decimals=2), # make sure to round the correlations to fit the values within the heatmap boxes
            vmin=-1.0, vmax=1.0, center = 0, # vmin and vmax set the minimum and maximum values for the scale
            cmap=cmap2, annot=True, linewidths=.5,
            cbar_kws={"shrink": .8},
            annot_kws={'fontsize': 10},
            xticklabels=FreqBand,
            yticklabels=Pipelines)
ax.set_title('Split-half reliability for power across all trials')
plt.tight_layout()

# Power - condition differences
Pow_corrs_r = pd.read_csv('SplitHalf_diff_rvals.csv', sep = ',', header=None) 
r_values = Pow_corrs_r.to_numpy()

cmap2 = sns.color_palette('Spectral', as_cmap=True)
Pipelines = ['Roche', 'MADE', 'BOND', 'HAPPE']
FreqBand = ['\u03B4', '\u03B8', '\u03B1', '\u03B2',]

# Set the theme for the plot and plot it (make sure to execute all the following lines together)
fig, ax = plt.subplots()
sns.set_theme(style="white", font='DejaVu Sans', font_scale=1.2, rc={'figure.figsize':(10, 7)})
sns.heatmap(np.round(r_values, decimals=2), # make sure to round the correlations to fit the values within the heatmap boxes
            vmin=-1.0, vmax=1.0, center = 0, # vmin and vmax set the minimum and maximum values for the scale
            cmap=cmap2, annot=True, linewidths=.5,
            cbar_kws={"shrink": .8},
            annot_kws={'fontsize': 10},
            xticklabels=FreqBand,
            yticklabels=Pipelines)
ax.set_title('Split-half reliability for power differences between conditions')
plt.tight_layout()
