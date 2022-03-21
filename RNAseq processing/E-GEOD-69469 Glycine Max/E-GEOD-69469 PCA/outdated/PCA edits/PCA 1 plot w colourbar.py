# -*- coding: utf-8 -*-
"""
Created on Tue May 22 16:24:38 2018

@author: weixiong001
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from sklearn import decomposition
from sklearn import preprocessing  # don't need this if data is already scaled
import csv
'''
This script takes one .txt file, performs PCA analysis on it and generatres a
PCA plot in .png format.
Requires colourmap_file, which contains sample's timepoint information along with
its associated value in the colourbar. This is needed for colourbar generation
'''

input_file = 'E-GEOD-69469_PCA_m1.txt'
colourmap_file = 'E-GEOD-69469_colourbar.csv'
file_name = input_file.split('.')[0]

# creates colormap_dict, which has as keys, day/night timepoint information (e.g.
# L1, L3, L5.5....), and has as values, the keys' corresponding number in a 24h day. 
# E.g. L1 is a value of 1, L3 has a value of 3, and L5.5 has a value of 5.5 
colormap_dict = {}
with open(colourmap_file, 'r', newline='') as csvfile:
    reader = csv.reader(csvfile, delimiter='\t')
    next(reader)
    for row in reader:
        new_row = row[0].split(',')
        colormap_dict[new_row[0]] = new_row[1]

###############################################################################
# data preprocessing and PCA calcuations
###############################################################################
gene_ex_df = pd.read_csv(input_file, delimiter="\t", index_col=0)
#gene_ex_df_row_names = gene_ex_df.index
# print(gene_ex_df.head())  # Gives first 5 rows in dataset
# print(gene_ex_df.shape)  # Gives dimensions of dataset
# below scale data, but it is not needed here since data is already scaled
#scaled_gene_ex_df = preprocessing.scale(gene_ex_df.T)
scaled_gene_ex_df = gene_ex_df.T

pca = decomposition.PCA() # n_components=2; optional argument to put in here
pca.fit(scaled_gene_ex_df)
pca_data = pca.transform(scaled_gene_ex_df)
# run %matplotlib qt in console (doesn't work in script) to show graphs in a
# separate window
per_var = np.round(pca.explained_variance_ratio_ * 100, decimals=1)
labels = ['PC' + str(x) for x in range(1, len(per_var)+1)]

###############################################################################
# creates a PCA plot
###############################################################################
pca_df = pd.DataFrame(pca_data, index=gene_ex_df.columns, columns=labels)
plt.rc('figure', figsize=(25, 25))
plt.rc('ytick', labelsize=25)    # fontsize of the tick labels
plt.rc('xtick', labelsize=25)    # fontsize of the tick labels
plt.title('PCA plot', fontsize=40)
plt.xlabel('PC1 - {0}%'.format(per_var[0]), fontsize=30)
plt.ylabel('PC2 - {0}%'.format(per_var[1]), fontsize=30)

# column names corresponds to timepoints
gene_ex_df_col_names = gene_ex_df.columns
# creates a list of int values for colormap, where each value corresponds to a
# point on the colormap's gradient
values_for_colourmap = [colormap_dict[x] for x in gene_ex_df_col_names]
# numpy array of phases (light/dark), where 'L' corresponds to light ad 'D'
# corresponds to dark
# numpy array is needed for subsetting pca_df.PC1 and pca_df.PC2 below 
day_night_phase = np.array([x[0] for x in pca_df.PC1.index])
# '^' corresponds to dark and 'o' corresponds to light
markers = ["^","o"]
# converts values_for_colourmap into a numpy array, as the c argument requires
# sequences to be of this type, otherwise a AttributeError: 'list' object has no attribute 'shape'
# will be produced
for i, phase in enumerate(np.unique(day_night_phase)):
    #cseq refers to the sequence used to map PCA values onto the colourmap
    phase_specific_cseq = np.array(values_for_colourmap)[day_night_phase == phase]
    phase_specific_PC1 = pca_df.PC1[day_night_phase == phase]
    phase_specific_PC2 = pca_df.PC2[day_night_phase == phase]
    sc = plt.scatter(phase_specific_PC1, phase_specific_PC2, marker=markers[i],
                     c=phase_specific_cseq, vmin=1, vmax=24, cmap='plasma_r', s=400)
cbar = plt.colorbar(sc, ticks=[1, 5, 10, 15, 20, 24], orientation='horizontal')
cbar.set_label('Timepoint (h)', size=30)
legend = plt.legend(labels=('Dark', 'Light'), fontsize=30)
legend.get_frame().set_linewidth(3.0)
plt.savefig('PCA_plot_' + file_name + '.png')
plt.close()
