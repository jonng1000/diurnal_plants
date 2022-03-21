# -*- coding: utf-8 -*-
"""
Created on Mon May 21 14:04:34 2018

@author: weixiong001
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from sklearn import decomposition
import os
import csv
'''
This script finds a maximum of 8 .txt files in the working directory, performs PCA
analysis on it, and generates PCA plots in .svg format. It also assigns different shapes
to each data point in the graph, based on whether the data point is measured
under day or night conditions.
Note: must ensure that all .txt files in the working directory are gene
expression files.
'''
# Creates master dictionary, cmap_master_dict, which is a dictonary of
# dictionaries. Species (sample) names are the keys for the top level dictionary, while
# timepoints from that species are keys for the bottom level dictionary, and
# colourmap scores corresponding to the timepoints are their values.

# This is done because the timepoints mean different things in each sample, and
# a common colourbar will be used to represent all of them. Hence the
# same timepoint between samples may correspond to different colourmap scores.
cmap_master_dict = {}
with open('colormap_values_PCA draft.csv', 'r', newline='') as csvfile:
    reader = csv.reader(csvfile, delimiter=',')
    next(reader)
    for row in reader:
        species = row[0]
        time = row[1]
        cmap_score = row[2]
        if species not in cmap_master_dict:
            cmap_master_dict[species] = {}
        cmap_master_dict[species][time] = {}
        cmap_master_dict[species][time] = cmap_score


def PCA_analysis(input_file, file_name, first_3):
    '''
    This function runs the codes to generate PCA plots.
    '''
    ###############################################################################
    # data preprocessing and PCA calcuations
    ###############################################################################
    gene_ex_df = pd.read_csv(input_file, delimiter="\t", index_col=0)
    # print(gene_ex_df.head())  # Gives first 5 rows in dataset
    # print(gene_ex_df.shape)  # Gives dimensions of dataset
    # below scale data, but it is not needed here since RNA-seq data is already scalled
    # scaled_gene_ex_df = preprocessing.scale(gene_ex_df.T)
    scaled_gene_ex_df = gene_ex_df.T

    pca = decomposition.PCA()  # n_components=2; optional argument to put in here
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
    # remove ticks and tick labels from x and y axes
    plt.tick_params(
        axis='both',        # changes apply to both axes
        which='both',       # both major and minor ticks are affected
        bottom=False,       # ticks along the bottom edge are off
        left=False,         # ticks along the top edge are off
        labelbottom=False,  # labels along the bottom edge are off
        labelleft=False)    # labels along the left edge are off
    plt.rc('axes', linewidth=1)
    plt.title(first_3, fontsize=9, y=0.96)
    # labelpad argument: space between label and axis
    plt.xlabel('PC1 - {0}%'.format(per_var[0]), labelpad=1.3, fontsize=7)
    plt.ylabel('PC2 - {0}%'.format(per_var[1]), labelpad=1.3, fontsize=7)
    # column names corresponds to timepoints
    gene_ex_df_col_names = gene_ex_df.columns

    # creates a list of int values for colormap, where each value corresponds to a
    # point on the colormap's gradient
    values_for_colourmap = [cmap_master_dict[first_3][x] for x in gene_ex_df_col_names]
    # numpy array of phases (light/dark), where 'L' corresponds to light ad 'D'
    # corresponds to dark
    # numpy array is needed for subsetting pca_df.PC1 and pca_df.PC2 below 
    day_night_phase = np.array([x[0] for x in pca_df.PC1.index])
    # '^' corresponds to dark and 'o' corresponds to light
    markers = ["^","o"]
    for i, phase in enumerate(np.unique(day_night_phase)):
        # cseq refers to the sequence used to map PCA values onto the colourmap
        # it needs to be a numpy array, as the c argument requires it,
        # otherwise a AttributeError: 'list' object has no attribute 'shape'
        # will be produced
        phase_specific_cseq = np.array(values_for_colourmap)[day_night_phase == phase]
        phase_specific_PC1 = pca_df.PC1[day_night_phase == phase]
        phase_specific_PC2 = pca_df.PC2[day_night_phase == phase]
        sc = plt.scatter(phase_specific_PC1, phase_specific_PC2, marker=markers[i], edgecolor='black', linewidth=0.3,
                         c=phase_specific_cseq, vmin=1, vmax=24, cmap='plasma_r', s=24)
    return None

# This section runs through all files in the current working directory, and
# creates a list of input files for PCA plotting
list_files = []
for file in os.listdir('./'):
    if file.endswith('.txt'):
        input_file = file
        file_name = file.split('.')[0]
        # first 3 letters of species name, will be used as PCA plot title
        first_3 = file_name[:3]
        list_files.append((input_file, file_name, first_3))

# This sections calls PCA_analysis(input_file, file_name, first_3), generates
# up to 8 PCA plots, and saves them in one .svg file
plt.rc('figure', figsize=(10, 4.4))
index = 0
fig, ax = plt.subplots(2, 4)
for i in range(1, 9):
    plt.subplot(2, 4, i)
    PCA_analysis(*list_files[index])
    index += 1
plt.subplots_adjust(left=0.397, bottom=0.3)
plt.subplots_adjust(hspace=0.35)
cbar = plt.colorbar(ax=ax, ticks=[1, 5, 10, 15, 20, 24], orientation='horizontal', shrink=0.3)
cbar.ax.tick_params(labelsize=9)
cbar.set_label('Timepoint (h)', size=9)
legend = plt.legend(labels=('Dark', 'Light'), loc=(-1, -0.6), fontsize=9)
legend.get_frame().set_linewidth(0.0)
plt.savefig('PCA_all_plots' + '.svg', bbox_inches='tight')
plt.close()
