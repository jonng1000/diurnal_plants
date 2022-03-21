# -*- coding: utf-8 -*-
"""
Created on Wed Nov 21 12:59:47 2018

@author: weixiong

This script takes all JTK output files, and does PCA on it, where each dot
represents the sample's time point. Gene expression is scaled on a time point by timepoint basis.
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from sklearn import decomposition
from sklearn import preprocessing
import os, csv

# Creates master dictionary, cmap_master_dict, which is a dictonary of
# dictionaries. Species (sample) names are the keys for the top level dictionary, while
# timepoints from that species are keys for the bottom level dictionary, and
# colourmap scores corresponding to the timepoints are their values.
# Uses the colormap_values_PCA file which provides the colour mapping for the
# creation of the master dictionary

# This is done because the timepoints mean different things in each sample, and
# a common colourbar will be used to represent all of them. Hence the
# same timepoint between samples may correspond to different colourmap scores.

# Adjusted figure size so that there isn't too much empty space in plots
# Fig size 10 chosen as smaller means dots go out of plot and bigger means too
# much empty space
# Removed tick labels from x and y axes
cmap_master_dict = {}
with open('colormap_values_PCA.csv', 'r', newline='') as csvfile:
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


for file in os.listdir():
    if file.startswith('JTK'):
        # kind of redunant, but this makes it convienent as less modification
        # downstream is needed
        input_file = file
        # file_name is used for saving the plot
        file_name = input_file.split('.')[1]
        gene_ex_df = pd.read_csv(input_file, delimiter="\t", index_col=0)
        cleaned = gene_ex_df.drop(['BH.Q', 'ADJ.P', 'PER', 'AMP', 'LAG'], axis=1)
        
        ###############################################################################
        # data preprocessing and PCA calcuations (modified from original one, dont
        # copy this into other scripts)
        ###############################################################################
        gene_ex_df_row_names = cleaned.index
        # below scales data
        scaled_gene_ex_df = preprocessing.scale(cleaned.T, axis=1)
        # NB: transposed the dataset, cleaned.T, as it is needed here, as I 
        # want PCA plot to be made using the samples' time points and not 
        # genes, hence the original shape, genes x tme points (rows x columns) 
        # is incorrect, as PCA will use row names to make the plot
        pca = decomposition.PCA() # n_components=2; optional argument to put in here
        pca.fit(scaled_gene_ex_df)
        pca_data = pca.transform(scaled_gene_ex_df)
        # run %matplotlib qt in console (doesn't work in script) to show graphs in a
        # separate window
        per_var = np.round(pca.explained_variance_ratio_ * 100, decimals=1)
        labels = ['PC' + str(x) for x in range(1, len(per_var)+1)]
        
        ###############################################################################
        # creates a PCA plot (modified from original one, dont
        # copy this into other scripts)
        ###############################################################################
        pca_df = pd.DataFrame(pca_data, index=cleaned.columns, columns=labels)
        plt.rc('figure', figsize=(10, 10))
        # remove ticks and tick labels from x and y axes
        plt.tick_params(
            axis='both',        # changes apply to both axes
            which='both',       # both major and minor ticks are affected
            bottom=False,       # ticks along the bottom edge are off
            left=False,         # ticks along the top edge are off
            labelbottom=False,  # labels along the bottom edge are off
            labelleft=False)    # labels along the left edge are off
        plt.title('PCA ' + file_name, fontsize=40)
        plt.xlabel('PC1 - {0}%'.format(per_var[0]), fontsize=30)
        plt.ylabel('PC2 - {0}%'.format(per_var[1]), fontsize=30)

        # time points are coloured according to colourbar
        values_for_colourmap = [cmap_master_dict[file_name][x] for x in pca_df.index]
        plt.scatter(pca_df.PC1, pca_df.PC2, c=values_for_colourmap,
                    vmin=0, vmax=23, cmap='plasma_r', s=1200)
        cbar = plt.colorbar(ticks=range(24), orientation='horizontal')
        cbar.set_label('LAG', size=30)
        
        plt.savefig('PCA_plot_z_' + file_name + '.pdf')
        plt.savefig('PCA_plot_z_' + file_name + '.png')
        plt.close()