# -*- coding: utf-8 -*-
"""
Created on Mon May 21 14:04:34 2018

@author: weixiong001
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd
from sklearn import decomposition
from sklearn import preprocessing  # don't need this if data is already scaled
import os
import csv
'''
This script finds all .txt files in the working directory, performs PCA
analysis on it, and generates PCA plots in .svg format,
Note: must ensure that all .txt files in the working directory are gene
expression files.
'''


def PCA_analysis(input_file, file_name, first_3):
    '''
    This function runs the codes to generate PCA plots.
    '''
    ###############################################################################
    # data preprocessing and PCA calcuations
    ###############################################################################
    gene_ex_df = pd.read_csv(input_file, delimiter="\t", index_col=0)
    gene_ex_df_row_names = gene_ex_df.index
    # print(gene_ex_df.head())  # Gives first 5 rows in dataset
    # print(gene_ex_df.shape)  # Gives dimensions of dataset
    # below scale data, but it is not needed here since RNA-seq data is already scalled
    # scaled_gene_ex_df = preprocessing.scale(gene_ex_df.T)
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
    use_colours = {"L": "green", "D": "cyan"}
    # remove ticks and tick labels from x and y axes
    plt.tick_params(
        axis='both',       # changes apply to both axes
        which='both',      # both major and minor ticks are affected
        bottom=False,      # ticks along the bottom edge are off
        left=False,        # ticks along the top edge are off
        labelbottom=False, # labels along the bottom edge are off
        labelleft=False)   # labels along the left edge are off
    plt.rc('axes', linewidth=1)
    plt.title(first_3, fontsize=9, y=0.96)
    # labelpad argument: space between label and axis
    plt.xlabel('PC1 - {0}%'.format(per_var[0]), labelpad=1.3, fontsize=7)
    plt.ylabel('PC2 - {0}%'.format(per_var[1]), labelpad=1.3, fontsize=7)
    # outputs a scatterplot and colours light and darkness samples differently
    plt.scatter(pca_df.PC1, pca_df.PC2, c=[use_colours[x[0]] for x in gene_ex_df.columns], s=24)
    return None

list_files = []
for file in os.listdir('./'):
    if file.endswith('.txt'):
        input_file = file
        file_name = file.split('.')[0]
        # first 3 letters of species name, will be used as PCA plot title
        first_3 = file_name[:3]
        list_files.append((input_file ,file_name, first_3)) 

plt.rc('figure', figsize=(10, 4))
# build the legend
green_patch = mpatches.Patch(color='green', label='Light')
cyan_patch = mpatches.Patch(color='cyan', label='Dark')
# set up for handles declaration
patches = [green_patch, cyan_patch]
# define and place the legend
legend = plt.figlegend(handles=patches, loc = (0.072, 0.01), fontsize=8)
# above function does no have arguments allowing legend's border width to be set,
# hence had to use the below code to do it
legend.get_frame().set_linewidth(0.0)
#legend.get_frame().set_edgecolor("black")  # optional
index = 0
for i in range(1, 9):
    plt.subplot(2, 4, i)
    PCA_analysis(*list_files[index])
    index += 1
plt.savefig('PCA_all_plots' + '.svg', bbox_inches='tight')
plt.close()
