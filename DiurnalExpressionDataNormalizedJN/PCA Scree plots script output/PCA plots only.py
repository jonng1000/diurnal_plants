# -*- coding: utf-8 -*-
"""
Created on Mon May 21 14:04:34 2018

@author: weixiong001
"""
import numpy as np
import matplotlib.pyplot as plt
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


def PCA_analysis(input_file, file_name):
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
    plt.rc('figure', figsize=(25, 20))
    # remove ticks and tick labels from x and y axes
    plt.tick_params(
        axis='both',       # changes apply to both axes
        which='both',      # both major and minor ticks are affected
        bottom=False,      # ticks along the bottom edge are off
        left=False,        # ticks along the top edge are off
        labelbottom=False, # labels along the bottom edge are off
        labelleft=False)   # labels along the left edge are off
    plt.rc('axes', linewidth=10)
    plt.title('PCA plot', fontsize=100)
    plt.xlabel('PC1 - {0}%'.format(per_var[0]), fontsize=80)
    plt.ylabel('PC2 - {0}%'.format(per_var[1]), fontsize=80)
    # outputs a scatterplot and colours light and darkness samples differently
    plt.scatter(pca_df.PC1, pca_df.PC2, c=[use_colours[x[0]] for x in gene_ex_df.columns], s=3600)
    plt.savefig('PCA_plot_' + file_name + '.svg')
    plt.close()
    return None


for file in os.listdir('./'):
    if file.endswith('.txt'):
        input_file = file
        file_name = file.split('.')[0]
        PCA_analysis(input_file, file_name)
