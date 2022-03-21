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
import os
import csv
'''
This script takes one .txt file and performs PCA
analysis on it. It generatres a Scree and PCA plot in .png format, and a .txt
file containing the top 10 genes with the biggest loading scores for PC1
Note: must ensure that all .txt files in the working directory are gene
expression files, except for those containing top 10 genes, as the script will
ignore those
'''

# uncomment when this script is run using only one input file
input_file = 'Smo_Mat.txt'
file_name = input_file.split('.')[0]
output_file = 'top10_' + file_name + '.txt'

###############################################################################
# data preprocessing and PCA calcuations
###############################################################################
gene_ex_df = pd.read_csv(input_file, delimiter="\t", index_col=0)
gene_ex_df_row_names = gene_ex_df.index
# print(gene_ex_df.head())  # Gives first 5 rows in dataset
# print(gene_ex_df.shape)  # Gives dimensions of dataset
# below scales data, but it is not needed here since RNA-seq data is already scalled
# scaled_gene_ex_df = preprocessing.scale(gene_ex_df.T)
# if above is uncommented to scale data, remember to comment out below line
scaled_gene_ex_df = gene_ex_df.T

pca = decomposition.PCA() # n_components=2; optional argument to put in here
pca.fit(scaled_gene_ex_df)
pca_data = pca.transform(scaled_gene_ex_df)
# run %matplotlib qt in console (doesn't work in script) to show graphs in a
# separate window
per_var = np.round(pca.explained_variance_ratio_ * 100, decimals=1)
labels = ['PC' + str(x) for x in range(1, len(per_var)+1)]

###############################################################################
# creates a Scree plot
###############################################################################
plt.rc('figure', figsize=(25, 20))
plt.rc('ytick', labelsize=25)    # fontsize of the tick labels
plt.rc('xtick', labelsize=10)    # fontsize of the tick labels
plt.bar(x=range(1, len(per_var)+1), height=per_var, tick_label=labels)
plt.ylabel('Percentage of Explained Variance', fontsize=30)
plt.xlabel('Principal Component', fontsize=30)
plt.title('Scree plot', fontsize=40)
# uncomment below to show Scree plot, note that Spyder by default runs iPython,
# which automatically this command even if its not entered
# plt.show()
plt.savefig('Scree_plot_' + file_name + '.png')
# clears figure from memory and closes window if it is displayed in a separate
# window
plt.close()

###############################################################################
# creates a PCA plot
###############################################################################
pca_df = pd.DataFrame(pca_data, index=gene_ex_df.columns, columns=labels)
use_colours = {"L": "green", "D": "cyan"}
plt.rc('xtick', labelsize=25)    # fontsize of the tick labels
plt.title('PCA plot', fontsize=40)
plt.xlabel('PC1 - {0}%'.format(per_var[0]), fontsize=30)
plt.ylabel('PC2 - {0}%'.format(per_var[1]), fontsize=30)
# outputs a scatterplot and colours light and darkness samples differently
plt.scatter(pca_df.PC1, pca_df.PC2, c=[use_colours[x[0]] for x in gene_ex_df.columns], s=400)
# see comments on plt.show() above code section to understand why it is
# commented out here
# plt.show()
plt.rc('font', size=20)  # controls default text sizes
for sample in pca_df.index:
    plt.annotate(sample, (pca_df.PC1.loc[sample], pca_df.PC2.loc[sample]))
plt.savefig('PCA_plot_' + file_name + '.png')
plt.close()

###############################################################################
# gets the name of the top 10 genes that contribute most to PC1 (largest loading
# scores) and writes it to a text fie.
###############################################################################
loading_scores = pd.Series(pca.components_[0], index=gene_ex_df_row_names)
sorted_loading_scores = loading_scores.abs().sort_values(ascending=False)
top_10_genes = sorted_loading_scores[0:10].index.values
with open(output_file, 'w', newline='') as f:
    writer = csv.writer(f, delimiter='\t')
    writer.writerow(['Top 10 genes with largest loading scores (absolute value) contributing to PC1'])
    for gene, value in loading_scores[top_10_genes].iteritems():
        writer.writerow([gene, value])
        