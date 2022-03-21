# -*- coding: utf-8 -*-
"""
Created on Wed Nov 21 12:59:47 2018

@author: weixiong

This script takes all JTK output files, and does PCA on it, where each dot
represents one gene. Gene expression is scaled on a gene by gene basis.
Genes not showing significant rhythmic are coloured grey, the rest are
coloured.
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from sklearn import decomposition
from sklearn import preprocessing
import os

for file in os.listdir():
    if file.startswith('JTK'):
        # kind of redunant, but this makes it convienent as less modification
        # downstream is needed
        input_file = file
        # file_name is used for saving the plot
        file_name = input_file.split('.')[1]
        
        gene_ex_df = pd.read_csv(input_file, delimiter="\t", index_col=0)
        lag = gene_ex_df['LAG']
        p_values = gene_ex_df['ADJ.P']
        cleaned = gene_ex_df.drop(['BH.Q', 'ADJ.P', 'PER', 'AMP', 'LAG'], axis=1)
        
        ###############################################################################
        # data preprocessing and PCA calcuations (modified from original one, dont
        # copy this into other scripts)
        ###############################################################################
        gene_ex_df_row_names = cleaned.index
        # below scales data
        scaled_gene_ex_df = preprocessing.scale(cleaned, axis=1)
        # NB: normally would transpose the dataset, cleaned.T, but this is not needed
        # here, as I want PCA plot to be made using the genes, hence the original shape,
        # genes x tme points (rows x columns) is correct, as PCA will use row names to
        # make the plot
        
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
        pca_df = pd.DataFrame(pca_data, index=cleaned.index, columns=labels)
        pca_df['LAG'] = lag
        pca_df['ADJ.P'] = p_values
        
        plt.rc('figure', figsize=(25, 25))
        plt.rc('ytick', labelsize=25)    # fontsize of the tick labels
        plt.rc('xtick', labelsize=25)    # fontsize of the tick labels
        plt.title('PCA ' + file_name, fontsize=40)
        plt.xlabel('PC1 - {0}%'.format(per_var[0]), fontsize=30)
        plt.ylabel('PC2 - {0}%'.format(per_var[1]), fontsize=30)
        
        # find all genes with p-value < 0.05
        criteria = pca_df['ADJ.P'] < 0.05
        lower_than_p = pca_df[criteria]
        # find all genes with p-value >= 0.05
        higher_than_p = pca_df[~criteria]
        # genes whoose p-value for rhythmicity >= 0.05 are coloured grey
        plt.scatter(higher_than_p.PC1, higher_than_p.PC2, c='grey', s=400)
        # genes whoose p-value for rhythmicity < 0.05 are coloured according
        # to colourbar
        plt.scatter(lower_than_p.PC1, lower_than_p.PC2, c=lower_than_p['LAG'],
                    vmin=0, vmax=23, cmap='plasma_r', s=400)
        cbar = plt.colorbar(ticks=range(24), orientation='horizontal')
        cbar.set_label('LAG', size=30)
        
        plt.savefig('PCA_plot_z_' + file_name + '.pdf')
        plt.savefig('PCA_plot_z_' + file_name + '.png')
        plt.close()