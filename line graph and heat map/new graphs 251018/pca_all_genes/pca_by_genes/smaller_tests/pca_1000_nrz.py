# -*- coding: utf-8 -*-
"""
Created on Wed Nov 21 12:59:47 2018

@author: weixiong

sp_file is used for to get LAG values as input_file doesnt have them, used for
colouring
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from sklearn import decomposition
from sklearn import preprocessing

input_file = 'Cre_Mat.txt'
sp_file = 'JTK.Cre_Mat.txt'
file_name = input_file.split('.')[0]

gene_ex_df = pd.read_csv(input_file, delimiter="\t", index_col=0)
sp_df = pd.read_csv(sp_file, delimiter="\t", index_col=0)
# selects for first 1000 rows
top_1000 = gene_ex_df.head(1000)

###############################################################################
# data preprocessing and PCA calcuations (modified from original one, dont
# copy this into other scrips)
###############################################################################
gene_ex_df_row_names = top_1000.index

# below scales data, but it is not needed here since RNA-seq data is already scalled
scaled_gene_ex_df = preprocessing.scale(top_1000, axis=1)
# if above is uncommented to scale data, remember to comment out below line
# NB: normally would transpose the dataset, cleaned.T, but this is not needed
# here, as I want PCA plot to be made using the genes, hence the original shape,
# genes x tme points (rows x columns) is correct, as PCA will use row names to
# make the plot
# hence scaled_gene_ex_df variable is redundant, but keeping it to reduce
# unnecessary changes to my downstream code
#scaled_gene_ex_df = top_1000

pca = decomposition.PCA() # n_components=2; optional argument to put in here
pca.fit(scaled_gene_ex_df)
pca_data = pca.transform(scaled_gene_ex_df)
# run %matplotlib qt in console (doesn't work in script) to show graphs in a
# separate window
per_var = np.round(pca.explained_variance_ratio_ * 100, decimals=1)
labels = ['PC' + str(x) for x in range(1, len(per_var)+1)]

###############################################################################
# creates a PCA plot (modified from original one, dont
# copy this into other scrips)
###############################################################################
pca_df = pd.DataFrame(pca_data, index=top_1000.index, columns=labels)
plt.rc('figure', figsize=(25, 25))
plt.rc('ytick', labelsize=25)    # fontsize of the tick labels
plt.rc('xtick', labelsize=25)    # fontsize of the tick labels
plt.title('PCA plot', fontsize=40)
plt.xlabel('PC1 - {0}%'.format(per_var[0]), fontsize=30)
plt.ylabel('PC2 - {0}%'.format(per_var[1]), fontsize=30)

# requires sp_file
# this maps gene IDs from pca_df onto their respective LAG values from sp_file
lag_1000 = sp_df['LAG'].loc[pca_df.index]
p_values = sp_df['ADJ.P'].loc[pca_df.index]
pca_df['LAG'] = lag_1000
pca_df['ADJ.P'] = p_values

# finds all genes with p-value < 0.05
criteria = pca_df['ADJ.P'] < 0.05
lower_than_p = pca_df[criteria]
# find all genes with p-value >= 0.05
higher_than_p = pca_df[~criteria]

plt.scatter(higher_than_p.PC1, higher_than_p.PC2, c='grey', s=400)
plt.scatter(lower_than_p.PC1, lower_than_p.PC2, c=lower_than_p['LAG'],
            vmin=0, vmax=23, cmap='plasma_r', s=400)
cbar = plt.colorbar(ticks=range(24), orientation='horizontal')
cbar.set_label('LAG', size=30)

plt.savefig('test_PCA_nr_z_' + file_name + '.pdf')
plt.close()