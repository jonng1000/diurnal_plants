# -*- coding: utf-8 -*-
"""
Created on Wed Nov 21 12:59:47 2018

@author: weixiong
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from sklearn import decomposition
from sklearn import preprocessing

input_file = 'JTK.Cre_Mat.txt'
file_name = input_file.split('.')[1]

gene_ex_df = pd.read_csv(input_file, delimiter="\t", index_col=0)
# eliminates all genes with p-value >= 0.05
criteria = gene_ex_df['ADJ.P'] < 0.05
lower_than_p = gene_ex_df[criteria]
# selects for first 1000 rows
top_1000 = lower_than_p.head(1000)
cleaned = top_1000.drop(['BH.Q', 'ADJ.P', 'PER', 'AMP', 'LAG'], axis=1)

###############################################################################
# data preprocessing and PCA calcuations (modified from original one, dont
# copy this into other scrips)
###############################################################################
gene_ex_df_row_names = cleaned.index
# below scales data, but it is not needed here since RNA-seq data is already scalled
scaled_gene_ex_df = preprocessing.scale(cleaned, axis=1)
# if above is uncommented to scale data, remember to comment out below line
# NB: normally would transpose the dataset, cleaned.T, but this is not needed
# here, as I want PCA plot to be made using the genes, hence the original shape,
# genes x tme points (rows x columns) is correct, as PCA will use row names to
# make the plot
# hence scaled_gene_ex_df variable is redundant, but keeping it to reduce
# unnecessary changes to my downstream code
# scaled_gene_ex_df = cleaned

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
pca_df = pd.DataFrame(pca_data, index=cleaned.index, columns=labels)
plt.rc('figure', figsize=(25, 25))
plt.rc('ytick', labelsize=25)    # fontsize of the tick labels
plt.rc('xtick', labelsize=25)    # fontsize of the tick labels
plt.title('PCA plot', fontsize=40)
plt.xlabel('PC1 - {0}%'.format(per_var[0]), fontsize=30)
plt.ylabel('PC2 - {0}%'.format(per_var[1]), fontsize=30)

plt.scatter(pca_df.PC1, pca_df.PC2, c=top_1000['LAG'], vmin=0, vmax=23, cmap='plasma_r', s=400)
'''
for sample in pca_df.index:
    plt.annotate(sample, (pca_df.PC1.loc[sample], pca_df.PC2.loc[sample]),
                 fontsize=30)
'''
cbar = plt.colorbar(ticks=range(24), orientation='horizontal')
cbar.set_label('LAG', size=30)

plt.savefig('123_' + file_name + '.pdf')
plt.close()