import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

meta = pd.read_csv("../ref/mirna_metadata.csv")

cpm_2 = pd.read_csv("../counts/cpm_over2_matrix.tsv", sep="\t", index_col=0)

# run PCA
input_df = cpm_2.transpose() # transpose the matrix so samples are in rows
pca = PCA(cpm_2.shape[1]) # compute the same # of PCs as samples
pca.fit(input_df)

############ Percent variance explained bar plot
# Make df to plot the percent variance explained
pct_explained_df = pd.DataFrame(data=pca.explained_variance_ratio_,
                                columns=['prop_variance'])

# compute percent variance from proportion variance
pct_explained_df['pct_variance'] = pct_explained_df.prop_variance*100
pct_explained_df['PC'] = list(range(1, cpm_2.shape[1]+1))
print(pct_explained_df)

# Set the figure size in inches
plt.figure(figsize=(4, 5))

# Use Seaborn library to make categorical bar plot
ax = sns.catplot(data=pct_explained_df, x='PC', y='pct_variance', kind='bar')

# Add nice x and y axis labels
ax.set(ylabel='Percent of variance explained', xlabel='Principal Component')

# Save the plot 
plt.savefig('../plots/pca_variance_explained.png', dpi=300)

############ actual PCA plot

# first calculate the PCA representation of the data
pca_data = pca.transform(input_df) 
cols = ['PC{} ({:.1f}%)'.format(pc, var) for pc, var in zip(pct_explained_df.PC, pct_explained_df.pct_variance)]
print(cols)
pca_df = pd.DataFrame(data=pca_data, columns=cols)
pca_df.set_index(input_df.index, inplace=True)

# Perform a merge based on matching indices
pca_df = pca_df.merge(meta, left_index=True, right_on='sampleID')

# Set the figure size in inches
plt.figure(figsize=(5, 4))

ax = sns.scatterplot(data=pca_df, x=cols[0], y=cols[1], 
                     s=150, # size of point
                     hue='stage', style='genotype')
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize='medium')

# Save the plot 
plt.savefig('../plots/pca.png', dpi=300, bbox_inches='tight')
