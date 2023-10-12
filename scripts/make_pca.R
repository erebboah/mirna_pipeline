library(tidyverse)

meta = read.csv("../ref/mirna_practice_metadata.csv")

cpm_2 = read.table("../counts/cpm_over2_matrix.tsv",sep="\t",header=TRUE, row.names=1)

# Perform PCA using prcomp
pca = prcomp(t(cpm_2))

# Extract percentage of variance explained from summary
pct_explained_df = data.frame(
  PC = seq_len(ncol(cpm_2)),
  pct_variance = summary(pca)$importance[2,]
)

# Make bar plot for percent variance explained
bar_plot = ggplot(pct_explained_df, aes(x = factor(PC), y = pct_variance*100)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  labs(y = "Percent of variance explained", x = "Principal Component") +
  theme_minimal()

ggsave("../plots/pca_variance_explained.png", plot = bar_plot, width = 4, height = 5, dpi = 300)

############ Actual PCA plot
# Create PCA data frame
pca_df = as.data.frame(pca$x)
pca_df$sampleID = colnames(cpm_2)

# Merge with metadata
pca_df = merge(pca_df, meta, by = "sampleID")

x_axis_labels = sprintf("PC%d (%.1f%%)", seq_len(ncol(cpm_2)), 100 * pct_explained_df$pct_variance)

color_palette = c('PND_14' = '#ec80eb', 'PNM_02' = '#800080')

# Scatter plot
scatter_plot = ggplot(pca_df, aes(x = PC1, y = PC2, color = timepoint, shape = technician)) +
  geom_point(size = 3)+
  scale_color_manual(values = color_palette) +
  labs(x = x_axis_labels[1], y = x_axis_labels[2]) +
  theme_minimal() 

# Save the plot
ggsave("../plots/pca.png", plot = scatter_plot, width = 5, height = 4, dpi = 300)