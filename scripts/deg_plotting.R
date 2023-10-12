library(tidyverse)
library(ggrepel)
library(ComplexHeatmap)
library(optparse)

# Define the command-line options
option_list = list(
  make_option(c("--fname"), type="character", default=NULL, help="Input CSV file path"),
  make_option(c("--l2fc"), type="numeric", default=NULL, help="Log2 Fold Change Cutoff"),
  make_option(c("--padj"), type="numeric", default=NULL, help="Adjusted P-value Cutoff"),
  make_option(c("--outliers"), type="character", default=NULL, action="store", help="List of outliers")
)

opt_parser = OptionParser(option_list=option_list)

opt = parse_args(opt_parser)

fname <- opt$fname
l2fc_cutoff <- opt$l2fc
padj_cutoff <- opt$padj
outliers <- strsplit(opt$outliers, " ")[[1]]

new_fname = do.call("rbind", strsplit(fname, ".csv"))[,1]
new_fname = do.call("rbind", strsplit(new_fname, "/degs/"))[,2]
title = gsub("_", " ", new_fname)
print(title)

result = read.csv(fname)

######## Volcano plot
# add diffexpressed column
result <- result %>%
  mutate(diffexpressed = ifelse(
    !is.na(padj) & padj < padj_cutoff & !is.na(log2FoldChange),
    ifelse(log2FoldChange > l2fc_cutoff, "UP", ifelse(-log2FoldChange > l2fc_cutoff, "DOWN", "NO")),
    "NO"
  ))

print(table(result$diffexpressed))

result$delabel = NA
result$delabel[result$diffexpressed != "NO"] = result$gene_name[result$diffexpressed != "NO"]

result$delabel_gene_id = NA
result$delabel_gene_id[result$diffexpressed != "NO"] = result$gene_id[result$diffexpressed != "NO"]

volcano = ggplot(result, aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed, label=delabel)) +
                geom_point() + 
                theme_minimal() +
                geom_text_repel() +
                scale_color_manual(values=c("blue", "black", "red")) +
                geom_vline(xintercept=c(-l2fc_cutoff, l2fc_cutoff), col="darkgrey") +
                geom_hline(yintercept=-log10(padj_cutoff), col="darkgrey") + ggtitle(title)
# Save the plot
ggsave(paste0("../plots/",new_fname,"_volcano.png"), plot = volcano, width = 8, height = 8, dpi = 300)

# Assuming 'meta' is your metadata dataframe
meta = read.csv("../ref/mirna_practice_metadata.csv")
meta = meta[!(meta$sampleID %in% outliers),]
meta$sex <- factor(meta$sex)
meta$timepoint <- factor(meta$timepoint)
meta$technician <- factor(meta$technician)

# Create heatmap annotations
sex_palette = c("Male" = "blue", "Female" = "pink")
tp_palette = c('PND_14' = '#ec80eb', 'PNM_02' = '#800080')
tech_palette = c('NM' = '#72c2a8', 'SB' = '#ecaa2e', 'RM' = '#de731c')

column_ha = HeatmapAnnotation(sex = setNames(meta$sex, meta$sex), 
                               timepoint = setNames(meta$timepoint, meta$timepoint),
                               technician = setNames(meta$technician, meta$technician),
                               col = list(sex = sex_palette,timepoint = tp_palette, technician = tech_palette))

# Get scaled CPM matrix to plot, filtered by DEGs
cpm_2 = read.table("../counts/cpm_over2_matrix.tsv",sep="\t",header=TRUE, row.names=1)
cpm_2 = as.matrix(cpm_2[rownames(cpm_2) %in% na.omit(result$delabel_gene_id),])
cpm_2 = cpm_2[match(na.omit(result$delabel_gene_id),rownames(cpm_2)),]
rownames(cpm_2) = na.omit(result$delabel)
cpm_scaled = t(scale(t(cpm_2)))

cpm_scaled = cpm_scaled[,!(colnames(cpm_scaled) %in% outliers)]

# Create the heatmap
heatmap = Heatmap(cpm_scaled, 
                  top_annotation = column_ha,
                  name = "row-scaled    \nCPM", 
                  cluster_rows = TRUE,
                  cluster_columns = TRUE,
                  show_column_names = TRUE,
                  show_row_names = TRUE)

# Save the plot to a PNG file
png(paste0("../plots/",new_fname,"_heatmap.png"), width = 8, height = 16, units = "in", res = 300)
heatmap
dev.off()

