# Visualizion 
library(dplyr)
library(tibble)
library(tidyverse)
library(DESeq2)
library(RColorBrewer)
library(pheatmap)
library(ggplot2)

# Create metadata tibble
metadata_tb <- metadata %>%
  rownames_to_column(var = "samplename") %>%
  as_tibble()

# Metadata with only Condition
metadata_c <- metadata[, "condition", drop = FALSE]

# Prepare normalized counts for visualization
normalized_counts_tb <- normalized_counts %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble()

head(rownames(dds))
# 
# # Choose a significant gene from your sigOE results
# # View the top significant genes (adjust or filter as you like)
# sigOE %>% arrange(padj) %>% View()
# 
# top_gene <- sigOE %>% arrange(padj) %>% slice(1) %>% pull(gene)
# top_gene
# 
# class(sigOE$padj)
# class(sigOE$gene)

# Heatmap
### Extract normalized expression for significant genes from the OE and control samples (2:4 and 7:9)
norm_OEsig <- normalized_counts_tb[,c(1:5,6:9)] %>% 
  dplyr::filter(gene %in% sigOE$gene)  

### Set a color palette
heat_colors <- brewer.pal(6, "YlOrRd")

### Run pheatmap using the metadata data frame for the annotation
pheatmap(norm_OEsig[2:9], 
         color = heat_colors, 
         cluster_rows = T, 
         show_rownames = F,
         annotation = metadata_c, 
         border_color = NA, 
         fontsize = 10, 
         scale = "row", 
         fontsize_row = 10, 
         height = 20)


# Heatmap of Top DEGs

library(dplyr)
library(RColorBrewer)
library(pheatmap)

### Define number of top genes to plot
top_n <- 50  # Change this to 20, 100, etc., as needed

### Select top N significant DEGs based on adjusted p-value
top_sigOE <- sigOE %>%
  arrange(padj) %>%
  head(top_n)

### Filter normalized counts for just those top genes
norm_top_OEsig <- normalized_counts_tb %>%
  dplyr::filter(gene %in% top_sigOE$gene)

# Convert tibble to data frame
norm_top_OEsig_df <- as.data.frame(norm_top_OEsig)

# Set gene names as rownames
rownames(norm_top_OEsig_df) <- norm_top_OEsig_df$gene

# Remove the gene column
norm_top_OEsig_df$gene <- NULL

# Convert to matrix
norm_top_OEsig_matrix <- as.matrix(norm_top_OEsig_df)

# Plot heatmap
pheatmap(norm_top_OEsig_matrix, 
         color = heat_colors, 
         cluster_rows = TRUE, 
         show_rownames = TRUE,
         annotation_col = metadata_c, 
         border_color = NA, 
         fontsize = 10, 
         scale = "row", 
         fontsize_row = 8, 
         height = 10)


## Obtain logical vector where TRUE values denote padj values < 0.05 and fold change > 1.5 in either direction
res_tableOE_tb <- res_tableOE_tb %>% 
  dplyr::mutate(threshold_OE = padj < 0.05 & abs(log2FoldChange) >= 0.58)

## Volcano plot
ggplot(res_tableOE_tb) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = threshold_OE)) +
  ggtitle("Expression") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  scale_x_continuous(limits = c(-10, 10)) +   # Customize x-axis
  scale_y_continuous(limits = c(0, 250)) +     # Customize y-axis
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))) 

