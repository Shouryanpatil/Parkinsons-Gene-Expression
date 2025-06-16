if (!requireNamespace("clusterProfiler")) install.packages("clusterProfiler")
if (!requireNamespace("org.Hs.eg.db")) install.packages("org.Hs.eg.db")
if (!requireNamespace("enrichplot")) BiocManager::install("enrichplot")

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

## Prapare Gene list
# Filter for significantly DEGs (padj < 0.05)
sig_genes <- sigOE_annotated[sigOE_annotated$padj < 0.05, ]

# Remove NAs from ENTREZID
sig_genes <- sig_genes[!is.na(sig_genes$ENTREZID), ]

# Get vector of ENTREZ IDs
entrez_ids <- unique(sig_genes$ENTREZID)

## GO Biological Process Enrichment
go_bp <- enrichGO(gene = entrez_ids,
                  OrgDb = org.Hs.eg.db,
                  keyType = "ENTREZID",
                  ont = "BP",       # BP: Biological Process
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  readable = TRUE)

# View top results
head(go_bp)

## Output results from GO analysis to a table
cluster_summary <- data.frame(go_bp)

write.csv(cluster_summary, "clusterProfiler.csv")

## Dotplot 
dotplot(go_bp, showCategory = 50, title = "GO: Biological Process") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_y_discrete(labels = label_wrap_gen(width = 100))


## KEGG Pathway Enrichment
kegg <- enrichKEGG(gene = entrez_ids,
                   organism = "hsa", # hsa = human
                   pvalueCutoff = 0.05)

# Make readable
kegg <- setReadable(kegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
head(kegg)

## Output results from Pathway Enrichment to a table
Ppathway_Enrichment <- data.frame(kegg)

write.csv(Ppathway_Enrichment, "pathway.csv")

# Barplot
barplot(kegg, showCategory = 50, title = "KEGG Pathways")+
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_y_discrete(labels = label_wrap_gen(width = 100))
