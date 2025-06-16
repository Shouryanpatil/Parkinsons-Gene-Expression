# Load the required libraries
library(org.Hs.eg.db)
library(AnnotationDbi)

# Extract unique gene symbols from your data
gene_symbols <- unique(sigOE$gene)

# Annotate gene symbols using org.Hs.eg.db
annot <- AnnotationDbi::select(org.Hs.eg.db,
                               keys = gene_symbols,
                               columns = c("ENSEMBL", "ENTREZID", "GENENAME"),
                               keytype = "SYMBOL")

# View the annotation table
head(annot)

# Merge annotation with sigOE on gene symbol
sigOE_annotated <- merge(sigOE, annot, by.x = "gene", by.y = "SYMBOL", all.x = TRUE)

# View result
head(sigOE_annotated)