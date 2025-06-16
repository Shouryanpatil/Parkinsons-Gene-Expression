#BiocManager::install("DEGreport")
library(DESeq2)
library(tximport)
library(readr)
library(DEGreport)


## Run analysis
dds <- DESeq(dds)

## Check the size factors
sizeFactors(dds)

## Total number of raw counts per sample
colSums(counts(dds))

## Total number of normalized counts per sample
colSums(counts(dds, normalized=T))

## Plot dispersion estimates
plotDispEsts(dds)

# Extract Results (Hypothesis Testing)
res <- results(dds)

## Define contrasts for MOV10 overexpression
contrast_oe <- c("condition", "Disease", "Control")

# Extract Results for That Contrast
res_tableOE <- results(dds, contrast = contrast_oe, alpha = 0.05)

# Check what type of object is returned
class(res_tableOE)

# View Results Table
res_tableOE %>%
  data.frame() %>%
  View()

# Check Column Descriptions (Optional but informative)
mcols(res_tableOE, use.names = TRUE)

# Filter genes by zero expression
res_tableOE[which(res_tableOE$baseMean == 0),] %>% 
  data.frame() %>% 
  View()

# Filter genes that have an extreme outlier
res_tableOE[which(is.na(res_tableOE$pvalue) &
                    is.na(res_tableOE$padj) &
                    res_tableOE$baseMean > 0), ] %>%
  data.frame() %>%
  View()

# Filter Genes Below the Low Mean Threshold (Independent Filtering)
res_tableOE[which(!is.na(res_tableOE$pvalue) &
                    is.na(res_tableOE$padj) &
                    res_tableOE$baseMean > 0), ] %>%
  data.frame() %>%
  View()

# Save the unshrunken results (for comparison later)
res_tableOE_unshrunken <- res_tableOE

# Check available coefficients for lfcShrink()
resultsNames(dds)

BiocManager::install("apeglm")
library(apeglm)

# Apply fold change shrinkage
res_tableOE <- lfcShrink(dds, coef="condition_Disease_vs_Control", type="apeglm")

# Generate MA plot using unshrunken results
plotMA(res_tableOE_unshrunken, ylim=c(-2,2))

# Generate MA plot using shrunken results
plotMA(res_tableOE, ylim=c(-2,2))

# Summarize the results table
summary(res_tableOE, alpha = 0.05)

# Set significance threshold
padj.cutoff <- 0.05

# Convert results to tibble
res_tableOE_tb <- res_tableOE %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble()

# Subset to significant genes (adjusted p-value < 0.05)
sigOE <- res_tableOE_tb %>%
  dplyr::filter(padj < padj.cutoff)

# View significant genes
sigOE

# Exercise 

# Count of DE genes
nrow(sigOE)

# Number of upregulated genes
sigOE %>% filter(log2FoldChange > 0) %>% nrow()

# Number of downregulated genes
sigOE %>% filter(log2FoldChange < 0) %>% nrow()







# Likelihood ratio test
dds_lrt <- DESeq(dds, test="LRT", reduced = ~ 1)

# Extract results for LRT
res_LRT <- results(dds_lrt)

# View results for LRT
res_LRT

# Create a tibble for LRT results
res_LRT_tb <- res_LRT %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

# Subset to return genes with padj < 0.05
sigLRT_genes <- res_LRT_tb %>% 
  dplyr::filter(padj < padj.cutoff)

# Get number of significant genes
nrow(sigLRT_genes)

# Compare to numbers we had from Wald test
nrow(sigOE)