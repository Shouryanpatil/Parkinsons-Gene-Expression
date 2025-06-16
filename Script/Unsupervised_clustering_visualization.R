## QC

#Variance Stabilization 
rld <- rlog(dds, blind=TRUE)

# PCA Plot (Sample Clustering Visualization)
plotPCA(rld, intgroup="condition")

# Explore More PCs (Optional but recommended if PCA unclear)
pca <- prcomp(t(assay(rld)))
df <- cbind(metadata, pca$x)
ggplot(df) + geom_point(aes(x=PC3, y=PC4, color=condition))

## Sample-to-Sample Distance Heatmap (Hierarchical Clustering)
library(pheatmap)

### Extract the rlog matrix from the object
rld_mat <- assay(rld)  

### Compute pairwise correlation values
rld_cor <- cor(rld_mat)    ## cor() is a base R function


### Plot heatmap using the correlation matrix and the metadata object
pheatmap(rld_cor, annotation = metadata)
