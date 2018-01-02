
#Call libraries
library(monocle)
library (reshape2)
library (plyr)

#Import datasets (Count matrix; Gene annotation and Phenodata)
count_matrix <- read.table("count.txt", sep = "\t", header = TRUE, row.names = 1)
gene_ann <- read.table("annotation.txt", sep = "\t", header = TRUE, row.names = 1)
sample_sheet <- read.table("pheno.txt", sep = "\t", header = TRUE, row.names = 1)

#Create CellDataSet object
pd <- new("AnnotatedDataFrame", data = sample_sheet)
fd <- new("AnnotatedDataFrame", data = gene_ann)
ATLAS <- new("CellDataSet",
        exprs = as.matrix(count_matrix, "sparseMatrix"),
        lowerDetectionLimit = 0.5,
        phenoData = pd, featureData = fd,
        expressionFamily=negbinomial.size()) #Use negative binomial with fixed variance function

#Estimate size factors
ATLAS <- estimateSizeFactors(ATLAS)
ATLAS <- estimateDispersions(ATLAS)

ATLAS <- detectGenes(ATLAS, min_expr = 50)
print(head(fData(ATLAS)))
expressed_genes <- row.names(subset(fData(ATLAS),
    num_cells_expressed >= 4))

# Log-transform each value in the expression matrix.
L <- log(exprs(ATLAS[expressed_genes,]))

# Standardize each gene, so that they are all on the same scale,
# Then melt the data with plyr so we can plot it easily
melted_dens_df <- melt(Matrix::t(scale(Matrix::t(L))))

# Plot the distribution of the standardized gene expression values.
qplot(value, geom = "density", data = melted_dens_df) +
stat_function(fun = dnorm, size = 0.5, color = 'red') +
xlab("Standardized log(Counts)") +
ylab("Density")




LGR5_id <- row.names(subset(fData(ATLAS), gene_short_name == "LGR5"))
OLFM4_id <- row.names(subset(fData(ATLAS), gene_short_name == "OLFM4"))

cth <- newCellTypeHierarchy()
cth <- addCellType(cth, "Fetal Tissue",
       classify_func = function(x) { x[LGR5_id,] >= 800 })
cth <- addCellType(cth, "Adult Tissue",
       classify_func = function(x)
       { x[LGR5_id,] < 8000 & x[OLFM4_id,] > 1000 })

ATLAS <- classifyCells(ATLAS, cth, 0.1)
table(pData(ATLAS)$CellType)

pie <- ggplot(pData(ATLAS),
aes(x = factor(1), fill = factor(CellType))) + geom_bar(width = 1)
pie + coord_polar(theta = "y") +
theme(axis.title.x = element_blank(), axis.title.y = element_blank())


marker_diff <- markerDiffTable(ATLAS[expressed_genes,],
            cth,
            residualModelFormulaStr = "~type + num_genes_expressed",
            cores = 1)

candidate_clustering_genes <- row.names(subset(marker_diff, qval < 0.5))
marker_spec <- calculateMarkerSpecificity(ATLAS[candidate_clustering_genes,], cth)
head(selectTopMarkers(marker_spec, 3))


semisup_clustering_genes <- unique(selectTopMarkers(marker_spec, 500)$gene_id)
ATLAS <- setOrderingFilter(ATLAS, semisup_clustering_genes)
plot_ordering_genes(ATLAS)

ATLAS <- reduceDimension(ATLAS, max_components = 2, num_dim = 3,
  norm_method = 'log',
  reduction_method = 'tSNE',
  residualModelFormulaStr = "~type + num_genes_expressed",
  verbose = T)
ATLAS <- clusterCells(ATLAS, num_clusters = 7)
plot_cell_clusters(ATLAS, 1, 2, color = "group")




###########################################################################################

semisup_clustering_genes <- unique(selectTopMarkers(marker_spec, 500)$gene_id)
HSMM <- setOrderingFilter(HSMM, semisup_clustering_genes)
plot_ordering_genes(HSMM)













disp_table <- dispersionTable(ATLAS)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 200)
ATLAS <- setOrderingFilter(ATLAS, unsup_clustering_genes$gene_id)
plot_ordering_genes(ATLAS)

ATLAS <- reduceDimension(ATLAS, max_components = 2, num_dim = 6,
                reduction_method = 'tSNE', verbose = T)
ATLAS <- clusterCells(ATLAS, num_clusters = 7)
plot_cell_clusters(ATLAS, 1, 2, color = "group",
    markers = c("CDX2", "APOA1"))
