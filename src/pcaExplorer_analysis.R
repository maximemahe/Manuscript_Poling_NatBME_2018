#Used packages
library(DESeq2)
library(pcaExplorer)
library(DEFormats)
library(ggplot2)

#Open Raw counts from featureCounts routine
count = read.table("~Data/meta_analysis/count.txt", header=T, row.names = 1, sep="\t")

#Group setup
Adult <- grep("Adult",colnames(count),ignore.case=T)
Infant <- grep("Infant",colnames(count),ignore.case=T)
Fetal <- grep("Fetal",colnames(count),ignore.case=T)
Transplanted <- grep("Transplant",colnames(count),ignore.case=T)
Spring <- grep("Spring",colnames(count),ignore.case=T)
Child <- grep("Child",colnames(count),ignore.case=T)
HIO_H9 <- grep("H9",colnames(count),ignore.case=T)
HIO_H1 <- grep("H1_HIO",colnames(count),ignore.case=F)

group <- as.factor(c(rep("tHIO+S",times=length(Spring)), rep("tHIO",times=length(Transplanted)),rep("HIO_H9",times=length(HIO_H9)), rep("HIO_H1",times=length(HIO_H1)),rep("Fetal",times=length(Fetal)),rep("Infant",times=length(Infant)),rep("Child",times=length(Child)),rep("Adult",times=length(Adult))))
pheno = read.table("pheno.txt", header=T, row.names = 1, sep="\t")

#Create Large DGEList
dge = DGEList(count, group = group)
dge$samples$batch=pheno$batch
dge

#Create Large DESeqDataSet
dds = as.DESeqDataSet(dge)
dds = as.DESeqDataSet(dge, design=~ batch + group)
vsd = vst(dds, blind=FALSE, fitType = "parametric") #Variance Stabilization Transformation
#rld = rlogTransformation(dds) #rlogTransformation on DESeqTransform (option)
dds <- DESeq(dds)
res <- results(dds)
annotation <- data.frame(gene_name = rownames(dds), row.names = rownames(dds), stringsAsFactors = FALSE)

#Compute PCA
pcaobj <- prcomp(t(assay(vsd))) #Use VST
pcascree(pcaobj,type="pev", title="Proportion of explained proportion of variance")
res_pc <- correlatePCs(pcaobj,colData(dds))

#Plot PCA (10000 genes) and export image
png('PCA_10000genes.png', width = 750, height = 600)
pcaplot(vsd, intgroup = "group", ntop = 10000, returnData = FALSE, pcX = 1, pcY = 2, title = "PCA on 10000 genes", text_labels = FALSE, point_size = 7, ellipse = TRUE, ellipse.prob = 0.95) +
    theme(axis.text.y   = element_text(size=12), axis.text.x   = element_text(size=12), axis.title.y  = element_text(size=14), axis.title.x  = element_text(size=14),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    panel.border = element_rect(colour = "black", fill=NA, size=1)) +
    scale_color_manual(values= c("#F5EA14", "#E6AA22", "#F48220", "#EC268F", "#DD57A4", "#FBF7C9", "#E5CBE2", "#8A4B9C")) +
    geom_point(shape=1, size=7, colour="black")
dev.off()

#Plot PCA (5000 genes) and export image
png('PCA_5000genes.png', width = 750, height = 600)
pcaplot(vsd, intgroup = "group", ntop = 5000, returnData = FALSE, pcX = 1, pcY = 2, title = "PCA on 5000 genes", text_labels = FALSE, point_size = 7, ellipse = TRUE, ellipse.prob = 0.95) +
    theme(axis.text.y   = element_text(size=12), axis.text.x   = element_text(size=12), axis.title.y  = element_text(size=14), axis.title.x  = element_text(size=14),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    panel.border = element_rect(colour = "black", fill=NA, size=1)) +
    scale_color_manual(values= c("#F5EA14", "#E6AA22", "#F48220", "#EC268F", "#DD57A4", "#FBF7C9", "#E5CBE2", "#8A4B9C")) +
    geom_point(shape=1, size=7, colour="black")
dev.off()

#Plot PCA (2500 genes) and export image
png('PCA_2500genes.png', width = 750, height = 600)
pcaplot(vsd, intgroup = "group", ntop = 2500, returnData = FALSE, pcX = 1, pcY = 2, title = "PCA on 2500 genes", text_labels = FALSE, point_size = 7, ellipse = TRUE, ellipse.prob = 0.95) +
    theme(axis.text.y   = element_text(size=12), axis.text.x   = element_text(size=12), axis.title.y  = element_text(size=14), axis.title.x  = element_text(size=14),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    panel.border = element_rect(colour = "black", fill=NA, size=1)) +
    scale_color_manual(values= c("#F5EA14", "#E6AA22", "#F48220", "#EC268F", "#DD57A4", "#FBF7C9", "#E5CBE2", "#8A4B9C")) +
    geom_point(shape=1, size=7, colour="black")
dev.off()

#Plot PCA (500 genes) and export image
png('PCA_500genes.png', width = 750, height = 600)
pcaplot(vsd, intgroup = "group", ntop = 500, returnData = FALSE, pcX = 1, pcY = 2, title = "PCA on 500 genes", text_labels = FALSE, point_size = 7, ellipse = TRUE, ellipse.prob = 0.95) +
    theme(axis.text.y   = element_text(size=12), axis.text.x   = element_text(size=12), axis.title.y  = element_text(size=14), axis.title.x  = element_text(size=14),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    panel.border = element_rect(colour = "black", fill=NA, size=1)) +
    scale_color_manual(values= c("#F5EA14", "#E6AA22", "#F48220", "#EC268F", "#DD57A4", "#FBF7C9", "#E5CBE2", "#8A4B9C")) +
    geom_point(shape=1, size=7, colour="black")
dev.off()

#Extract top 10 genes per top/bottom loadings in each PCs
#Fix(hi-loadings) i.e. barplot(geneloadings_extreme, las = 2, col = c(rep("#BD202E",
        topN), rep("#2E368E", topN)), main = paste0(title, "PC",
        whichpc))
#Plot and export image for PC1
png('Loadings_PC1.png', width = 800, height = 400)
hi_loadings(pcaobj, whichpc = 1, topN = 10,annotation = annotation)
dev.off()
#Plot and export image for PC2
png('Loadings_PC2.png', width = 800, height = 400)
hi_loadings(pcaobj, whichpc = 2, topN = 10,annotation = annotation)
dev.off()
#Plot and export image for PC3
png('Loadings_PC3.png', width = 800, height = 400)
hi_loadings(pcaobj, whichpc = 3, topN = 10,annotation = annotation)
dev.off()

#Compute limmaquickpca2go
library(org.Hs.eg.db)
goquick <- limmaquickpca2go(vsd,
                                   pca_ngenes = 10000,
                                   inputType = "SYMBOL",
                                   organism = "Hs")

head(goquick$PC1$posLoad)
head(goquick$PC1$negLoad)
#repeat for each PCs






***********************************************************************************************
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
