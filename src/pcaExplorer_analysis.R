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
HIO_H1 <- grep("HIO_H1",colnames(count),ignore.case=F)

group <- as.factor(c(rep("tHIO+S",times=length(Spring)), rep("tHIO",times=length(Transplanted)), rep("HIO_H1",times=length(HIO_H1)),rep("Fetal",times=length(Fetal)),rep("Infant",times=length(Infant)),rep("Child",times=length(Child)),rep("Adult",times=length(Adult))))
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
#Lunch Shinny app #pcaExplorer(dds = dds, rlt = vsd, annotation = annotation)

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

disp_table <- dispersionTable(ATLAS)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 1)
ATLAS <- setOrderingFilter(ATLAS, unsup_clustering_genes$gene_id)
plot_ordering_genes(ATLAS)

ATLAS <- reduceDimension(ATLAS, max_components = 2, num_dim = 2,
            reduction_method = 'tSNE',
            residualModelFormulaStr = "~type + num_genes_expressed",
            verbose = T)


*************************************************
SEURAT

nbt.data=read.table("combat.txt",sep="\t",header=TRUE,row.names=1)
nbt.data=log(nbt.data+1)
corner(nbt.data)
nbt=new("seurat",raw.data=nbt.data)
nbt=setup(nbt,project="NBT",min.cells = 1,names.field = 2,names.delim = "_",min.genes = 1000,is.expr=0.5,)
nbt=mean.var.plot(nbt,y.cutoff = 0.1,x.low.cutoff = 0.1,fxn.x = expMean,fxn.y = logVarDivMean)
length(nbt@var.genes)
nbt=pca(nbt,do.print=FALSE)
pca.plot(nbt,1,2,pt.size = 2)

viz.pca(nbt,1:2)
pcHeatmap(nbt,pc.use = 1,do.balanced = FALSE)

nbt=jackStraw(nbt,num.replicate = 200,do.print = FALSE)
jackStrawPlot(nbt,PCs = 1:12)

nbt=project.pca(nbt,do.print=FALSE)
pcHeatmap(nbt,pc.use = 1,use.full = TRUE,do.balanced = TRUE,remove.key = TRUE)

nbt.sig.genes=pca.sig.genes(nbt,1:5,pval.cut = 1e-4,max.per.pc = 200)
length(nbt.sig.genes)
nbt=pca(nbt,pc.genes=nbt.sig.genes,do.print = FALSE)
nbt=jackStraw(nbt,num.replicate = 200,do.print = FALSE)
jackStrawPlot(nbt,PCs = 1:15)

nbt=run_tsne(nbt,dims.use = 1:5,max_iter=2000)
tsne.plot(nbt,pt.size = 4) #Run on PCs



*****
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
ATLAS <- clusterCells(ATLAS, num_clusters = 8)
plot_cell_clusters(ATLAS, 1, 2, color = "group",
    markers = c("CDX2", "APOA1"))


#######################################################################################

count = read.table("combat_data.txt", header=T, row.names = 1, sep="\t")
> #Group setup
> Adult <- grep("Adult",colnames(count),ignore.case=T)
> Infant <- grep("Infant",colnames(count),ignore.case=T)
> Fetal <- grep("Fetal",colnames(count),ignore.case=T)
> Transplanted <- grep("Transplant",colnames(count),ignore.case=T)
> Spring <- grep("Spring",colnames(count),ignore.case=T)
> Child <- grep("Child",colnames(count),ignore.case=T)
> HIO_H1 <- grep("HIO_H1",colnames(count),ignore.case=F)
>
> group <- as.factor(c(rep("tHIO+S",times=length(Spring)), rep("tHIO",times=length(Transplanted)), rep("HIO_H1",times=length(HIO_H1)),rep("Fetal",times=length(Fetal)),rep("Infant",times=length(Infant)),rep("Child",times=length(Child)),rep("Adult",times=length(Adult))))
> pheno = read.table("pheno.txt", header=T, row.names = 1, sep="\t")
> pca <- prcomp(t(count),scale=TRUE,center=TRUE)
> scores <- data.frame(colnames(count), pca$x[,1:ncol(pca$x)],group)
> theme <- theme(legend.position="bottom",
+                legend.title=element_blank(),
+                legend.background = element_rect(fill="white", size=.5, linetype="dotted"),
+                panel.background = element_rect(fill = "white", colour = "black"),
+                panel.grid.minor=element_blank(),
+                panel.grid.major=element_blank())
> qplot(x=PC1, y=PC2, data=scores) +
+     theme +
+     scale_fill_manual(values=c("#F5EA14", "#FBF7C9", "#8A4B9C", "#E5CBE2")) +
+     geom_point(shape=21,aes(fill=factor(group)), size=5)  +
+     xlim(-150, 150) +
+     ylim(-150, 150) +
+     labs(x = "PC1 (36.14%)", y = "PC2 (17.33%)")
Error: Insufficient values in manual scale. 7 needed but only 4 provided.
> qplot(x=PC1, y=PC2, data=scores) +
+     theme +
+     scale_fill_manual(values=c("#F5EA14", "#FBF7C9", "#8A4B9C", "#E5CBE2", "#F5EA14", "#FBF7C9", "#8A4B9C", "#E5CBE2")) +
+     geom_point(shape=21,aes(fill=factor(group)), size=5)  +
+     xlim(-150, 150) +
+     ylim(-150, 150) +
+     labs(x = "PC1 (36.14%)", y = "PC2 (17.33%)")

** Expression to Zscore
#Converts expression value of a gene in each sample to z-score based on expression across the samples
file <- data
data <- as.matrix(read.table("normalized_counstable",sep="\t",header=T,row.names=1))
cat("\n")
head(data)
zscore_data <- c()
for(i in 1:nrow(data))
{
  mean_row <- mean(as.numeric(data[i,]),na.rm=TRUE)
  sd_row <- sd(as.numeric(data[i,]),na.rm=TRUE)
  zscore_row <- c()
  for(j in 1:ncol(data))
  {
   zscore_each <- ((data[i,j]-mean_row)/sd_row)
   zscore_row <-append(zscore_row,zscore_each)
  }
  zscore_data<-rbind(zscore_data,zscore_row)
  Sys.sleep(0.00005)
  cat(paste0("Zscore Calculations for Gene: ", i, " Completed", sep="\r"))
}
col <- cbind(row.names(data),zscore_data)
write.table(col,"zscore_data.txt",sep="\t",row.names=F,quote=F)

#Plot PCA
data = read.table("zscore_data.txt", header=T, row.names = 1, sep="\t")

#Group setup
Adult <- grep("Adult",colnames(data),ignore.case=T)
Infant <- grep("Infant",colnames(data),ignore.case=T)
Fetal <- grep("Fetal",colnames(data),ignore.case=T)
Transplanted <- grep("Transplant",colnames(data),ignore.case=T)
Spring <- grep("Spring",colnames(data),ignore.case=T)
Child <- grep("Child",colnames(data),ignore.case=T)
HIO_H1 <- grep("HIO_H1",colnames(data),ignore.case=F)

group <- as.factor(c(rep("tHIO+S",times=length(Spring)), rep("tHIO",times=length(Transplanted)), rep("HIO_H1",times=length(HIO_H1)),rep("Fetal",times=length(Fetal)),rep("Infant",times=length(Infant)),rep("Child",times=length(Child)),rep("Adult",times=length(Adult))))
pheno = read.table("pheno_combat.txt", header=T, row.names = 1, sep="\t")
pca <- prcomp(t(data),scale=TRUE,center=TRUE)
scores <- data.frame(colnames(data), pca$x[,1:ncol(pca$x)],group)
theme <- theme(legend.position="bottom",
              legend.title=element_blank(),
              legend.background = element_rect(fill="white", size=.5, linetype="dotted"),
              panel.background = element_rect(fill = "white", colour = "black"),
              panel.grid.minor=element_blank(),
              panel.grid.major=element_blank())

#Get proportion of variance explained per component
summary(pca) #Add PC1 and PC2 proportion of variance explained to PCA axis titles

qplot(x=PC1, y=PC2, data=scores) +
     theme +
     scale_fill_manual(values=c("#F5EA14", "#FBF7C9", "#8A4B9C", "#E5CBE2", "#F5EA14", "#FBF7C9", "#8A4B9C", "#E5CBE2")) +
     geom_point(shape=21,aes(fill=factor(group)), size=5)  +
     xlim(-150, 150) +
     ylim(-150, 150) +
     labs(x = "PC1 (36.14%)", y = "PC2 (17.33%)")
