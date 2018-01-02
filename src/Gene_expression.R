#Used packages
library(DESeq2)
library(pcaExplorer)
library(DEFormats)
library(ggplot2)

setwd("~Data/gene_analysis/Raw_Counts.txt") # Set working directory
count = read.table("Raw_Counts.txt", header=T, row.names = 1, sep="\t", check.names = FALSE) # Open Raw counts from featureCounts routine

#Group setup
Adult <- grep("Adult",colnames(count),ignore.case=T)
Infant <- grep("Infant",colnames(count),ignore.case=T)
Spring <- grep("Spring",colnames(count),ignore.case=T)
tHIO <- grep("tHIO",colnames(count),ignore.case=T)
group <- as.factor(c(rep("Spring",times=length(Spring)), rep("tHIO",times=length(tHIO)), rep("Infant",times=length(Infant)), (rep("Adult",times=length(Adult)))))
pheno = read.table("pheno.txt", header=T, row.names = 1, sep="\t")

#Create Large DGEList
dge = DGEList(count, group = group)
dge$samples$batch=pheno$batch
dge

#Create Large DESeqDataSet
dds = as.DESeqDataSet(dge)
dds = as.DESeqDataSet(dge, design=~ batch + group)
#vsd = vst(dds, blind=FALSE, fitType = "parametric") #Variance Stabilization Transformation
rld = rlogTransformation(dds) #rlogTransformation on DESeqTransform
dds <- DESeq(dds)
res <- results(dds)
annotation <- data.frame(gene_name = rownames(dds), row.names = rownames(dds), stringsAsFactors = FALSE)

#Compute PCA
pcaobj <- prcomp(t(assay(rld)))
pcascree(pcaobj,type="pev", title="Proportion of explained proportion of variance")
res_pc <- correlatePCs(pcaobj,colData(dds))

#Extract top 10 genes per top/bottom loadings in each PCs
#Fix(hi-loadings) i.e. barplot(geneloadings_extreme, las = 2, col = c(rep("#BD202E",
        topN), rep("#2E368E", topN)), main = paste0(title, "PC",
        whichpc))
#Plot and export image for PC1
png('Loadings_PC1.png', width = 800, height = 350)
hi_loadings(pcaobj, whichpc = 1, topN = 10,annotation = annotation)
dev.off()
#Plot and export image for PC2
png('Loadings_PC2.png', width = 800, height = 350)
hi_loadings(pcaobj, whichpc = 2, topN = 10,annotation = annotation)
dev.off()
#Plot and export image for PC3
png('Loadings_PC2.png', width = 800, height = 350)
hi_loadings(pcaobj, whichpc = 3, topN = 10,annotation = annotation)
dev.off()

#Compute limmaquickpca2go
library(org.Hs.eg.db)
goquick <- limmaquickpca2go(rld,
                                   pca_ngenes = 10000,
                                   inputType = "SYMBOL",
                                   organism = "Hs")

head(goquick$PC1$posLoad)
head(goquick$PC1$negLoad)
#repeat for each PCs
