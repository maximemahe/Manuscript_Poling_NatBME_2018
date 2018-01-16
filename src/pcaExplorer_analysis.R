#Used packages
library(DESeq2)
library(pcaExplorer)
library(DEFormats)
library(ggplot2)

#Open Raw counts from featureCounts routine
count = read.table("~Data/meta_analysis/raw_countstable.txt", header=T, row.names = 1, sep="\t")

#Group setup
Adult <- grep("Adult",colnames(count),ignore.case=T)
Infant <- grep("Infant",colnames(count),ignore.case=T)
Fetal <- grep("Fetal",colnames(count),ignore.case=T)
Transplanted <- grep("Transplant",colnames(count),ignore.case=T)
Spring <- grep("Spring",colnames(count),ignore.case=T)
Child <- grep("Child",colnames(count),ignore.case=T)
HIO_H1 <- grep("H1_HIO",colnames(count),ignore.case=F)
HIO_H9 <- grep("H9_HIO",colnames(count),ignore.case=F)

group <- as.factor(c(rep("tHIO+S",times=length(Spring)), rep("tHIO",times=length(Transplanted)), rep("HIO_H9",times=length(HIO_H9)), rep("HIO_H1",times=length(HIO_H1)), rep("Fetal",times=length(Fetal)),rep("Infant",times=length(Infant)),rep("Child",times=length(Child)),rep("Adult",times=length(Adult))))
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
#repeat for each PC
