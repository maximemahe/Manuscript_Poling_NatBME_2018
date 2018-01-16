#ComBat correction
library(devtools)
library(Biobase)
library(sva)

count = read.table("~Data/meta_analysis/log_data.txt", header=T, row.names = 1, sep="\t")
pheno = read.table("~Data/meta_analysis/pheno_combat.txt", header=T, row.names = 1, sep="\t")

#getVar <- apply(count, 1, var)
#param <- 2
#data= count[getVar > param & !is.na(getVar), ]

#Select the 10000 genes with the highest variance (To keep consistence with PCA related to Figure 3g)
vars <- apply(count,1, var,na.rm=TRUE)
vars <- sort(vars,decreasing=TRUE)
vars <- vars[1:10000] #Select 10000 genes with highest variance
data <- count[names(vars),]

batch = pheno$batch
modcombat = model.matrix(~1, data=pheno)
modgroup = model.matrix(~group, data=pheno)
combat_edata = ComBat(dat=data, batch=batch, mod=modgroup, par.prior=TRUE, prior.plots=FALSE )

Adult <- grep("Adult",colnames(count),ignore.case=T)
Infant <- grep("Infant",colnames(count),ignore.case=T)
Fetal <- grep("Fetal",colnames(count),ignore.case=T)
Transplanted <- grep("Transplant",colnames(count),ignore.case=T)
Spring <- grep("Spring",colnames(count),ignore.case=T)
Child <- grep("Child",colnames(count),ignore.case=T)
HIO_H9 <- grep("H9",colnames(count),ignore.case=T)
HIO_H1 <- grep("H1_HIO",colnames(count),ignore.case=F)

group <- as.factor(c(rep("tHIO+S",times=length(Spring)),
                     rep("tHIO",times=length(Transplanted)),
                     rep("HIO_H9",times=length(HIO_H9)),
                     rep("HIO_H1",times=length(HIO_H1)),
                     rep("Fetal",times=length(Fetal)),
                     rep("Infant",times=length(Infant)),
                     rep("Child",times=length(Child)),
                     rep("Adult",times=length(Adult))))

pca <- prcomp(t(combat_edata),scale=TRUE,center=TRUE)
scores <- data.frame(colnames(combat_edata), pca$x[,1:ncol(pca$x)],group)

theme <- theme(legend.position="bottom",
               legend.title=element_blank(),
               legend.background = element_rect(fill="white", size=.5, linetype="dotted"),
               panel.background = element_rect(fill = "white", colour = "black"),
               panel.grid.minor=element_blank(),
               panel.grid.major=element_blank())

#Get proportion of variance explained per component
summary(pca) #Add PC1 and PC2 proportion of variance explained to PCA axis titles

PCAplot <-ggplot(scores, aes(x = PC1, y= PC2)) +
          geom_point(shape=21,aes(fill=factor(group)), size=5)  +
          #stat_ellipse(mapping = aes(fill=factor(group))) +
          theme +
          scale_fill_manual(values=c("#F5EA14", "#E6AA22", "#F48220", "#EC268F", "#DD57A4", "#FBF7C9", "#E5CBE2", "#8A4B9C")) +
          xlim(-100, 150) +
          ylim(-175, 100) +
          labs(x = "PC1 (58.26%)", y = "PC2 (46.50%)")
plot(PCAplot)

#Save and export plot
pdf(file = "ComBat_PCA_Plot.pdf", width = 10, height = 10) # defaults to 7 x 7 inches
plot(PCAplot)
dev.off()

png(file = "ComBat_PCA_Plot.png", width = 600, height = 600) # defaults to 7 x 7 inches
plot(PCAplot)
dev.off()
