#z-zcore normalization
#Converts expression value of a gene in each sample to z-score based on expression across the samples
data <- as.matrix(read.table("~Data/meta_analysis/log_data.txt",sep="\t",header=T,row.names=1))
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

#Compute and plot PCA of z-score based expression
data2 = read.table("zscore_data.txt", header=T, row.names = 1, sep="\t")
#Renqame header
names(data2) <- c("Spring_6441",	"Spring_6455",	"Spring_6734",	"Spring_6776",
                  "Transplant_6434",	"Transplant_6453",	"Transplant_6902",	"Transplant_H1_1_JM",	"Transplant_H1_2_JM",	"Transplant_H1_1_JM.1",
                  "H9_HIO_1",	"H9_HIO_2",	"H9_HIO_3",
                  "H1_HIO_1",	"H1_HIO_2",	"H1_HIO_3",	"H1_HIO_4_JM",	"H1_HIO_5_JM",
                  "Fetal_SI_2",	"Fetal_SI_3",	"Fetal_SI_4",	"Fetal_SI_5",	"Fetal_SI_6",
                  "Infant_SI_1",	"Infant_SI_14078",	"Infant_SI_2",	"Infant_SI_3",	"Infant_SI_4",	"Infant_SI_5",
                  "Child_SI_1",
                  "Adult_DUO_1",	"Adult_DUO_2",	"Adult_SI_14012",	"Adult_SI_14021",	"Adult_SI_14097",	"Adult_SI_2",	"Adult_SI_3",	"Adult_SI_4",	"Adult_SI_5")

#Select the 10000 genes with the highest variance (To keep consistence with PCA related to Figure 3g)
vars <- apply(data2,1, var,na.rm=TRUE)
vars <- sort(vars,decreasing=TRUE)
vars <- vars[1:10000] #Select 10000 genes with highest variance
data2 <- data2[names(vars),]

#Group setup
Adult <- grep("Adult",colnames(data),ignore.case=T)
Infant <- grep("Infant",colnames(data),ignore.case=T)
Fetal <- grep("Fetal",colnames(data),ignore.case=T)
Transplanted <- grep("Transplant",colnames(data),ignore.case=T)
Spring <- grep("Spring",colnames(data),ignore.case=T)
Child <- grep("Child",colnames(data),ignore.case=T)
HIO_H9 <- grep("H9",colnames(data),ignore.case=T)
HIO_H1 <- grep("H1_HIO",colnames(data),ignore.case=F)


group <- as.factor(c(rep("tHIO+S",times=length(Spring)),
                     rep("tHIO",times=length(Transplanted)),
                     rep("HIO_H9",times=length(HIO_H9)),
                     rep("HIO_H1",times=length(HIO_H1)),
                     rep("Fetal",times=length(Fetal)),
                     rep("Infant",times=length(Infant)),
                     rep("Child",times=length(Child)),
                     rep("Adult",times=length(Adult))))

pheno = read.table("~Data/meta_analysis/pheno_combat.txt", header=T, row.names = 1, sep="\t")
pca <- prcomp(t(data2),scale=TRUE,center=TRUE)
scores <- data.frame(colnames(data2), pca$x[,1:ncol(pca$x)],group)
theme <- theme(legend.position="bottom",
               legend.title=element_blank(),
               legend.background = element_rect(fill="white", size=.5, linetype="dotted"),
               panel.background = element_rect(fill = "white", colour = "black"),
               panel.grid.minor=element_blank(),
               panel.grid.major=element_blank())

#Get proportion of variance explained per component
summary(pca) #Add PC1 and PC2 proportion of variance explained to PCA axis titles

#Plot PCA with predefined theme

PCAplot <-ggplot(scores, aes(x = PC1, y= PC2)) +
  geom_point(shape=21,aes(fill=factor(group)), size=5)  +
  #stat_ellipse(mapping = aes(fill=factor(group))) +
  theme +
  scale_fill_manual(values=c("#F5EA14", "#E6AA22", "#F48220", "#EC268F", "#DD57A4", "#FBF7C9", "#E5CBE2", "#8A4B9C")) +
  xlim(-150, 150) +
  ylim(-100, 100) +
  labs(x = "PC1 (53.48%)", y = "PC2 (39.38%)")
plot(PCAplot)

#Save and export plot
pdf(file = "zScore_PCA_Plot.pdf", width = 10, height = 10) # defaults to 7 x 7 inches
plot(PCAplot)
dev.off()

png(file = "zScore_PCA_Plot.png", width = 600, height = 600) # defaults to 7 x 7 inches
plot(PCAplot)
dev.off()
