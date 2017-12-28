#Used packages
library(ggplot2)

#Open filtered counts table
count = read.table("FilteredValues.txt", header=T, row.names = 1, sep="\t", check.names = FALSE)

#Group setup
Adult <- grep("Adult",colnames(count),ignore.case=T)
Infant <- grep("Infant",colnames(count),ignore.case=T)
Spring <- grep("Spring",colnames(count),ignore.case=T)
tHIO <- grep("tHIO",colnames(count),ignore.case=T)
group <- as.factor(c(rep("Adult",times=length(Adult)), rep("Infant",times=length(Infant)), rep("Spring",times=length(Spring)), rep("tHIO",times=length(tHIO))))

#Compute PCA
pca <- prcomp(t(count),scale=TRUE,center=TRUE)
scores <- data.frame(colnames(count), pca$x[,1:ncol(pca$x)],group)

#Get proportion of variance explained per component
summary(pca) #Add PC1 and PC2 proportion of variance explained to PCA axis titles

#Set plot theme
theme <- theme(legend.position="bottom",
        legend.title=element_blank(),
        legend.background = element_rect(fill="white", size=.5, linetype="dotted"),
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank())

#Plot PCA with pre-defined theme
qplot(x=PC1, y=PC2, data=scores) +
        theme +
        scale_fill_manual(values=c("#F5EA14", "#FBF7C9", "#8A4B9C", "#E5CBE2")) +
        geom_point(shape=21,aes(fill=factor(group)), size=5)  +
        xlim(-150, 150) +
        ylim(-150, 150) +
        labs(x = "PC1 (36.14%)", y = "PC2 (17.33%)")

#Save and export plot
ggsave("PCA_Plot.png", plot = last_plot(), width = 10, height = 10, units = "cm",
          dpi = 300)
