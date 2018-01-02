#Used packages
library(ggplot2)
library(reshape2)
library(ggdendrogram)

#Open filtered counts table
cormat <- read.table("~Data/gene_analysis/CorMatrix.txt", header=T, sep="\t", check.names = FALSE, row.names = 1)
cormat <- as.matrix(cormat)
cormat<- round(cormat, 2)

#Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
    cormat[lower.tri(cormat)]<- NA
    return(cormat)
  }

upper_tri <- get_upper_tri(cormat)
melted_cormat <- melt(upper_tri, na.rm = TRUE)

#Heatmap
ggheatmap <- ggplot(data = melted_cormat, aes(Var1, Var2, fill = value)) +
    geom_raster() +
    scale_fill_gradient2(low = "#712B8D", high = "#F5EA14", mid = "#030203",
      midpoint = 0.7, limit = c(0.5,1), space = "Lab",name="Coefficient value") +
    theme_minimal() + # minimal theme
    theme(axis.text.x = element_blank()) +
    coord_fixed()

#Add correlation coefficients on the heatmap
corheatmap <- ggheatmap +
              geom_text(aes(Var1, Var2, label = value), color = "grey50", size = 3) +
              theme(
                    axis.title.x = element_blank(),
                    axis.title.y = element_blank(),
                    panel.grid.major = element_blank(),
                    panel.border = element_blank(),
                    panel.background = element_blank(),
                    axis.ticks = element_blank(),
                    legend.position = c(0.75, 0.1),
                    legend.direction = "horizontal") +
                    guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                                  title.position = "top", title.hjust = 0.5)) +
                    labs(title="Spearman correlation matrix")

#Save and export plot
ggsave("Corr_Plot.png", plot = last_plot(), width = 10, height = 10, units = "cm", dpi = 300)

#Compute dendrogram
dend <- hclust(dist(cormat), method = "complete")
dend_plot <- ggdendrogram(dend) +
            labs(title="Cluster dendrogram")
ggsave("Dend_Plot.png", plot = dend_plot, width = 10, height = 10, units = "cm", dpi = 300)
