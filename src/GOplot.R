library(GOplot)
#Import gene list and GO enrichment data from toppgene suite
#Pathways upregulated in tHIO+S
GO <- read.table("GO-Pathway.txt", sep = "\t", header = TRUE)
genelist <- read.table("genelist_thiospring-thio.txt", sep = "\t", header = TRUE)
circ <- circle_dat(GO, genelist)
genes <- read.table("selected_genes.txt", sep = "\t", header = TRUE)
process <- read.table("selected_process.txt", sep = "\t", header = TRUE)
process <- as.matrix(process)
chord <- chord_dat(circ, genes, process)
GOChord(chord, limit = c(2, 0),
              gene.order = 'logFC',
              ribbon.col = c("#F79420", "#FAAF41", "#F7EC33", "#FBF8CD", "#FFFFFF", "#EECDE1", "#6B52A1", "#D94397"),
              gene.space = 0.25,
              gene.size = 8,
              space = 0.02)
ggsave("GOPathway.png", plot = last_plot(), width = 40, height = 45, units = "cm", dpi = 200)

#Pathways upregulated in tHIO+S
GO <- read.table("tHIOS_unique_all.txt", sep = "\t", header = TRUE)
genelist <- read.table("genelist_tHIOSunique.txt", sep = "\t", header = TRUE)
circ <- circle_dat(GO, genelist)
genes <- read.table("genes.txt", sep = "\t", header = TRUE)
process <- read.table("Processes.txt", sep = "\t", header = TRUE)
process <- as.matrix(process)
chord <- chord_dat(circ, genes, process)
GOChord(chord, limit = c(2, 0),
              gene.order = 'logFC',
              ribbon.col = c("#F79420", "#FAAF41", "#F7EC33", "#FBF8CD", "#FFFFFF", "#EECDE1", "#6B52A1", "#D94397"),
              gene.space = 0.25,
              gene.size = 8,
              space = 0.02)
ggsave("GOtHIOS.png", plot = last_plot(), width = 40, height = 45, units = "cm", dpi = 200)
