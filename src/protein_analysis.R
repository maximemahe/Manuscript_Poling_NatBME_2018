library("factoextra")

data <- read.table("~/Data/ProteinArray.txt", sep = "\t", header = T, row.names = 1)
df <- data[1:8, 1:32]
> df <- scale(df)
