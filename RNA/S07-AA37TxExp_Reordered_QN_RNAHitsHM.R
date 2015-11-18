library("plyr")
library("ggplot2")
library("reshape")
library("gplots")
library("preprocessCore")
library("coin")
library("magic")
library("RColorBrewer")
library("limma")
library("LDheatmap")





# Category Info.
Cat <- c("DISAA", "DISNM", "DISAA", "DISNM", "DISAA", "DISNM", "DISAA", "DISNM", "DISAA", "DISNM", "DISAA", "DISNM", "DISAA", "DISNM", "DISAA", "DISNM", "DISAA", "DISNM", "DISAA", "DISNM", "DISAA", "DISNM", "DISAA", "DISNM", "CONAA", "CONAA", "CONAA", "CONAA", "CONAA", "CONAA", "CONAA", "CONAA", "CONAA", "CONAA", "CONAA", "CONAA", "CONAA", "CONAA", "CONAA", "CONAA", "CONAA", "CONAA", "CONAA", "CONAA", "CONAA", "CONAA", "CONAA", "CONAA", "CONAA", "CONAA", "CONNM", "CONNM", "CONNM", "CONNM", "CONNM", "CONNM", "CONNM", "CONNM", "CONNM", "CONNM", "CONNM", "CONNM", "CONNM", "CONNM", "CONNM", "CONNM", "CONNM", "CONNM", "CONNM", "CONNM", "CONNM", "CONNM", "CONNM", "CONNM")
Cat2 <- c("AA", "NM", "AA", "NM", "AA", "NM", "AA", "NM", "AA", "NM", "AA", "NM", "AA", "NM", "AA", "NM", "AA", "NM", "AA", "NM", "AA", "NM", "AA", "NM", "AA", "AA", "AA", "AA", "AA", "AA", "AA", "AA", "AA", "AA", "AA", "AA", "AA", "AA", "AA", "AA", "AA", "AA", "AA", "AA", "AA", "AA", "AA", "AA", "AA", "AA", "NM", "NM", "NM", "NM", "NM", "NM", "NM", "NM", "NM", "NM", "NM", "NM", "NM", "NM", "NM", "NM", "NM", "NM", "NM", "NM", "NM", "NM", "NM", "NM")
Cat3 <- c("DISAA", "DISNM", "DISAA", "DISNM", "DISAA", "DISNM", "DISAA", "DISNM", "DISAA", "DISNM", "DISAA", "DISNM", "DISAA", "DISNM", "DISAA", "DISNM", "DISAA", "DISNM", "DISAA", "DISNM", "DISAA", "DISNM", "DISAA", "DISNM", "CON1", "CON2", "CON1", "CON2", "CON1", "CON2", "CON1", "CON2", "CON1", "CON2", "CON1", "CON2", "CON1", "CON2", "CON1", "CON2", "CON1", "CON2", "CON1", "CON2", "CON1", "CON2", "CON1", "CON2", "CON1", "CON2", "CON1", "CON2", "CON1", "CON2", "CON1", "CON2", "CON1", "CON2", "CON1", "CON2", "CON1", "CON2", "CON1", "CON2", "CON1", "CON2", "CON1", "CON2", "CON1", "CON2", "CON1", "CON2", "CON1", "CON2")

# RNA-SEQ HITS.
dat <- read.table("/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/AA/AA_RNA_Seq_20140630/20140630-NewestRound/Analyses/Combined/AA37TxExp_RNA_Hits_Reordered_QN_GeneSymbol.txt", header = TRUE, sep = "\t")
dat_data <- as.matrix(dat[,3:length(dat[1,])])
row.names(dat_data) <- dat[,1]

# RNA-Seq Hits Heatmap.
my_palette <- colorRampPalette(c("turquoise", "gray12", "gold"))(n = 500)
hc_ts <- hclust(dist(t(dat_data), method = "euclidean"), method = "ward.D2")
total_heatmap <- heatmap.2(dat_data, Colv=as.dendrogram(hc_ts), Rowv=FALSE, dendrogram = "column", col = my_palette, margins=c(5,16), density.info = "none", trace = "none", scale = "none", keysize = 0.5, ColSideColors = c(rep(c("gold","turquoise"),12), rep("red",26), rep("black",24)), cexRow = 0.7)
total_heatmap <- heatmap.2(dat_data, Colv=as.dendrogram(hc_ts), Rowv=FALSE, dendrogram = "column", col = my_palette, margins=c(5,16), density.info = "none", trace = "none", scale = "row", keysize = 0.5, ColSideColors = c(rep(c("gold","turquoise"),12), rep("red",26), rep("black",24)), cexRow = 0.7)

# Heatmap of differences.
dat_diff <- matrix(nrow = 112, ncol = 74)
colnames(dat_diff) <- colnames(dat_data)
row.names(dat_diff) <- row.names(dat_data)
for (i in seq(1, 73, by = 2)) {
	dat_diff[,i] <- (dat_data[,i] - dat_data[,i + 1])/2
	dat_diff[,i + 1] <- - (dat_data[,i] - dat_data[,i + 1])/2	
}
hc_ts <- hclust(dist(t(dat_diff), method = "euclidean"), method = "ward.D2")
total_heatmap <- heatmap.2(dat_diff, Colv=as.dendrogram(hc_ts), Rowv=FALSE, dendrogram = "column", col = my_palette, margins=c(5,16), density.info = "none", trace = "none", scale = "none", keysize = 0.5, ColSideColors = c(rep(c("gold","turquoise"),12), rep("red",26), rep("black",24)), cexRow = 0.7)




