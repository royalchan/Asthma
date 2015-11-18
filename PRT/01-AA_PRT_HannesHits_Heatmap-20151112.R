############################################################
# TO-DO LIST.
############################################################
DONE. # 1. Plot heatmaps of Hannes Hits as of 11/12/2015.





############################################################
# A. Library Used.
############################################################
library("gplots")
library("ggplot2")
library("RColorBrewer")





############################################################
# B. PLOT HEATMAPS FOR THE HITS (HANNES ANOVA HITS).
############################################################
############################################################
# B.1. Raw Data Loading.
############################################################
# Load data.
dat <- read.table("/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/AA/AA_Proteomics_20150729/Results/20151112-HitsBasedonHannesAnalysis/AA_PRT_Hits_ADJANOVAP0_05_20151112.txt",  header =TRUE, sep = "\t")

# Accession	pval_adj_anova	KN0006	KN0005	KN0008	KN0007	KN0016	KN0015	KN0020	KN0019	REF1R1	KN1034	KN1033	KN1076	KN1075	KN1294	KN1293	KN1496	KN1495	REF2R1	KN1518	KN1517	KN1611	KN1610	KN1804	KN1803	KN1892	KN1891	REF3R1	KN4750	KN4749	KN4784	KN4783	KN4818	KN4817	KN4906	KN4905	REF6R1	KN6116	KN6115	KN6202	KN6201	KN6296	KN6295	KN6506	KN6505	REF8R1_1	KN6116_2	KN6115_2	KN6202_2	KN6201_2	KN6296_2	KN6295_2	KN6506_2	KN6505_2	REF8R1_2
# Diagnosis	pval_adj_anova	DA_A	DA_N	CA_A	CA_A	DA_A	DA_N	CH_N	CH_N	REF1R1	CH_N	CH_N	DA_N	DA_A	CA_A	CA_A	CA_A	CA_A	REF2R1	CA_A	CA_A	CA_A	CA_A	CH_N	CH_N	DA_A	DA_N	REF3R1	CA_A	CA_A	CA_A	CA_A	CA_A	CA_A	DA_N	DA_A	REF6R1	DA_N	DA_A	CH_N	CH_N	CH_N	CH_N	DA_A	DA_N	REF8R1_1	DA_N	DA_A	CH_N	CH_N	CH_N	CH_N	DA_A	DA_N	REF8R1_2

row.names(dat) <- dat$Accession
dat_data <- as.matrix(dat[2:length(dat[,1]),3:length(dat[1,])])
dat_data <- apply(dat_data, 2, as.numeric)
dat_data[is.na(dat_data)] <- 0
nrow <- length(dat_data[,1])
ncol <- length(dat_data[1,])
datrownames <- row.names(dat)[2:length(dat[,1])]
row.names(dat_data) <- datrownames
dat_data[,37] <- (dat_data[,37] + dat_data[,46]) / 2
dat_data[,38] <- (dat_data[,38] + dat_data[,47]) / 2
dat_data[,39] <- (dat_data[,39] + dat_data[,48]) / 2
dat_data[,40] <- (dat_data[,40] + dat_data[,49]) / 2
dat_data[,41] <- (dat_data[,41] + dat_data[,50]) / 2
dat_data[,42] <- (dat_data[,42] + dat_data[,51]) / 2
dat_data[,43] <- (dat_data[,43] + dat_data[,52]) / 2
dat_data[,44] <- (dat_data[,44] + dat_data[,53]) / 2
dat_data <- dat_data[,c(1,2,5,6,13,12,25,26,35,34,38,37,43,44,3,4,14,15,16,17,19,20,21,22,28,29,30,31,32,33,7,8,10,11,23,24,39,40,41,42)]
dim(dat_data)
[1] 214  40

# Heatmap of pval_adj_anova < 0.05.
my_palette <- colorRampPalette(c("turquoise", "gray12", "gold"))(n = 500)
hc_ts <- hclust(dist(dat_data, method = "euclidean"), method = "ward.D")
total_heatmap <- heatmap.2(dat_data, Rowv=as.dendrogram(hc_ts), Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(7,7), density.info = "none", trace = "none", scale = "row", symkey = FALSE, keysize = 0.5, cexRow = 0.6, ColSideColor = c(rep(c("gold","turquoise"), 7), rep(c("red","red"), 8), rep(c("black","black"), 5)))
total_heatmap <- heatmap.2(dat_data[,c(1,3,5,7,9,11,13,2,4,6,8,10,12,14,15:40)], Rowv=as.dendrogram(hc_ts), Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(7,7), density.info = "none", trace = "none", scale = "row", symkey = FALSE, keysize = 0.5, cexRow = 0.6, ColSideColor = c(rep(c("gold"), 7), rep(c("turquoise"), 7), rep(c("red","red"), 8), rep(c("black","black"), 5)))
hc_ts2 <- hclust(dist(t(dat_data), method = "euclidean"), method = "ward.D")
total_heatmap <- heatmap.2(dat_data, Colv=as.dendrogram(hc_ts2), Rowv=as.dendrogram(hc_ts), dendrogram = "both", col = my_palette, margins=c(7,7), density.info = "none", trace = "none", scale = "row", symkey = FALSE, keysize = 0.5, cexRow = 0.6, ColSideColor = c(rep(c("gold","turquoise"), 7), rep(c("red","red"), 8), rep(c("black","black"), 5)))

# Heatmap of pval_adj_anova < 0.01.
n <- length(dat_data[,1])
pv <- as.numeric(as.character(dat[3:length(dat[,1]),2]))
dat_data <- dat_data[(1:n)[pv < 0.01],]
dim(dat_data)
[1] 134  40
my_palette <- colorRampPalette(c("turquoise", "gray12", "gold"))(n = 500)
hc_ts <- hclust(dist(dat_data, method = "euclidean"), method = "ward.D")
total_heatmap <- heatmap.2(dat_data, Rowv=as.dendrogram(hc_ts), Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(7,7), density.info = "none", trace = "none", scale = "row", symkey = FALSE, keysize = 0.5, cexRow = 0.6, ColSideColor = c(rep(c("gold","turquoise"), 7), rep(c("red","red"), 8), rep(c("black","black"), 5)))
total_heatmap <- heatmap.2(dat_data[,c(1,3,5,7,9,11,13,2,4,6,8,10,12,14,15:40)], Rowv=as.dendrogram(hc_ts), Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(7,7), density.info = "none", trace = "none", scale = "row", symkey = FALSE, keysize = 0.5, cexRow = 0.6, ColSideColor = c(rep(c("gold"), 7), rep(c("turquoise"), 7), rep(c("red","red"), 8), rep(c("black","black"), 5)))
hc_ts2 <- hclust(dist(t(dat_data), method = "euclidean"), method = "ward.D")
total_heatmap <- heatmap.2(dat_data, Colv=as.dendrogram(hc_ts2), Rowv=as.dendrogram(hc_ts), dendrogram = "both", col = my_palette, margins=c(7,7), density.info = "none", trace = "none", scale = "row", symkey = FALSE, keysize = 0.5, cexRow = 0.6, ColSideColor = c(rep(c("gold","turquoise"), 7), rep(c("red","red"), 8), rep(c("black","black"), 5)))




