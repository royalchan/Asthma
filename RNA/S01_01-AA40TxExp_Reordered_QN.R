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


dat <- read.table("/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/AA/AA_Methyl_Seq/Analyses/20141208-From_George/20150301-MetHitGene_Expression/AA40TxExp_Reordered.txt", sep = "\t", header = TRUE)
# Column Names.
# Transcript_ID_Chr_L_R	KNAA3B	KNAA3A	KNAA4A	KNAA4B	KNAA6A	KNAA6B	KNAA8B	KNAA8A	KN-C0001	KN-C0002	KN6	KN5	KN0016	KN0015	KN1075	KN1076	KN1892	KN1891	KN4002	KN4001	KN4905	KN4906	KN6115	KN6116	KN6506	KN6505	KN6857	KN6858	KN0007	KN0008	KN1293	KN1294	KN1495	KN1496	KN1517	KN1518	KN1610	KN1611	KN1711	KN1712	KN2027	KN2028	KN4091	KN4092	KN4149	KN4150	KN4749	KN4750	KN4783	KN4784	KN4817	KN4818	KN5235	KN5236	KN6335	KN6336	KN0019	KN0020	KN1033	KN1034	KN1803	KN1804	KN2195	KN2196	KN4117	KN4118	KN4175	KN4176	KN4473	KN4474	KN5261	KN5262	KN6035	KN6036	KN6201	KN6202	KN6295	KN6296	KN6709	KN6710

row.names(dat) <- dat$Transcript_ID_Chr_L_R
dat_data <- as.matrix(dat[,2:length(dat[1,])])
nrow <- length(dat_data[,1])
ncol <- length(dat_data[1,])
datrownames <- row.names(dat)
datcolnames <- colnames(dat)[2:length(colnames(dat))]
dat_data_log <- log2(dat_data + 0.01)


# Quantile Normalization.
total_matrix <- data.matrix(dat_data_log)
total_matrix_qn <- normalize.quantiles(total_matrix)
row.names(total_matrix_qn) <- row.names(dat_data)
colnames(total_matrix_qn) <- colnames(dat_data)

min(total_matrix_qn)
[1] -6.643856
max(total_matrix_qn)
[1] 13.5963

# Output QN Log2 Data.
write.table(total_matrix_qn, file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/AA/AA_Methyl_Seq/Analyses/20141208-From_George/20150301-MetHitGene_Expression/AA40TxExp_Reordered_QN.txt", sep = "\t", row.names = TRUE, col.names = TRUE, append = FALSE)




