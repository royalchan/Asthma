############################################################
# TO-DO LIST.
############################################################
# Calculate Spearman Correlation Coefficient of Combined AA RNA Data.





############################################################
# A. Library Used.
############################################################
############################################################
# A.1. Full Libraries.
############################################################
library("gplots")
library("ggplot2")
library("reshape2")
library("preprocessCore")
library("coin")
library("magic")
library("RColorBrewer")
library("limma")
library("LDheatmap")
library("MASS")
library("lomb")
library("cts")
library("Sushi")
library("heatmap3")
library("geneplotter")





############################################################
# B. QUANTILE NORMALIZED LOG2 DATA LOADING AND RESCALING.
############################################################
# Load data.
dat <- read.table("/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/AA/AA_Methyl_Seq/Analyses/20141208-From_George/20150507-AA_RNA_DASampleSpearmanCorr_Cytoscape/AA_DA_RNA_Reordered_QN.txt",  header =TRUE, sep = "\t")

# SAMPLE Numbers.
# GeneID	KNAA4A	KNAA4B	KNAA6A	KNAA6B	KNAA8B	KNAA8A	KN.C0001	KN.C0002	KN6	KN5	KN1075	KN1076	KN1892	KN1891	KN4002	KN4001	KN4905	KN4906	KN6115	KN6116	KN6506	KN6505	KN6857	KN6858	KN0007	KN0008	KN1293	KN1294	KN1495	KN1496	KN1517	KN1518	KN1610	KN1611	KN1711	KN1712	KN2027	KN2028	KN4149	KN4150	KN4749	KN4750	KN4783	KN4784	KN4817	KN4818	KN5235	KN5236	KN6335	KN6336	KN0019	KN0020	KN1033	KN1034	KN1803	KN1804	KN2195	KN2196	KN4117	KN4118	KN4175	KN4176	KN4473	KN4474	KN5261	KN5262	KN6035	KN6036	KN6201	KN6202	KN6295	KN6296	KN6709	KN6710


row.names(dat) <- dat$GeneID
dat_data <- as.matrix(dat[,2:length(dat[1,])])
nrow <- length(dat_data[,1])
ncol <- length(dat_data[1,])
datrownames <- row.names(dat)
datcolnames <- colnames(dat)[2:length(dat[1,])]

# Standardization by row.
dat_data_st <- dat_data
for (i in 1:length(dat_data_st[,1])) {
	dat_data_st[i,] <- (dat_data_st[i,] - mean(dat_data_st[i,]))/sqrt(var(dat_data_st[i,]))
}


############################################################
# C. SPEARMAN CORRELATION
############################################################
# Spearman Correlation.
scor <- cor(dat_data_st, use="everything", method="spearman")
# Output Spearman Correlation Coefficient.
write.table(scor, file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/AA/AA_Methyl_Seq/Analyses/20141208-From_George/20150507-AA_RNA_DASampleSpearmanCorr_Cytoscape/AA_DA_RNA_Spearman_Rho.txt", sep = "\t", row.names = TRUE, col.names = TRUE, append = FALSE)

scor2 <- scor^2
LDheatmap(scor2, genetic.distance=NULL, distances="physical", LDmeasure="r", title="AA RNA-Spearman Correlation", add.map=FALSE, add.key=TRUE, geneMapLocation=0.15, geneMapLabelX=0.5, geneMapLabelY=0.3, SNP.name=colnames(dat_data), color=rainbow(1000), newpage=TRUE, name="ldheatmap", vp.name=NULL, pop=FALSE)
# Output Spearman Correlation Coefficient Squared.
write.table(scor2, file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/AA/AA_Methyl_Seq/Analyses/20141208-From_George/20150507-AA_RNA_DASampleSpearmanCorr_Cytoscape/AA_DA_RNA_Spearman_RhoSquared.txt", sep = "\t", row.names = TRUE, col.names = TRUE, append = FALSE)




