library("gplots")
library("preprocessCore")
library("coin")
library("magic")
library("RColorBrewer")
library("limma")
library("LDheatmap")
library(affy)





dat <- read.table("/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/AA/AA_RNA_Seq_20140630/20140630-NewestRound/Analyses/Combined/AA37TxExp_Reordered_QN.txt", header = TRUE, sep = "\t")
dat.m <- as.matrix(dat[,2:length(dat[1,])])
row.names(dat.m) <- dat[,1]
col.name <- colnames (dat.m)

plotname <- c()
par(mfrow = c(5,8), mar = c(2,2,2,2))
for (i in seq(1, 74, by = 2)) {
	plotname[i] <- paste(col.name[i], "vs", col.name[i+1], sep = " ", collapse = NULL)
	M <- dat.m[,i] - dat.m[,i+1]
	A <- (dat.m[,i] + dat.m[,i+1])/2
	ma.plot(A, M, cex=1, plot.method="smoothScatter") 
	title(plotname[i])
}




