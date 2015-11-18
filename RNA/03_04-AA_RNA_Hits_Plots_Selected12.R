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
# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


# In V2 the gene symbols are put first in the row names.
dat <- read.table("/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/AA/AA_Methyl_Seq/Analyses/20141208-From_George/20150317-RNA_Seq_EpigenomicMarks/AA_RNA_Hits_Reordered_QN_Selected12.txt", sep = "\t", header = TRUE)
# Column Names.
# GeneID	Transcript_ID_Chr_L_R	KNAA3B	KNAA3A	KNAA4A	KNAA4B	KNAA6A	KNAA6B	KNAA8B	KNAA8A	KN-C0001	KN-C0002	KN6	KN5	KN0016	KN0015	KN1075	KN1076	KN1892	KN1891	KN4002	KN4001	KN4905	KN4906	KN6115	KN6116	KN6506	KN6505	KN6857	KN6858	KN0007	KN0008	KN1293	KN1294	KN1495	KN1496	KN1517	KN1518	KN1610	KN1611	KN1711	KN1712	KN2027	KN2028	KN4091	KN4092	KN4149	KN4150	KN4749	KN4750	KN4783	KN4784	KN4817	KN4818	KN5235	KN5236	KN6335	KN6336	KN0019	KN0020	KN1033	KN1034	KN1803	KN1804	KN2195	KN2196	KN4117	KN4118	KN4175	KN4176	KN4473	KN4474	KN5261	KN5262	KN6035	KN6036	KN6201	KN6202	KN6295	KN6296	KN6709	KN6710


dat <- dat[,-2]
row.names(dat) <- make.names(dat$GeneID, unique = TRUE)
dat_data <- as.matrix(dat[,2:length(dat[1,])])
# # Clean the low end of dat_data by resetting all FPKM < 3 to 0.
# # Cutoff value is log2(3) = 1.584963.
# dat_data[dat_data < log2(3)] <- 0
nrow <- length(dat_data[,1])
ncol <- length(dat_data[1,])
datrownames <- row.names(dat)
datcolnames <- colnames(dat)[2:length(colnames(dat))]
Cat <- c("DISAA", "DISNM", "DISAA", "DISNM", "DISAA", "DISNM", "DISAA", "DISNM", "DISAA", "DISNM", "DISAA", "DISNM", "DISAA", "DISNM", "DISAA", "DISNM", "DISAA", "DISNM", "DISAA", "DISNM", "DISAA", "DISNM", "DISAA", "DISNM", "DISAA", "DISNM", "DISAA", "DISNM", "CONAA", "CONAA", "CONAA", "CONAA", "CONAA", "CONAA", "CONAA", "CONAA", "CONAA", "CONAA", "CONAA", "CONAA", "CONAA", "CONAA", "CONAA", "CONAA", "CONAA", "CONAA", "CONAA", "CONAA", "CONAA", "CONAA", "CONAA", "CONAA", "CONAA", "CONAA", "CONAA", "CONAA", "CONNM", "CONNM", "CONNM", "CONNM", "CONNM", "CONNM", "CONNM", "CONNM", "CONNM", "CONNM", "CONNM", "CONNM", "CONNM", "CONNM", "CONNM", "CONNM", "CONNM", "CONNM", "CONNM", "CONNM", "CONNM", "CONNM", "CONNM", "CONNM")
Cat2 <- c("AA", "NM", "AA", "NM", "AA", "NM", "AA", "NM", "AA", "NM", "AA", "NM", "AA", "NM", "AA", "NM", "AA", "NM", "AA", "NM", "AA", "NM", "AA", "NM", "AA", "NM", "AA", "NM", "AA", "AA", "AA", "AA", "AA", "AA", "AA", "AA", "AA", "AA", "AA", "AA", "AA", "AA", "AA", "AA", "AA", "AA", "AA", "AA", "AA", "AA", "AA", "AA", "AA", "AA", "AA", "AA", "NM", "NM", "NM", "NM", "NM", "NM", "NM", "NM", "NM", "NM", "NM", "NM", "NM", "NM", "NM", "NM", "NM", "NM", "NM", "NM", "NM", "NM", "NM", "NM")
Cat3 <- c("DISAA", "DISNM", "DISAA", "DISNM", "DISAA", "DISNM", "DISAA", "DISNM", "DISAA", "DISNM", "DISAA", "DISNM", "DISAA", "DISNM", "DISAA", "DISNM", "DISAA", "DISNM", "DISAA", "DISNM", "DISAA", "DISNM", "DISAA", "DISNM", "DISAA", "DISNM", "DISAA", "DISNM", "CON1", "CON2", "CON1", "CON2", "CON1", "CON2", "CON1", "CON2", "CON1", "CON2", "CON1", "CON2", "CON1", "CON2", "CON1", "CON2", "CON1", "CON2", "CON1", "CON2", "CON1", "CON2", "CON1", "CON2", "CON1", "CON2", "CON1", "CON2", "CON1", "CON2", "CON1", "CON2", "CON1", "CON2", "CON1", "CON2", "CON1", "CON2", "CON1", "CON2", "CON1", "CON2", "CON1", "CON2", "CON1", "CON2", "CON1", "CON2", "CON1", "CON2", "CON1", "CON2")


# Heatmap.
dat_data_da <- dat_data
dat_data_mean <- matrix(0,12,24)
for (j in 1:12) {
	for (i in 1:12) {
		dat_data_mean[j,i] <- mean(c(dat_data_da[j,i], dat_data_da[j,i+12]))
		dat_data_mean[j,i+12] <- mean(c(dat_data_da[j,i], dat_data_da[j,i+12]))
	}
}
dat_data_centered <- dat_data_da - dat_data_mean
my_palette <- colorRampPalette(c("turquoise", "gray12", "gold"))(n = 500)
hc_ts <- hclust(dist(t(dat_data_centered), method = "euclidean"), method = "ward.D")
total_heatmap <- heatmap.2(dat_data_centered, Rowv=FALSE, Colv=as.dendrogram(hc_ts), dendrogram = "column", col = my_palette, margins=c(5,18), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(-12,-0.1,length=167), seq(-0.1,0.1,length=167), seq(0.1,12,length=167)), ColSideColors = c(rep("red", 12), rep("green", 12)), keysize = 0.5, cexRow = 1)
total_heatmap <- heatmap.2(dat_data_centered, Colv=FALSE, Rowv=FALSE, col = my_palette, margins=c(5,18), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(-12,-0.1,length=167), seq(-0.1,0.1,length=167), seq(0.1,12,length=167)), ColSideColors = c(rep("red", 12), rep("green", 12)), keysize = 0.5, cexRow = 1)


# Boxplot of the groups.
# Four Groups.
par(mfrow = c(5,6))
pv <- c()
for (i in 1:30) {
	dat_row <- cbind(dat_data[i,], Cat)
	dat_row_df <- as.data.frame(dat_row)
	test <- kruskal.test(as.numeric(as.character(dat_row_df[,1])), dat_row_df[,2])
	pv[i] <- test$p.value
	pvprint <- round(pv[i], digits=3)
	boxplot(as.numeric(as.character(dat_row[,1]))~dat_row[,2], main = paste(datrownames[i], ", p = ", pvprint))
	stripchart(as.numeric(dat_row[,1])~dat_row[,2], vertical = TRUE, method = "jitter", pch = 16, col = "blue", bg = "bisque", add = TRUE) 
}
par(mfrow = c(5,6))
pv <- c()
for (i in 31:60) {
	dat_row <- cbind(dat_data[i,], Cat)
	dat_row_df <- as.data.frame(dat_row)
	test <- kruskal.test(as.numeric(as.character(dat_row_df[,1])), dat_row_df[,2])
	pv[i] <- test$p.value
	pvprint <- round(pv[i], digits=3)
	boxplot(as.numeric(as.character(dat_row[,1]))~dat_row[,2], main = paste(datrownames[i], ", p = ", pvprint))
	stripchart(as.numeric(dat_row[,1])~dat_row[,2], vertical = TRUE, method = "jitter", pch = 16, col = "blue", bg = "bisque", add = TRUE) 
}
par(mfrow = c(5,6))
pv <- c()
for (i in 61:90) {
	dat_row <- cbind(dat_data[i,], Cat)
	dat_row_df <- as.data.frame(dat_row)
	test <- kruskal.test(as.numeric(as.character(dat_row_df[,1])), dat_row_df[,2])
	pv[i] <- test$p.value
	pvprint <- round(pv[i], digits=3)
	boxplot(as.numeric(as.character(dat_row[,1]))~dat_row[,2], main = paste(datrownames[i], ", p = ", pvprint))
	stripchart(as.numeric(dat_row[,1])~dat_row[,2], vertical = TRUE, method = "jitter", pch = 16, col = "blue", bg = "bisque", add = TRUE) 
}
par(mfrow = c(5,6))
pv <- c()
for (i in 91:length(dat_data[,1])) {
	dat_row <- cbind(dat_data[i,], Cat)
	dat_row_df <- as.data.frame(dat_row)
	test <- kruskal.test(as.numeric(as.character(dat_row_df[,1])), dat_row_df[,2])
	pv[i] <- test$p.value
	pvprint <- round(pv[i], digits=3)
	boxplot(as.numeric(as.character(dat_row[,1]))~dat_row[,2], main = paste(datrownames[i], ", p = ", pvprint))
	stripchart(as.numeric(dat_row[,1])~dat_row[,2], vertical = TRUE, method = "jitter", pch = 16, col = "blue", bg = "bisque", add = TRUE) 
}
# Two Groups.
par(mfrow = c(5,6))
pv <- c()
for (i in 1:30) {
	dat_row <- cbind(dat_data[i,], Cat2)
	dat_row_df <- as.data.frame(dat_row)
	dat.split <- split(dat_row_df, dat_row_df[,2])
	t1 <- dat.split$AA[,1]
	t2 <- dat.split$NM[,1]
	test <- wilcox.test(as.numeric(as.character(t1)), as.numeric(as.character(t2)), paired = FALSE, alternative = "two.sided")
	pv[i] <- test$p.value
	pvprint <- round(pv[i], digits=3)
	boxplot(as.numeric(dat_row[,1])~dat_row[,2], main = paste(datrownames[i], ", p = ", pvprint))
	if (mean(as.numeric(as.character(t1))) < mean(as.numeric(as.character(t2)))) {
		colr <- c("turquoise")
	} else {
		colr <- c("gold")
	}
	stripchart(as.numeric(dat_row[,1])~dat_row[,2], vertical = TRUE, method = "jitter", pch = 16, col = colr, bg = "bisque", add = TRUE) 
}
par(mfrow = c(5,6))
pv <- c()
for (i in 31:60) {
	dat_row <- cbind(dat_data[i,], Cat2)
	dat_row_df <- as.data.frame(dat_row)
	dat.split <- split(dat_row_df, dat_row_df[,2])
	t1 <- dat.split$AA[,1]
	t2 <- dat.split$NM[,1]
	test <- wilcox.test(as.numeric(as.character(t1)), as.numeric(as.character(t2)), paired = FALSE, alternative = "two.sided")
	pv[i] <- test$p.value
	pvprint <- round(pv[i], digits=3)
	boxplot(as.numeric(dat_row[,1])~dat_row[,2], main = paste(datrownames[i], ", p = ", pvprint))
	if (mean(as.numeric(as.character(t1))) < mean(as.numeric(as.character(t2)))) {
		colr <- c("turquoise")
	} else {
		colr <- c("gold")
	}
	stripchart(as.numeric(dat_row[,1])~dat_row[,2], vertical = TRUE, method = "jitter", pch = 16, col = colr, bg = "bisque", add = TRUE) 
}
par(mfrow = c(5,6))
pv <- c()
for (i in 61:90) {
	dat_row <- cbind(dat_data[i,], Cat2)
	dat_row_df <- as.data.frame(dat_row)
	dat.split <- split(dat_row_df, dat_row_df[,2])
	t1 <- dat.split$AA[,1]
	t2 <- dat.split$NM[,1]
	test <- wilcox.test(as.numeric(as.character(t1)), as.numeric(as.character(t2)), paired = FALSE, alternative = "two.sided")
	pv[i] <- test$p.value
	pvprint <- round(pv[i], digits=3)
	boxplot(as.numeric(dat_row[,1])~dat_row[,2], main = paste(datrownames[i], ", p = ", pvprint))
	if (mean(as.numeric(as.character(t1))) < mean(as.numeric(as.character(t2)))) {
		colr <- c("turquoise")
	} else {
		colr <- c("gold")
	}
	stripchart(as.numeric(dat_row[,1])~dat_row[,2], vertical = TRUE, method = "jitter", pch = 16, col = colr, bg = "bisque", add = TRUE) 
}
par(mfrow = c(5,6))
pv <- c()
for (i in 91:length(dat_data[,1])) {
	dat_row <- cbind(dat_data[i,], Cat2)
	dat_row_df <- as.data.frame(dat_row)
	dat.split <- split(dat_row_df, dat_row_df[,2])
	t1 <- dat.split$AA[,1]
	t2 <- dat.split$NM[,1]
	test <- wilcox.test(as.numeric(as.character(t1)), as.numeric(as.character(t2)), paired = FALSE, alternative = "two.sided")
	pv[i] <- test$p.value
	pvprint <- round(pv[i], digits=3)
	boxplot(as.numeric(dat_row[,1])~dat_row[,2], main = paste(datrownames[i], ", p = ", pvprint))
	if (mean(as.numeric(as.character(t1))) < mean(as.numeric(as.character(t2)))) {
		colr <- c("turquoise")
	} else {
		colr <- c("gold")
	}
	stripchart(as.numeric(dat_row[,1])~dat_row[,2], vertical = TRUE, method = "jitter", pch = 16, col = colr, bg = "bisque", add = TRUE) 
}


# UNPaired Test in Discordant Twins.
# Twelve Analytes -- T-Test.
par(mfrow = c(3,4))
pv <- c()
for (i in 1:length(dat_data[,1])) {
	t1 <- dat_data[i,c(1:14, 29:56)]
	t2 <- dat_data[i,c(15:28, 57:80)]
	test <- t.test(as.numeric(as.character(t1)), as.numeric(as.character(t2)), paired = FALSE, alternative = "two.sided")
	pv[i] <- round(test$p.value, digit = 3)
	tplot <- cbind(c(t1, t2), c(rep("Asthma", 12), rep("Healthy", 12)))
	boxplot(as.numeric(as.character(tplot[,1]))~as.character(tplot[,2]), main = paste(datrownames[i], ", p = ", pv[i]))
	if (mean(as.numeric(as.character(t1))) < mean(as.numeric(as.character(t2)))) {
		colr <- c("turquoise")
	} else {
		colr <- c("gold")
	}
	stripchart(as.numeric(as.character(tplot[,1]))~as.character(tplot[,2]), vertical = TRUE, method = "overplot", pch = 16, col = colr, bg = "bisque", add = TRUE) 
}
# Twelve Analytes -- U-Test.
par(mfrow = c(3,4))
pv <- c()
for (i in 1:length(dat_data[,1])) {
	t1 <- dat_data[i,c(1:14, 29:56)]
	t2 <- dat_data[i,c(15:28, 57:80)]
	test <- wilcox.test(as.numeric(as.character(t1)), as.numeric(as.character(t2)), paired = FALSE, alternative = "two.sided")
	pv[i] <- round(test$p.value, digit = 3)
	tplot <- cbind(c(t1, t2), c(rep("Asthma", 12), rep("Healthy", 12)))
	boxplot(as.numeric(as.character(tplot[,1]))~as.character(tplot[,2]), main = paste(datrownames[i], ", p = ", pv[i]))
	if (mean(as.numeric(as.character(t1))) < mean(as.numeric(as.character(t2)))) {
		colr <- c("turquoise")
	} else {
		colr <- c("gold")
	}
	stripchart(as.numeric(as.character(tplot[,1]))~as.character(tplot[,2]), vertical = TRUE, method = "overplot", pch = 16, col = colr, bg = "bisque", add = TRUE) 
}




