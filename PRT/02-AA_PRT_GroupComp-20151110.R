############################################################
# TO-DO LIST.
############################################################
DONE. # 1. Analyze AA MS Proteome Discordant Data as of 11/10/2015.





############################################################
# A. Library Used.
############################################################
library("gplots")
library("ggplot2")
library("preprocessCore")
library("coin")
library("magic")
library("RColorBrewer")
library("limma")
library("LDheatmap")
library("devtools")
# install_github("ggbiplot", "vqv")
library("ggbiplot")
library("plyr")
library("reshape")
library("MASS")
library("plotrix")
library("seriation")

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
  library(grid)

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

# PLOTCOR MODIFIED.
plotcormod <- function(r, addtext=TRUE, atcex=NULL, incdiag=FALSE,
  rorder=FALSE, ...) {

  # round to the nearest hundredth
	rr <- round(r, 2)

	dimr <- dim(r)
	sym <- isSymmetric(r)

	# get rid of diagonal numbers
	if (!incdiag & sym) {
    diag(rr) <- NA
	}

	rrf <- format(rr)
	rrf[grep("NA", rrf)] <- ""
	rra <- abs(rr)
	nx <- dimr[2]
	ny <- dimr[1]
	if (is.null(atcex)) {
    atcex <- 8/max(nx, ny)
	}
	namzx <- dimnames(rr)[[2]]
	namzy <- dimnames(rr)[[1]]

	# order rows/columns
	ordx <- 1:nx
	ordy <- 1:ny
	if (rorder) {
		# the seriate() function prints out % explained variance for method="PCA"
		# I used capture.output to avoid having this print to the screen
		dummy <- capture.output(ser <- seriate((1-r)/2, method="PCA"))
		ordy <- rev(get_order(ser, 1))
		ordx <- rev(get_order(ser, 2))
	}
	if (sym) {
    ordx <- rev(ordy)
	}

	# categorize correlations from -1 to 1 by 0.01
	brks <- seq(-1, 1, 0.01)
	rcat <- apply(rr, 2, cut, breaks=brks, labels=FALSE)

	# ORI: assign colors on the cyan-magenta scale
	# assign colors on the turquoise-gray12-gold scale
	# my_palette <- colorRampPalette(c("turquoise", "gray12", "gold"))(n = length(brks))
	my_palette <- colorRampPalette(c("turquoise", "turquoise", "gold", "gold"))(n = length(brks))
	colz <- apply(rcat, 2, function(x) my_palette[x])
	par(xaxs="i", yaxs="i", mar=c(0.1, 6, 6, 0.1), ...)
	eqscplot(1, 1, type="n", xlim=c(0.5, nx+0.5), ylim=c(0.5, ny+0.5),
    xlab="", ylab="", axes=FALSE, tol = 0)
	for(i in 1:nx) {
	for(j in 1:ny) {
		io <- ordx[i]
		jo <- ordy[j]
		draw.ellipse(i, j, a=0.5, b=0.5, col=colz[jo, io], border=NA)
		draw.ellipse(i, j, a=0.5, b=(1-rr[jo, io]^2)/2,
      angle=45*c(-1, 1)[1+(rr[jo, io]>0)], col="white", border=NA)
		if (addtext & !is.na(rra[jo, io])) {
      text(i, j, rrf[jo, io], cex=atcex, col=rgb(0, 0, 0, alpha=rra[jo, io]))
		}
	}}
	axis(3, at=1:nx, labels=namzx[ordx], las=2, tick=FALSE, cex=0.5)
	axis(2, at=1:ny, labels=namzy[ordy], las=2, tick=FALSE, cex=0.5)

	list(rev(ordy), ordx)
}





############################################################
# B. U TEST HITS UPDATED WITH UNIFIED CUT OFF (P < 0.005 OR FDR-P < 0.05).
############################################################
############################################################
# B.1. Raw Data Loading (All Quantiles).
############################################################
# Load data.
dat <- read.table("/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/AA/AA_Proteomics_20150729/Results/20151106_Pass02/AA_PRT_Combined_Results_ANNO_DA_20151110.txt",  header =TRUE, sep = "\t")

# Accession	KN0006	KN0005	KN0016	KN0015	KN1076	KN1075	KN1892	KN1891	KN4906	KN4905	KN6116	KN6115	KN6506	KN6505	KN6116_2	KN6115_2	KN6506_2	KN6505_2	REF1R1	REF2R1	REF3R1	REF6R1	REF8R1_1	REF8R1_2
# Diagnosis	DA_A	DA_N	DA_A	DA_N	DA_N	DA_A	DA_A	DA_N	DA_N	DA_A	DA_N	DA_A	DA_A	DA_N	DA_N	DA_A	DA_A	DA_N	REF1R1	REF2R1	REF3R1	REF6R1	REF8R1_1	REF8R1_2

row.names(dat) <- dat$Accession
dat_data <- as.matrix(dat[2:length(dat[,1]),2:length(dat[1,])])
dat_data <- apply(dat_data, 2, as.numeric)
dat_data[is.na(dat_data)] <- 0
nrow <- length(dat_data[,1])
ncol <- length(dat_data[1,])
datrownames <- row.names(dat)[2:length(dat[,1])]
row.names(dat_data) <- datrownames
datcolnames <- colnames(dat)[2:length(dat[1,])]
dim(dat_data)
[1] 7853   24

# Heatmap of Raw Data.
my_palette <- colorRampPalette(c("turquoise", "gray12", "gold"))(n = 500)
hc_ts <- hclust(dist(dat_data, method = "euclidean"), method = "ward.D")
# total_heatmap <- heatmap.2(dat_data, Rowv=as.dendrogram(hc_ts), Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(7,7), density.info = "none", trace = "none", scale = "none", symkey = FALSE, keysize = 0.5, cexRow = 0.6)
total_heatmap <- heatmap.2(dat_data, Rowv=as.dendrogram(hc_ts), Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(7,7), density.info = "none", trace = "none", scale = "row", symkey = FALSE, keysize = 0.5, cexRow = 0.6)
hc_ts2 <- hclust(dist(t(dat_data), method = "euclidean"), method = "ward.D")
# total_heatmap <- heatmap.2(dat_data, Colv=as.dendrogram(hc_ts2), Rowv=as.dendrogram(hc_ts), dendrogram = "both", col = my_palette, margins=c(7,7), density.info = "none", trace = "none", scale = "none", symkey = FALSE, keysize = 0.5, cexRow = 0.6)
total_heatmap <- heatmap.2(dat_data, Colv=as.dendrogram(hc_ts2), Rowv=as.dendrogram(hc_ts), dendrogram = "both", col = my_palette, margins=c(7,7), density.info = "none", trace = "none", scale = "row", symkey = FALSE, keysize = 0.5, cexRow = 0.6)


############################################################
# B.2. SPEARMAN CORRELATION
############################################################
# Correlation between samples.
scor <- cor(dat_data, method = "spearman")
plotcormod(scor, mar=c(2, 10, 10, 2))

# Correlation between in-twin differences (A-N).
andiff <- matrix(nrow = 7853, ncol = 7)
andiff[,1] <- dat_data[,1] - dat_data[,2]
andiff[,2] <- dat_data[,3] - dat_data[,4]
andiff[,3] <- dat_data[,6] - dat_data[,5]
andiff[,4] <- dat_data[,7] - dat_data[,8]
andiff[,5] <- dat_data[,10] - dat_data[,9]
andiff[,6] <- dat_data[,12] - dat_data[,11]
andiff[,7] <- dat_data[,13] - dat_data[,14]
colnames(andiff) <- c("KN0006-KN0005", "KN0016-KN0015", "KN1075-KN1076", "KN1892-KN1891", "KN4905-KN4906", "KN6115-KN6116", "KN6506-KN6505")
scordiff <- cor(andiff, method = "spearman")
plotcormod(scordiff, mar=c(2, 10, 10, 2))









############################################################
# B.3. U TEST
############################################################
# Use average ratio for the replicated samples.
dat_data_ave <- dat_data
dat_data_ave[,5] <- (dat_data[,5] + dat_data[,6]) / 2
dat_data_ave[,6] <- (dat_data[,5] + dat_data[,6]) / 2
dat_data_ave[,7] <- (dat_data[,7] + dat_data[,8]) / 2
dat_data_ave[,8] <- (dat_data[,7] + dat_data[,8]) / 2
dat_data_ave[,13] <- (dat_data[,13] + dat_data[,14]) / 2
dat_data_ave[,14] <- (dat_data[,13] + dat_data[,14]) / 2
dat_data_ave[,15] <- (dat_data[,15] + dat_data[,16]) / 2
dat_data_ave[,16] <- (dat_data[,15] + dat_data[,16]) / 2
n <- length(dat_data_ave[,1])

pv <- c()
for (i in 1:n)
   {
    	t1 <- c(dat_data[i,c(1:5,7)])
    	t2 <- c(dat_data[i,c(9:13,15)])
		test <- wilcox.test(t1, t2, paired = TRUE, alternative = "two.sided")
		pv[i] <- test$p.value
   }
pv[is.na(pv)] <- 1   

summary(pv)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
0.03125 0.21880 0.56250 0.52240 0.84380 1.00000 

# FDR-Adjusted p-value.
new_p <- p.adjust(pv, method = "fdr", n = length(pv))
summary(new_p)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
 0.8004  0.8190  0.9236  0.9047  0.9827  1.0000 

# Export Top Hits.
length((1:n)[pv < 0.05])
[1] 54
length((1:n)[new_p < 0.05])
[1] 0
# Export hits with p-value < 0.05.
writecontentall <- cbind(pv[(1:n)[pv < 0.05]], dat_data[(1:n)[pv < 0.05],])
write(c("Accession","p-value",colnames(dat_data)), file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/AA/AA_Proteomics_20150729/Results/20151029_Pass01/AA_PRT_Combined_Results_NA_p0_05.txt", sep = "\t", ncolumn = 23)
write.table(writecontentall, file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/AA/AA_Proteomics_20150729/Results/20151029_Pass01/AA_PRT_Combined_Results_NA_p0_05.txt", sep = "\t", row.names = TRUE, col.names = FALSE, append = TRUE, quote = FALSE)
# Heatmap of Hits.
my_palette <- colorRampPalette(c("turquoise", "gray12", "gold"))(n = 500)
hc_ts <- hclust(dist(writecontentall[,2:(length(writecontentall[1,]))], method = "euclidean"), method = "ward.D")
hc_ts2 <- hclust(dist(t(writecontentall[,2:(length(writecontentall[1,]))]), method = "euclidean"), method = "ward.D")
total_heatmap <- heatmap.2(writecontentall[,2:(length(writecontentall[1,]))], Rowv=as.dendrogram(hc_ts), Colv=as.dendrogram(hc_ts2), dendrogram = "both", col = my_palette, margins=c(7,7), density.info = "none", trace = "none", scale = "row", symkey = FALSE, keysize = 0.5, cexRow = 0.6, ColSideColors = c(rep("red",8), rep("blue", 8), rep("black", 5)))


############################################################
# B.4. DENSITY PLOTS
############################################################
# Density Plots of Samples by Pair.
densityplot <- c(list())
for (i in 1:8) {
	allsignals <- c(as.vector(dat_data[,i]), as.vector(dat_data[,i+8]))
	dendat <- data.frame(RATIO_TO_REF = allsignals, CATEGORIES = c(rep("ASTHMA", length(as.vector(dat_data[,i]))), rep("NORMAL", length(as.vector(dat_data[,i+8])))))
	densityplot[[i]] <- ggplot(dendat, aes(x = RATIO_TO_REF, fill = CATEGORIES)) + geom_density(alpha = 0.5) + ggtitle(paste("COMPARISON OF TWIN PROTEIN RATIOS, ASTHMA VS NORMAL\n(",colnames(dat_data)[i]," vs ", colnames(dat_data)[i+8], ")"))
}
multiplot(densityplot[[1]], densityplot[[2]], densityplot[[3]], densityplot[[4]], densityplot[[5]], densityplot[[6]], densityplot[[7]], densityplot[[8]], cols=2)

# Density Plots of Refs.
allsignals <- c(as.vector(dat_data[,17]), as.vector(dat_data[,18]),as.vector(dat_data[,19]), as.vector(dat_data[,20]), as.vector(dat_data[,21]))
dendat <- data.frame(RATIO_TO_REF126 = allsignals, CATEGORIES = c(rep("REF1R1", length(as.vector(dat_data[,17]))), rep("REF2R1", length(as.vector(dat_data[,18]))), rep("REF6R1", length(as.vector(dat_data[,19]))), rep("REF8R1_1", length(as.vector(dat_data[,20]))), rep("REF8R1_2", length(as.vector(dat_data[,21])))))
ggplot(dendat, aes(x = RATIO_TO_REF126, fill = CATEGORIES)) + geom_density(alpha = 0.5) + ggtitle(paste("TWIN PROTEIN REF131/126 RATIOS (5 REPLICATES)"))

# Density Plots of Samples by Pair (Log10).
densityplot <- c(list())
for (i in 1:8) {
	allsignals <- c(as.vector(log10(dat_data[,i])), as.vector(log10(dat_data[,i+8])))
	dendat <- data.frame(LOG10_RATIO_TO_REF = allsignals, CATEGORIES = c(rep("ASTHMA", length(as.vector(dat_data[,i]))), rep("NORMAL", length(as.vector(dat_data[,i+8])))))
	densityplot[[i]] <- ggplot(dendat, aes(x = LOG10_RATIO_TO_REF, fill = CATEGORIES)) + geom_density(alpha = 0.5) + ggtitle(paste("COMPARISON OF TWIN PROTEIN LOG10 RATIOS, ASTHMA VS NORMAL\n(",colnames(dat_data)[i]," vs ", colnames(dat_data)[i+8], ")"))
}
multiplot(densityplot[[1]], densityplot[[2]], densityplot[[3]], densityplot[[4]], densityplot[[5]], densityplot[[6]], densityplot[[7]], densityplot[[8]], cols=2)

# Density Plots of Refs (Log10).
allsignals <- c(as.vector(log10(dat_data[,17])), as.vector(log10(dat_data[,18])),as.vector(log10(dat_data[,19])), as.vector(log10(dat_data[,20])), as.vector(log10(dat_data[,21])))
dendat <- data.frame(LOG10_RATIO_TO_REF126 = allsignals, CATEGORIES = c(rep("REF1R1", length(as.vector(dat_data[,17]))), rep("REF2R1", length(as.vector(dat_data[,18]))), rep("REF6R1", length(as.vector(dat_data[,19]))), rep("REF8R1_1", length(as.vector(dat_data[,20]))), rep("REF8R1_2", length(as.vector(dat_data[,21])))))
ggplot(dendat, aes(x = LOG10_RATIO_TO_REF126, fill = CATEGORIES)) + geom_density(alpha = 0.5) + ggtitle(paste("TWIN PROTEIN REF131/126 LOG10 RATIOS (5 REPLICATES)"))




