############################################################
# TO-DO LIST.
############################################################
DONE. # 1. Analyze AA MS Proteome Discordant Data as of 12/10/2015.





############################################################
# A. Library Used.
############################################################
# Libraries used.
library("gplots")


# Libraries not used.
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
dat1 <- read.table("/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/AA/AA_Proteomics_20150729/Results/20151206_1stRepHitsBasedonHannes/AA_Allreplicate1_reloaded_12092015_genecbd_p0_01_daonly.txt",  header =TRUE, sep = "\t")
dat2 <- read.table("/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/AA/AA_Proteomics_20150729/Results/20151206_1stRepHitsBasedonHannes/AA_Allreplicate1_reloaded_12092015_genecbd_p0_01_full.txt",  header =TRUE, sep = "\t")

# Gene Name	KN0006	KN0005	KN0008	KN0007	KN0016	KN0015	KN0020	KN0019	REF1R1	KN1034	KN1033	KN1076	KN1075	KN1294	KN1293	KN1496	KN1495	REF2R1	KN1518	KN1517	KN1611	KN1610	KN1804	KN1803	KN1892	KN1891	REF3R1	KN2024	KN2023	KN2028	KN2027	KN2196	KN2195	KN3740	KN3739	REF4R1	KN4092	KN4091	KN4118	KN4117	KN4150	KN4149	KN4176	KN4175	REF5R1	KN4750	KN4749	KN4784	KN4783	KN4818	KN4817	KN4906	KN4905	REF6R1	KN4964	KN4963	KN5236	KN5235	KN5274	KN5273	KN6036	KN6035	REF7R1	KN6116	KN6115	KN6202	KN6201	KN6296	KN6295	KN6506	KN6505	REF8R1	KN6709	KN6710	KN6857	KN6858	REF9R1
# Gene Name	DA_A	DA_N	CA_A	CA_A	DA_A	DA_N	CH_N	CH_N	REF1R1	CH_N	CH_N	DA_N	DA_A	CA_A	CA_A	CA_A	CA_A	REF2R1	CA_A	CA_A	CA_A	CA_A	CH_N	CH_N	DA_A	DA_N	REF3R1	CA_A	CA_A	CA_A	CA_A	CH_N	CH_N	CA_A	CA_A	REF4R1	CA_A	CA_A	CH_N	CH_N	CA_A	CA_A	CH_N	CH_N	REF5R1	CA_A	CA_A	CA_A	CA_A	CA_A	CA_A	DA_N	DA_A	REF6R1	CA_A	CA_A	CA_A	CA_A	CH_N	CH_N	CH_N	CH_N	REF7R1	DA_N	DA_A	CH_N	CH_N	CH_N	CH_N	DA_A	DA_N	REF8R1	CH_N	CH_N	DA_A	DA_N	REF9R1

status <- dat1[1,]
dat1 <- dat1[-1,]
dat2 <- dat2[-1,]

row.names(dat1) <- dat1[,1]
dat_data1 <- as.matrix(dat1[,2:length(dat1[1,])])
class(dat_data1) <- "numeric"
dim(dat_data1)
[1] 135  77

row.names(dat2) <- dat2[,1]
dat_data2 <- as.matrix(dat2[,2:length(dat2[1,])])
class(dat_data2) <- "numeric"
dim(dat_data2)
[1] 151  77

# Define groups.
da_a <- c(1, 5, 13, 25, 53, 65, 70, 75)
da_n <- c(2, 6, 12, 26, 52, 64, 71, 76)
ca_a <- c(3, 4, 14, 15, 16, 17, 19, 20, 21, 22, 28, 29, 30, 31, 34, 35, 37, 38, 41, 42, 46, 47, 48, 49, 50, 51, 55, 56, 57, 58)
ch_n <- c(7, 8, 10, 11, 23, 24, 32, 33, 39, 40, 43, 44, 59, 60, 61, 62, 66, 67, 68, 69, 73, 74)
ref <- c(9, 18, 27, 36, 45, 54, 63, 72, 77)


# Heatmap of Raw Data (DA Hits).
my_palette <- colorRampPalette(c("turquoise", "gray12", "gold"))(n = 500)
# Row Cluster.
hc_ts <- hclust(dist(dat_data1, method = "euclidean"), method = "ward.D")
# Row Scaled.
total_heatmap <- heatmap.2(dat_data1[,c(da_a,da_n,ca_a,ch_n)], Rowv=as.dendrogram(hc_ts), Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(5,7), density.info = "none", trace = "none", scale = "row", symkey = FALSE, keysize = 0.5, cexRow = 0.8, ColSideColor = c(rep("gold", 8), rep("turquoise", 8), rep("red", 30), rep("black", 22)))
# Unscaled.
total_heatmap <- heatmap.2(dat_data1[,c(da_a,da_n,ca_a,ch_n)], Rowv=as.dendrogram(hc_ts), Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(5,7), density.info = "none", trace = "none", scale = "none", symkey = FALSE, keysize = 0.5, cexRow = 0.8, ColSideColor = c(rep("gold", 8), rep("turquoise", 8), rep("red", 30), rep("black", 22)), breaks = c(seq(0,0.5,length=167), seq(0.51,2,length=167), seq(2.01,15,length=167)))
# Row + Col Cluster.
hc_ts2 <- hclust(dist(t(dat_data1[,c(da_a,da_n,ca_a,ch_n)]), method = "euclidean"), method = "ward.D")
# Row Scaled.
total_heatmap <- heatmap.2(dat_data1[,c(da_a,da_n,ca_a,ch_n)], Rowv=as.dendrogram(hc_ts), Colv=as.dendrogram(hc_ts2), dendrogram = "both", col = my_palette, margins=c(5,7), density.info = "none", trace = "none", scale = "row", symkey = FALSE, keysize = 0.5, cexRow = 0.8, ColSideColor = c(rep("gold", 8), rep("turquoise", 8), rep("red", 30), rep("black", 22)))
# Unscaled.
total_heatmap <- heatmap.2(dat_data1[,c(da_a,da_n,ca_a,ch_n)], Rowv=as.dendrogram(hc_ts), Colv=as.dendrogram(hc_ts2), dendrogram = "both", col = my_palette, margins=c(5,7), density.info = "none", trace = "none", scale = "none", symkey = FALSE, keysize = 0.5, cexRow = 0.8, ColSideColor = c(rep("gold", 8), rep("turquoise", 8), rep("red", 30), rep("black", 22)), breaks = c(seq(0,0.5,length=167), seq(0.51,2,length=167), seq(2.01,15,length=167)))
# Diff.
dat_tmp <- dat_data1[,c(da_a,da_n)]
dat_tmp[,c(1:8)] <- (dat_data1[,c(da_a)] - dat_data1[,c(da_n)]) / 2
dat_tmp[,c(9:16)] <- -(dat_data1[,c(da_a)] - dat_data1[,c(da_n)]) / 2
dat_tmp_da_a_mean <- apply(dat_tmp[,c(1:8)], 1, mean)
n <- length(dat_tmp[,1])
dat_tmp <- dat_tmp[(1:n)[order(dat_tmp_da_a_mean)],]
hc_ts <- hclust(dist(dat_tmp, method = "euclidean"), method = "ward.D")
# Row Scaled.
total_heatmap <- heatmap.2(dat_tmp, Rowv=as.dendrogram(hc_ts), Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(7,7), density.info = "none", trace = "none", scale = "row", symkey = FALSE, keysize = 0.5, cexRow = 0.8, ColSideColor = c(rep("gold", 8), rep("turquoise", 8)))
# Unscaled.
total_heatmap <- heatmap.2(dat_tmp, Rowv=as.dendrogram(hc_ts), Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(7,7), density.info = "none", trace = "none", scale = "none", symkey = FALSE, keysize = 0.5, cexRow = 0.8, ColSideColor = c(rep("gold", 8), rep("turquoise", 8)), breaks = c(seq(-5,-0.5,length=167), seq(-0.49,0.49,length=167), seq(0.5,5,length=167)))


# Heatmap of Raw Data (Full Sample Set Hits).
my_palette <- colorRampPalette(c("turquoise", "gray12", "gold"))(n = 500)
# Row Cluster.
hc_ts <- hclust(dist(dat_data2, method = "euclidean"), method = "ward.D")
# Row Scaled.
total_heatmap <- heatmap.2(dat_data2[,c(da_a,da_n,ca_a,ch_n)], Rowv=as.dendrogram(hc_ts), Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(5,7), density.info = "none", trace = "none", scale = "row", symkey = FALSE, keysize = 0.5, cexRow = 0.8, ColSideColor = c(rep("gold", 8), rep("turquoise", 8), rep("red", 30), rep("black", 22)))
# Unscaled.
total_heatmap <- heatmap.2(dat_data2[,c(da_a,da_n,ca_a,ch_n)], Rowv=as.dendrogram(hc_ts), Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(5,7), density.info = "none", trace = "none", scale = "none", symkey = FALSE, keysize = 0.5, cexRow = 0.8, ColSideColor = c(rep("gold", 8), rep("turquoise", 8), rep("red", 30), rep("black", 22)), breaks = c(seq(0,0.5,length=167), seq(0.51,2,length=167), seq(2.01,15,length=167)))
# Row + Col Cluster.
hc_ts2 <- hclust(dist(t(dat_data2[,c(da_a,da_n,ca_a,ch_n)]), method = "euclidean"), method = "ward.D")
# Row Scaled.
total_heatmap <- heatmap.2(dat_data2[,c(da_a,da_n,ca_a,ch_n)], Rowv=as.dendrogram(hc_ts), Colv=as.dendrogram(hc_ts2), dendrogram = "both", col = my_palette, margins=c(5,7), density.info = "none", trace = "none", scale = "row", symkey = FALSE, keysize = 0.5, cexRow = 0.8, ColSideColor = c(rep("gold", 8), rep("turquoise", 8), rep("red", 30), rep("black", 22)))
# Unscaled.
total_heatmap <- heatmap.2(dat_data2[,c(da_a,da_n,ca_a,ch_n)], Rowv=as.dendrogram(hc_ts), Colv=as.dendrogram(hc_ts2), dendrogram = "both", col = my_palette, margins=c(5,7), density.info = "none", trace = "none", scale = "none", symkey = FALSE, keysize = 0.5, cexRow = 0.8, ColSideColor = c(rep("gold", 8), rep("turquoise", 8), rep("red", 30), rep("black", 22)), breaks = c(seq(0,0.5,length=167), seq(0.51,2,length=167), seq(2.01,15,length=167)))
# Diff.
dat_tmp <- dat_data2[,c(da_a,da_n)]
dat_tmp[,c(1:8)] <- (dat_data2[,c(da_a)] - dat_data2[,c(da_n)]) / 2
dat_tmp[,c(9:16)] <- -(dat_data2[,c(da_a)] - dat_data2[,c(da_n)]) / 2
dat_tmp_da_a_mean <- apply(dat_tmp[,c(1:8)], 1, mean)
n <- length(dat_tmp[,1])
dat_tmp <- dat_tmp[(1:n)[order(dat_tmp_da_a_mean)],]
hc_ts <- hclust(dist(dat_tmp, method = "euclidean"), method = "ward.D")
# Row Scaled.
total_heatmap <- heatmap.2(dat_tmp, Rowv=as.dendrogram(hc_ts), Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(7,7), density.info = "none", trace = "none", scale = "row", symkey = FALSE, keysize = 0.5, cexRow = 0.8, ColSideColor = c(rep("gold", 8), rep("turquoise", 8)))
# Unscaled.
total_heatmap <- heatmap.2(dat_tmp, Rowv=as.dendrogram(hc_ts), Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(7,7), density.info = "none", trace = "none", scale = "none", symkey = FALSE, keysize = 0.5, cexRow = 0.8, ColSideColor = c(rep("gold", 8), rep("turquoise", 8)), breaks = c(seq(-5,-0.5,length=167), seq(-0.49,0.49,length=167), seq(0.5,5,length=167)))




