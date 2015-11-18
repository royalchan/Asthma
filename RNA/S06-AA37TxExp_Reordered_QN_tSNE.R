library("tsne")
library("scatterplot3d")





# Category Info.
Cat <- c("DISAA", "DISNM", "DISAA", "DISNM", "DISAA", "DISNM", "DISAA", "DISNM", "DISAA", "DISNM", "DISAA", "DISNM", "DISAA", "DISNM", "DISAA", "DISNM", "DISAA", "DISNM", "DISAA", "DISNM", "DISAA", "DISNM", "DISAA", "DISNM", "CONAA", "CONAA", "CONAA", "CONAA", "CONAA", "CONAA", "CONAA", "CONAA", "CONAA", "CONAA", "CONAA", "CONAA", "CONAA", "CONAA", "CONAA", "CONAA", "CONAA", "CONAA", "CONAA", "CONAA", "CONAA", "CONAA", "CONAA", "CONAA", "CONAA", "CONAA", "CONNM", "CONNM", "CONNM", "CONNM", "CONNM", "CONNM", "CONNM", "CONNM", "CONNM", "CONNM", "CONNM", "CONNM", "CONNM", "CONNM", "CONNM", "CONNM", "CONNM", "CONNM", "CONNM", "CONNM", "CONNM", "CONNM", "CONNM", "CONNM")
Cat2 <- c("AA", "NM", "AA", "NM", "AA", "NM", "AA", "NM", "AA", "NM", "AA", "NM", "AA", "NM", "AA", "NM", "AA", "NM", "AA", "NM", "AA", "NM", "AA", "NM", "AA", "AA", "AA", "AA", "AA", "AA", "AA", "AA", "AA", "AA", "AA", "AA", "AA", "AA", "AA", "AA", "AA", "AA", "AA", "AA", "AA", "AA", "AA", "AA", "AA", "AA", "NM", "NM", "NM", "NM", "NM", "NM", "NM", "NM", "NM", "NM", "NM", "NM", "NM", "NM", "NM", "NM", "NM", "NM", "NM", "NM", "NM", "NM", "NM", "NM")
Cat3 <- c("DISAA", "DISNM", "DISAA", "DISNM", "DISAA", "DISNM", "DISAA", "DISNM", "DISAA", "DISNM", "DISAA", "DISNM", "DISAA", "DISNM", "DISAA", "DISNM", "DISAA", "DISNM", "DISAA", "DISNM", "DISAA", "DISNM", "DISAA", "DISNM", "CON1", "CON2", "CON1", "CON2", "CON1", "CON2", "CON1", "CON2", "CON1", "CON2", "CON1", "CON2", "CON1", "CON2", "CON1", "CON2", "CON1", "CON2", "CON1", "CON2", "CON1", "CON2", "CON1", "CON2", "CON1", "CON2", "CON1", "CON2", "CON1", "CON2", "CON1", "CON2", "CON1", "CON2", "CON1", "CON2", "CON1", "CON2", "CON1", "CON2", "CON1", "CON2", "CON1", "CON2", "CON1", "CON2", "CON1", "CON2", "CON1", "CON2")


# RNA-SEQ HITS.
dat <- read.table("/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/AA/AA_RNA_Seq_20140630/20140630-NewestRound/Analyses/Combined/AA37TxExp_RNA_Hits_Reordered_QN.txt", header = TRUE, sep = "\t")
dat.m <- t(as.matrix(dat[,2:length(dat[1,])]))
colnames(dat.m) <- dat[,1]
row.names(dat.m) <- colnames (dat[,2:length(dat[1,])])

par(mfrow = c(1,1))
colors = c("gold", "turquoise", "red", "black")
names(colors) = unique(Cat)
ecb = function(x,y){plot(x,t='n'); text(x,labels=row.names(dat.m), col=colors[Cat])}
tsne_dat = tsne(dat.m, epoch_callback = ecb, perplexity=50)
write.table(tsne_dat, file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/RESULTS/AA/AA_RNA_Seq_20140630/20140630-NewestRound/Analyses/Combined/AA37TxExp_RNA_Hits_Reordered_QN_tsne.txt", sep = "\t")
# k = 3.
tsne_dat = tsne(dat.m, k = 3, epoch_callback = ecb, perplexity=50)
s3d <- scatterplot3d(tsne_dat[,1], tsne_dat[,2], tsne_dat[,3], pch = 16, color = colors[Cat], xlab = "t-SNE_First_Dimension", ylab = "t-SNE_Second_Dimension", zlab = "t-SNE_Third_Dimension", main = "Asthma Twin RNA-Seq t-SNE")
s3d.coords <- s3d$xyz.convert(tsne_dat[,1], tsne_dat[,2], tsne_dat[,3])
text(s3d.coords$x, s3d.coords$y, labels=row.names(dat.m), pos=4, cex=.5, col = colors[Cat])
# Add drop lines.
s3d <- scatterplot3d(tsne_dat[,1], tsne_dat[,2], tsne_dat[,3], pch = 16, color = colors[Cat], type = "h", xlab = "t-SNE_First_Dimension", ylab = "t-SNE_Second_Dimension", zlab = "t-SNE_Third_Dimension", main = "Asthma Twin RNA-Seq t-SNE")
s3d.coords <- s3d$xyz.convert(tsne_dat[,1], tsne_dat[,2], tsne_dat[,3])
text(s3d.coords$x, s3d.coords$y, labels=row.names(dat.m), pos=4, cex=.5, col = colors[Cat])
write.table(tsne_dat, file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/RESULTS/AA/AA_RNA_Seq_20140630/20140630-NewestRound/Analyses/Combined/AA37TxExp_RNA_Hits_Reordered_QN_k3_tsne.txt", sep = "\t")


# RNA-SEQ HITS IN DISCORDANT TWINS ONLY.
dat <- read.table("/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/AA/AA_RNA_Seq_20140630/20140630-NewestRound/Analyses/Combined/AA_DA_RNA_Hits_Reordered_QN.txt", header = TRUE, sep = "\t")
dat.m <- t(as.matrix(dat[,2:length(dat[1,])]))
colnames(dat.m) <- dat[,1]
row.names(dat.m) <- colnames (dat[,2:length(dat[1,])])

par(mfrow = c(1,1))
colors = c("gold", "turquoise")
names(colors) = unique(Cat[1:24])
ecb = function(x,y){plot(x,t='n'); text(x,labels=row.names(dat.m), col=colors[Cat[1:24]])}
tsne_dat = tsne(dat.m, epoch_callback = ecb, perplexity=50)
write.table(tsne_dat, file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/RESULTS/AA/AA_RNA_Seq_20140630/20140630-NewestRound/Analyses/Combined/AA37TxExp_RNA_Hits_DA_Reordered_QN_tsne.txt", sep = "\t")
# k = 3.
tsne_dat = tsne(dat.m, k = 3, epoch_callback = ecb, perplexity=50)
s3d <- scatterplot3d(tsne_dat[,1], tsne_dat[,2], tsne_dat[,3], pch = 16, color = colors[Cat[1:24]], xlab = "t-SNE_First_Dimension", ylab = "t-SNE_Second_Dimension", zlab = "t-SNE_Third_Dimension", main = "Asthma Twin RNA-Seq t-SNE")
s3d.coords <- s3d$xyz.convert(tsne_dat[,1], tsne_dat[,2], tsne_dat[,3])
text(s3d.coords$x, s3d.coords$y, labels=row.names(dat.m), pos=4, cex=.5, col = colors[Cat[1:24]])
# Add drop lines.
s3d <- scatterplot3d(tsne_dat[,1], tsne_dat[,2], tsne_dat[,3], pch = 16, color = colors[Cat[1:24]], type = "h", xlab = "t-SNE_First_Dimension", ylab = "t-SNE_Second_Dimension", zlab = "t-SNE_Third_Dimension", main = "Asthma Twin RNA-Seq t-SNE")
s3d.coords <- s3d$xyz.convert(tsne_dat[,1], tsne_dat[,2], tsne_dat[,3])
text(s3d.coords$x, s3d.coords$y, labels=row.names(dat.m), pos=4, cex=.5, col = colors[Cat[1:24]])
write.table(tsne_dat, file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/RESULTS/AA/AA_RNA_Seq_20140630/20140630-NewestRound/Analyses/Combined/AA37TxExp_RNA_Hits_DA_Reordered_QN_k3_tsne.txt", sep = "\t")


# CIBERSORT RESULTS.
datcs <- read.table("/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/RESULTS/AA/AA_RNA_Seq_20140630/20140630-NewestRound/Analyses/Combined/00-CIBERSORT_Results/CIBERSORT.Output_Job2.txt", header = TRUE, sep = "\t")
datcs.m <- as.matrix(datcs[,2:length(datcs[1,])])
row.names(datcs.m) <- datcs[,1]
col.name <- colnames (datcs.m)

par(mfrow = c(1,1))
colors = c("gold", "turquoise", "red", "black")
names(colors) = unique(Cat)
ecb = function(x,y){plot(x,t='n'); text(x,labels=row.names(datcs.m), col=colors[Cat])}
tsne_datcs = tsne(datcs.m, epoch_callback = ecb, perplexity=50)
write.table(tsne_datcs, file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/RESULTS/AA/AA_RNA_Seq_20140630/20140630-NewestRound/Analyses/Combined/AA37TxExp_CIBERSORT_tsne.txt", sep = "\t")
# k = 3.
tsne_dat = tsne(datcs.m, k = 3, epoch_callback = ecb, perplexity=50)
s3d <- scatterplot3d(tsne_dat[,1], tsne_dat[,2], tsne_dat[,3], pch = 16, color = colors[Cat], xlab = "t-SNE_First_Dimension", ylab = "t-SNE_Second_Dimension", zlab = "t-SNE_Third_Dimension", main = "Asthma Twin CIBERSORT t-SNE")
s3d.coords <- s3d$xyz.convert(tsne_dat[,1], tsne_dat[,2], tsne_dat[,3])
text(s3d.coords$x, s3d.coords$y, labels=row.names(datcs.m), pos=4, cex=.5, col = colors[Cat])
# Add drop lines.
s3d <- scatterplot3d(tsne_dat[,1], tsne_dat[,2], tsne_dat[,3], pch = 16, color = colors[Cat], type = "h", xlab = "t-SNE_First_Dimension", ylab = "t-SNE_Second_Dimension", zlab = "t-SNE_Third_Dimension", main = "Asthma Twin CIBERSORT t-SNE")
s3d.coords <- s3d$xyz.convert(tsne_dat[,1], tsne_dat[,2], tsne_dat[,3])
text(s3d.coords$x, s3d.coords$y, labels=row.names(datcs.m), pos=4, cex=.5, col = colors[Cat])
write.table(tsne_dat, file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/RESULTS/AA/AA_RNA_Seq_20140630/20140630-NewestRound/Analyses/Combined/AA37TxExp_CIBERSORT_k3_tsne.txt", sep = "\t")


# RNA-SEQ FULL RESULTS.
# TOO MANY DIMENSIONS, CANNOT RUN EVEN ON THE MACPRO.
dat <- read.table("/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/AA/AA_RNA_Seq_20140630/20140630-NewestRound/Analyses/Combined/AA37TxExp_Reordered_QN.txt", header = TRUE, sep = "\t")
dat.m <- t(as.matrix(dat[,2:length(dat[1,])]))
colnames(dat.m) <- dat[,1]
row.names(dat.m) <- colnames (dat[,2:length(dat[1,])])

par(mfrow = c(1,1))
colors = c("gold", "turquoise", "red", "black")
names(colors) = unique(Cat)
ecb = function(x,y){plot(x,t='n'); text(x,labels=row.names(dat.m), col=colors[Cat])}
tsne_dat = tsne(dat.m, epoch_callback = ecb, perplexity=50)
write.table(tsne_dat, file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/RESULTS/AA/AA_RNA_Seq_20140630/20140630-NewestRound/Analyses/Combined/AA37TxExp_Reordered_QN_tsne.txt", sep = "\t")




