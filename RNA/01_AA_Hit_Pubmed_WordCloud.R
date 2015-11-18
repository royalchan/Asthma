library("PubMedWordcloud")

# Gene Names in PubMed.
dat <- read.table("/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/AA/AA_Methyl_Seq/Analyses/20141208-From_George/20150528-AA_RNAMetHit_LitScraping/AA_Hit_Genes_Symbol_Cleaned.txt")
dat <- as.character(dat[,1])

# Empirically needs to get rid of the following entries by hand as they break getPMIDsByKeyWords by having no hits:
# C1ORF43	C8ORF82	UBE2NL	C4ORF33	LOC100130298	ANXA2P3	METTL12	ZNF561	ZNF765	ZFP41	ZNF416	HEATR5B	YLPM1	FBXO46	ZNF564	ZNF586	ZNF609	ZNF766	KIAA1191	OARD1	ZNF134	CCDC169	OR10P1	SECISBP2L.
for (i in 1:length(dat)) {
	aa_pmid <- c()
	aa_abstract <- c()
	aa_abstract_cleaned <- c()
	aa_pmid <- getPMIDsByKeyWords(keys = dat[i])
	if(length(aa_pmid) > 900) {k = 900} else {k = length(aa_pmid)}
	aa_abstract <- getAbstracts(aa_pmid[1:k])
	aa_abstract_cleaned <- cleanAbstracts(aa_abstract)
	filename <- paste("/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/AA/AA_Methyl_Seq/Analyses/20141208-From_George/20150528-AA_RNAMetHit_LitScraping/word_clouds/", dat[i], ".png", sep = "")
	png(filename, width = 960, height = 960, units = "px")
		plotWordCloud(aa_abstract_cleaned, scale = c(10, 0.5))
		text(0.5, 1, dat[i], cex = 5)
	dev.off()
}


# Gene Names and Asthma in PubMed.
dat2 <- read.table("/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/AA/AA_Methyl_Seq/Analyses/20141208-From_George/20150528-AA_RNAMetHit_LitScraping/AA_Hit_Genes_Symbol_Asthma.txt")
dat2 <- as.character(dat2[,1])

for (i in 1:length(dat2)) {
	aa_pmid <- c()
	aa_abstract <- c()
	aa_abstract_cleaned <- c()
	aa_pmid <- getPMIDsByKeyWords(keys = paste(dat2[i], "asthma", sep = " "))
	if(length(aa_pmid) > 900) {k = 900} else {k = length(aa_pmid)}
	aa_abstract <- getAbstracts(aa_pmid[1:k])
	aa_abstract_cleaned <- cleanAbstracts(aa_abstract)
	filename <- paste("/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/AA/AA_Methyl_Seq/Analyses/20141208-From_George/20150528-AA_RNAMetHit_LitScraping/word_clouds/", dat2[i], "_asthma.png", sep = "")
	png(filename, width = 960, height = 960, units = "px")
		plotWordCloud(aa_abstract_cleaned, scale = c(5, 0.5))
		text(0.5, 1, dat[i], cex = 5)
	dev.off()
}




