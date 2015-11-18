############################################################
# TO-DO LIST.
############################################################
DONE. # 1. Determine REF cutoff by distribution.





############################################################
# A. Library Used.
############################################################
library("gplots")
library("ggplot2")





############################################################
# B. DETERMINE CUTOFF BY REF DISTRIBUTION.
############################################################
############################################################
# B.1. Raw Data Loading.
############################################################
# Load data.
dat <- read.table("/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/AA/AA_Proteomics_20150729/Results/20151112-REF_Cutoff/AA_PRT_Combined_Results_ANNO_DA_20151110.txt",  header =TRUE, sep = "\t")

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


############################################################
# B.2. REF DISTRIBUTION
############################################################
refall <- c(dat_data[,19], dat_data[,20], dat_data[,21], dat_data[,22], dat_data[,23], dat_data[,24])

# Get rid of all the zero numbers as they are all missing values.
refall <- refall[refall != 0]
summary(refall)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  0.124   0.935   1.000   1.023   1.063  47.330 

quantile(refall, probs = seq(0,1,0.01))
      0%       1%       2%       3%       4%       5%       6%       7%       8%       9%      10%      11%      12% 
 0.12400  0.61200  0.69000  0.73200  0.76400  0.78820  0.80600  0.82300  0.83600  0.84800  0.85900  0.86900  0.87400 
     13%      14%      15%      16%      17%      18%      19%      20%      21%      22%      23%      24%      25% 
 0.88200  0.88900  0.89400  0.90000  0.90500  0.90900  0.91400  0.91900  0.92100  0.92400  0.92800  0.93200  0.93500 
     26%      27%      28%      29%      30%      31%      32%      33%      34%      35%      36%      37%      38% 
 0.93800  0.94200  0.94500  0.94900  0.95200  0.95500  0.95800  0.96000  0.96300  0.96600  0.96800  0.97100  0.97300 
     39%      40%      41%      42%      43%      44%      45%      46%      47%      48%      49%      50%      51% 
 0.97600  0.97800  0.98000  0.98300  0.98500  0.98700  0.98900  0.99100  0.99300  0.99500  0.99800  1.00000  1.00200 
     52%      53%      54%      55%      56%      57%      58%      59%      60%      61%      62%      63%      64% 
 1.00400  1.00700  1.00900  1.01100  1.01300  1.01500  1.01800  1.02000  1.02200  1.02500  1.02700  1.03000  1.03200 
     65%      66%      67%      68%      69%      70%      71%      72%      73%      74%      75%      76%      77% 
 1.03500  1.03700  1.03900  1.04200  1.04500  1.04800  1.05100  1.05400  1.05700  1.06000  1.06300  1.06700  1.07100 
     78%      79%      80%      81%      82%      83%      84%      85%      86%      87%      88%      89%      90% 
 1.07500  1.07800  1.08220  1.08800  1.09300  1.09800  1.10400  1.11000  1.11700  1.12500  1.13400  1.14400  1.15600 
     91%      92%      93%      94%      95%      96%      97%      98%      99%     100% 
 1.16900  1.18500  1.20100  1.22300  1.25200  1.29000  1.35008  1.43800  1.67072 47.32800 

# Mean and SD of logged values.
mean_refall <- mean(log10(refall))
mean_refall
[1] -0.0002543774
stdv <- sd(log10(refall))
stdv
[1] 0.07859842

# Mean +/- 3 * SD.
mean_refall-3* stdv
[1] -0.2360496
mean_refall+3* stdv
[1] 0.2355409
# Mean +/- 2 * SD.
mean_refall-2* stdv
[1] -0.1574512
mean_refall+2* stdv
[1] 0.1569425

# 10 ^ (Mean +/- 3 * SD).
10^(mean_refall-3* stdv)
[1] 0.580698
10^(mean_refall+3* stdv)
[1] 1.720049
# 10 ^ (Mean +/- 2 * SD).
10^(mean_refall-2* stdv)
[1] 0.6959031
10^(mean_refall+2* stdv)
[1] 1.435299

# Histogram.
hist(log10(refall), nclass = 100, main = "Log10(AA Proteome Reference Ratios)", xlab = "Log10(Ref_Ratio)")
abline(v = mean_refall-3* stdv, col = "red", lty = 2)
abline(v = mean_refall+3* stdv, col = "red", lty = 2)
abline(v = mean_refall-2* stdv, col = "green", lty = 2)
abline(v = mean_refall+2* stdv, col = "green", lty = 2)




