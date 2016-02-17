setwd('/home/fortelny/OtherAnalysis/2015_07_03_mRNAvsProtein/DraftMapOfHumanProteome/final')

library(gdata)
library(ggplot2)
library(reshape)
library(reshape2)

load("3_Mats.RData")

n <- 2000

# HOW STRONG IS CORRELATION OF THIS mRNA VS ALL OTHERS?  -------------------
# get all correlations (parallelized)
library(doMC)
registerDoMC(cores = 11)
idces <- 1:length(proteins)
corVecs <- foreach(p.idx=idces) %dopar% {
  # correlate protein with all mRNAs
  return(sapply(c(p.idx, sample(idces[-p.idx], n)), function(mRNA.idx) 
    cor(mats$prots[p.idx,], mats$mRNAs[mRNA.idx,], method="spearman", use="pairwise.complete.obs")
  ))
}
save(corVecs, file="12_Cors.RData")

# analyze
load("12_Cors.RData")
str(corVecs)
corVecs2 <- corVecs[which(sapply(corVecs, function(v) sum(abs(round(v,2)) != 1, na.rm=T)) > 0)] 
pDat <- data.frame(Fraction=sapply(corVecs2, function(v) sum(v[-1]>v[1], na.rm=T)/length(v[-1])))
ggplot(pDat, aes(x=Fraction)) + stat_ecdf() + 
  theme_bw(24) + 
  xlab(paste0("Fraction of random mRNAs (n=", n, " draws) \n with higher correlation")) + 
  xlim(0,1) + ylab("ECDF") + 
  ggtitle(paste0(sum(pDat$Fraction < 0.05), " of ", length(corVecs2), " are < 0.05"))
ggsave("12_Pvals_ECDF.pdf")
