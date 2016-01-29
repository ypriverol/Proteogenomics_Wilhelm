library(gdata)
library(ggplot2)
library(reshape)
library(reshape2)

load("3_Mats.RData")

genes <- proteins

# within gene variance
mean(apply(mats$prots, 1, var,na.rm=T), na.rm=T)

# between gene variance
mean(apply(mats$prots, 2, var,na.rm=T), na.rm=T)


# NAs ---------------------------------------------------------------------

table(values.av <- sapply(genes, function(p){
  return(sum(!is.na(mats$prots[p,]) & !is.na(mats$mRNA[p,])))
}))
ggplot(data.frame(Values=factor(values.av)), aes(x=Values)) + geom_histogram() + theme_bw(24)
ggsave("8_Values.pdf", width=7, height=5)
names(values.av) <- genes
str(genes2 <- genes[which(values.av >= 3)])
str(aDat <- data.frame(Values = values.av[genes]))
rownames(aDat) <- genes

# SPEARMAN ----------------------------------------------------------------

aDat$Correlation <- sapply(genes, function(p){ return(cor(mats$prots[p,], mats$mRNAs[p,], 
                                                  use="pairwise.complete.obs", 
                                                  method="spearman"))})
# Boxplots
ggplot(aDat, aes(x=factor(Values), y=Correlation)) + geom_boxplot() + 
  theme_bw(24) + xlab("Values available") + ylab("Correlation (R) per gene")
ggsave("8_Spearman_boxplots.pdf", width=7, height=5)

ggplot(subset(aDat, Values >= 9), aes(x=Correlation)) + geom_histogram() + 
  theme_bw(24) + xlab("Correlation (R) per gene") + ylab("Count") + xlim(-1, 1)
ggsave("8_Spearman_histogram.pdf", width=5, height=5)
median(subset(aDat, Values >= 9)$Correlation, na.rm=T)


# Manual R2 -------------------------------------------------------------------

str(aDat$Error <- apply((mats$prots - mats$Ratio)**2, 1, sum, na.rm=T))
str(aDat$Variance <- apply((mats$prots - apply(mats$prots, 1, mean, na.rm=T))**2, 1, sum, na.rm=T))
# below is just to check that the whole "apply"-mess worked correctly
for(i in 1:length(proteins)){
  stopifnot(aDat$Variance[i] == sum((mats$prots[i,] - mean(mats$prots[i,], na.rm=T))**2, na.rm=T))  
}
x <- aDat$Error / aDat$Variance
quantile(x[which(aDat$Values >2)], na.rm=T)
# x[x > 1] = 1
aDat$R2 <- 1 - x

# Looking at one example
# subset(aDat, Values == 12 & R2 < 0)
# plot(mats$mRNAs["ENSG00000197111",], mats$prots["ENSG00000197111", ], xlim=c(0,0.0006), ylim=c(0,0.0006))
# aDat["ENSG00000197111",]
# summary(lm(mats$prots["ENSG00000197111", ] ~ mats$mRNAs["ENSG00000197111", ] - 1))
# y <- mats$prots["ENSG00000197111", ]
# x <- mats$mRNAs["ENSG00000197111", ]
# str(yhat <- median(y/x, na.rm=F) * x)
# r2 <- 1-(sum((y-yhat)**2)/sum((y-mean(y))**2))
# r2
# abline(a=0, b=median(y/x, na.rm=F))
# abline(a=0, b=lm(y~x-1)$coef, col="red")
# abline(lm(y~x), col="blue")
# 
# ls.1<-lm(y~x-1)
# summary(ls.1)
# yhat.ls<-ls.1$fit
# 
# ls.2<-lm(y~x)
# 
# r2.ls <- 1-(sum((y-yhat.ls)**2)/sum((y-0)**2))
# r2.ls

# Boxplots
ggplot(aDat, aes(x=factor(Values), y=R2)) + geom_boxplot() + 
  theme_bw(24) + xlab("R2") + ylab("Count")

ggplot(aDat, aes(x=factor(Values), y=R2)) + geom_boxplot() + 
  theme_bw(24) + xlab("R2") + ylab("Count")


ggplot(aDat, aes(x=factor(Values), y=R2)) + geom_boxplot() + 
  theme_bw(24) + xlab("Values available") + ylab("R squared") + ylim(-2, 1)
ggsave("8_R2_boxplots.pdf", width=7, height=5)

sapply(split(aDat$R2, factor(aDat$Values)), min, na.rm=T)

# Histogram
ggplot(subset(aDat, Values >= 9), aes(x=R2)) + geom_histogram(binwidth=0.05) + 
  theme_bw(24) + xlab("R squared") + ylab("Count") + xlim(-2, 1)
ggsave("8_R2_histogram.pdf", width=5, height=5)
median(subset(aDat, Values >= 9)$R2, na.rm=T)
sum(subset(aDat, Values >= 9)$R2 < 0, na.rm=T)/sum(!is.na(subset(aDat, Values >= 9)$R2))

# A proper LS model (slope only) FOR EACH GENE -------------------------------------------------------

pModels <- lapply(genes2, function(p){ return(lm(mats$prots[p,] ~ mats$mRNAs[p,] - 1)) })
names(pModels) <- genes2
# get R2
aDat$R2_LS <- sapply(pModels[which(!is.na(pModels))], function(m){return(summary(m)$r.squared)})[genes]

# Boxplots
ggplot(aDat, aes(x=factor(Values), y=R2_LS)) + geom_boxplot() + 
  theme_bw(24) + xlab("Available values") + ylab("R squared")
ggsave("8_R2_LS_boxplots.pdf", width=7, height=5)

# Histogram
ggplot(subset(aDat, Values >= 9), aes(x=R2_LS)) + geom_histogram(binwidth=0.01) + 
  theme_bw(24) + xlab("R squared") + ylab("Count") + xlim(0, 1)
ggsave("8_R2_LS_histogram.pdf", width=5, height=5)
median(subset(aDat, Values >= 9)$R2_LS, na.rm=T)



# A proper LS model (intercept and slope) FOR EACH GENE -------------------------------------------------------

pModels2 <- lapply(genes2, function(p){ return(lm(mats$prots[p,] ~ mats$mRNAs[p,])) })
names(pModels2) <- genes2
summary(pModels2[[1]])
# get R2
aDat$R2_LS2 <- sapply(pModels2[which(!is.na(pModels2))], function(m){return(summary(m)$r.squared)})[genes]

# Boxplots
ggplot(aDat, aes(x=factor(Values), y=R2_LS2)) + geom_boxplot() + 
  theme_bw(24) + xlab("Available values") + ylab("R squared")
ggsave("8_R2_LS2_boxplots.pdf", width=7, height=5)

# Histogram
ggplot(subset(aDat, Values >= 9), aes(x=R2_LS2)) + geom_histogram(binwidth=0.01) + 
  theme_bw(24) + xlab("R2") + ylab("Count") + xlim(0, 1)
ggsave("8_R2_LS2_histogram.pdf", width=5, height=5)
median(subset(aDat, Values >= 9)$R2_LS2, na.rm=T)

# # OVERFIT EFFECT ----------------------------------------------------------

xx <- data.frame(
  claimed_MSE = apply((mats$prots-mats$Ratio)**2, 1, function(row){return(mean(row, na.rm=T))}),
  overfit_cor_MSE = apply((mats$prots-mats$RatioLOO)**2, 1, function(row){return(mean(row, na.rm=T))}),
  cnt_prot = apply(mats$prots, 1, function(r) sum(!is.na(r))),
  cnt_gene = apply(mats$mRNAs, 1, function(r) sum(!is.na(r))),
  cnt = apply(!is.na(mats$mRNAs) & !is.na(mats$prots), 1 ,sum)
)
quantile(xx$overfit_effect <- (xx$overfit_cor_MSE/xx$claimed_MSE), na.rm=T)

ggplot(xx, aes(y=overfit_effect, group=cnt, x=factor(cnt))) + geom_boxplot() + 
  scale_y_log10() + 
  xlab("Available values") + ylab("Overfit") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme_bw(24)
ggsave(paste0("8_OverfitEffect",label,".pdf"), height=5, width=7)
quantile(subset(xx, cnt==2)$overfit_effect, na.rm=T)
head(subset(xx, cnt==2, select=c(claimed_MSE, overfit_cor_MSE, overfit_effect)))
sapply(split(xx$overfit_effect, factor(xx$cnt)), function(v) return(median(v, na.rm=T)))

# # CORRELATIONS BY PROTEIN ACROSS SAMPLES ------------------------------------
# pCorMats <- list()
# for(matNam in c("mRNAs", "Ratio")){
#   pCorMats[[matNam]] <- sapply(1:length(proteins), function(p.idx){cor(mats$prots[p.idx,], mats[[matNam]][p.idx,], method="spearman", use="pairwise.complete.obs")})
# }
# sapply(pCorMats, function(vec) sum(!is.na(vec)))
# str(pDat <- melt(pCorMats))
# ggplot(pDat, aes(x=L1, y=value)) + geom_boxplot(fill = "lightgrey") + 
#   xlab("Method") + ylab("Correlation (R)") + 
#   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#   theme_bw(24)
# ggsave(paste0("4_AcrossSamples",label,".pdf"), width=5, height=5)