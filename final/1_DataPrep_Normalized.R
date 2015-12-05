setwd("~/OtherAnalysis/2015_07_03_mRNAvsProtein/DraftMapOfHumanProteome/")

# predicting protein simply (only from the mRNA in question)

load("0_0_data.RData")

tissues <- colnames(geneUN)
table(apply(!is.na(geneUN), 1, sum), apply(!is.na(protUN), 1, sum))
str(proteins <- rownames(geneUN)) # using all genes
# proteins2 is for lm but not used currently
str(proteins2 <- rownames(geneUN)[which(apply(!is.na(geneUN), 1, sum) != 0 & apply(!is.na(protUN), 1, sum) != 0)])

set.seed(8945)

# MATRICES ----------------------------------------------------------------

mats <- list(
  mRNAs = geneUN[proteins, tissues],
  prots = protUN[proteins, tissues]
  )

library(limma)
boxplot(t(mats$mRNA)[,1:6])
apply(t(mats$mRNA)[,1:6], 2, mean)
apply(t(mats$mRNA)[,1:6], 2, sd)
boxplot(scale(t(mats$mRNA)[,1:6]))
apply(scale(t(mats$mRNA)[,1:6]), 2, mean)
apply(scale(t(mats$mRNA)[,1:6]), 2, sd)
boxplot(limma::normalizeQuantiles(t(mats$mRNA)[,1:6]))

mats <- lapply(mats, function(m){
  return(t(scale(t(m))))
})

apply(scale(t(mats$mRNA)[,1:6]), 2, mean)
apply(scale(t(mats$mRNA)[,1:6]), 2, sd)

apply(scale(t(mats$prots)[,1:6]), 2, mean, na.rm=T)
apply(scale(t(mats$prots)[,1:6]), 2, sd, na.rm=T)

# RATIO -------------------------------------------------------------------

mats$ratio <- matrix(NA, nrow=length(proteins), ncol=length(tissues))
for(p.idx in 1:length(proteins)){
  predDat <- data.frame(prot=mats$prots[p.idx, ], mRNA=mats$mRNAs[p.idx, ])
  ratio <- median(mats$prots[p.idx, ]/mats$mRNAs[p.idx, ], na.rm=T)
  for(t.idx in 1:length(tissues)){
    mats$ratio[p.idx, t.idx] <- mats$mRNAs[p.idx, t.idx] * ratio
  }
}

# RATIO with Leave one out strategy-------------------------------------------------------------------

mats$ratio_LOO <- matrix(NA, nrow=length(proteins), ncol=length(tissues))
for(p.idx in 1:length(proteins)){
  predDat <- data.frame(prot=mats$prots[p.idx, ], mRNA=mats$mRNAs[p.idx, ])
  for(t.idx in 1:length(tissues)){
    ratio <- median(mats$prots[p.idx, -t.idx]/mats$mRNAs[p.idx, -t.idx], na.rm=T)
    mats$ratio_LOO[p.idx, t.idx] <- mats$mRNAs[p.idx, t.idx] * ratio
  }
}

# 
# # LINEAR MODEL WITHOUT INTERCEPT ------------------------------------------------------------
# 
# mats$lm_NI <- matrix(NA, nrow=length(proteins), ncol=length(tissues))
# for(p.idx in 1:length(proteins)){
#   predDat <- data.frame(prot=mats$prots[p.idx, ], mRNA=mats$mRNAs[p.idx, ])
#   for(t.idx in 1:length(tissues)){
#     lmFit <- lm(prot ~ 0 + mRNA, data = predDat[-t.idx,])
#     mats$lm_NI[p.idx, t.idx] <- predict(lmFit, predDat[t.idx,])
#   }
# }
# 
# # LINEAR MODEL WTH INTERCEPT ------------------------------------------------------------
# 
# mats$lm <- matrix(NA, nrow=length(proteins), ncol=length(tissues))
# for(p.idx in 1:length(proteins)){
#   predDat <- data.frame(prot=mats$prots[p.idx, ], mRNA=mats$mRNAs[p.idx, ])
#   for(t.idx in 1:length(tissues)){
#     lmFit <- lm(prot ~ 1 + mRNA, data = predDat[-t.idx,])
#     mats$lm[p.idx, t.idx] <- predict(lmFit, predDat[t.idx,])
#   }
# }
# 
# # Only intercept -------------------------------------------------------------------
# 
# mats$interceptOnly <- matrix(NA, nrow=length(proteins), ncol=length(tissues))
# for(p.idx in 1:length(proteins)){
#   prot <- mats$prots[p.idx, ]
#   for(t.idx in 1:length(tissues)){
#     mats$interceptOnly[p.idx, t.idx] <- mean(prot[-t.idx], na.rm=T)
#   }
# }

save(mats, tissues, proteins, file="5_mats.RData")
