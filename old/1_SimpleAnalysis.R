setwd("~/OtherAnalysis/2015_07_03_mRNAvsProtein/DraftMapOfHumanProteome/")

# this file predicts protein from mRNA with a variety of simple methods, simple meaning that these methods only use the mRNA for that protein as predictor

load("0_0_data.RData")

tissues <- colnames(geneUN)
table(apply(!is.na(geneUN), 1, sum), apply(!is.na(protUN), 1, sum))
str(proteins <- rownames(geneUN)[which(apply(!is.na(geneUN), 1, sum) == 12 & apply(!is.na(protUN), 1, sum) == 12)])

set.seed(8945)

# MATRICES ----------------------------------------------------------------

mats <- list(
  mRNAs = geneUN[proteins, tissues],
  prots = protUN[proteins, tissues]
  )


# RATIO -------------------------------------------------------------------

mats$ratio <- matrix(NA, nrow=length(proteins), ncol=length(tissues))
for(p.idx in 1:length(proteins)){
  predDat <- data.frame(prot=mats$prots[p.idx, ], mRNA=mats$mRNAs[p.idx, ])
  ratio <- median(mats$prots[p.idx, ]/mats$mRNAs[p.idx, ])
  for(t.idx in 1:length(tissues)){
    mats$ratio[p.idx, t.idx] <- mats$mRNAs[p.idx, t.idx] * ratio
  }
}

# RATIO with Leave one out strategy-------------------------------------------------------------------

mats$ratio_LOO <- matrix(NA, nrow=length(proteins), ncol=length(tissues))
for(p.idx in 1:length(proteins)){
  predDat <- data.frame(prot=mats$prots[p.idx, ], mRNA=mats$mRNAs[p.idx, ])
  for(t.idx in 1:length(tissues)){
    ratio <- median(mats$prots[p.idx, -t.idx]/mats$mRNAs[p.idx, -t.idx])
    mats$ratio_LOO[p.idx, t.idx] <- mats$mRNAs[p.idx, t.idx] * ratio
  }
}


# LINEAR MODEL WITHOUT INTERCEPT ------------------------------------------------------------

mats$lm_NI <- matrix(NA, nrow=length(proteins), ncol=length(tissues))
for(p.idx in 1:length(proteins)){
  predDat <- data.frame(prot=mats$prots[p.idx, ], mRNA=mats$mRNAs[p.idx, ])
  for(t.idx in 1:length(tissues)){
    lmFit <- lm(prot ~ 0 + mRNA, data = predDat[-t.idx,])
    mats$lm_NI[p.idx, t.idx] <- predict(lmFit, predDat[t.idx,])
  }
}

# LINEAR MODEL WTH INTERCEPT ------------------------------------------------------------

mats$lm <- matrix(NA, nrow=length(proteins), ncol=length(tissues))
for(p.idx in 1:length(proteins)){
  predDat <- data.frame(prot=mats$prots[p.idx, ], mRNA=mats$mRNAs[p.idx, ])
  for(t.idx in 1:length(tissues)){
    lmFit <- lm(prot ~ 1 + mRNA, data = predDat[-t.idx,])
    mats$lm[p.idx, t.idx] <- predict(lmFit, predDat[t.idx,])
  }
}

# Only intercept -------------------------------------------------------------------

mats$interceptOnly <- matrix(NA, nrow=length(proteins), ncol=length(tissues))
for(p.idx in 1:length(proteins)){
  prot <- mats$prots[p.idx, ]
  for(t.idx in 1:length(tissues)){
    mats$interceptOnly[p.idx, t.idx] <- mean(prot[-t.idx])
  }
}


save(mats, tissues, proteins, file="1_mats.RData")
