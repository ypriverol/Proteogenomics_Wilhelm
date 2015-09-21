setwd("~/OtherAnalysis/2015_07_03_mRNAvsProtein/DraftMapOfHumanProteome/")

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

# # ELASTIC NET MODEL - NO INTERCEPT ------------------------------------------------------------
# library(glmnet)
# mats$elastic_NI <- matrix(NA, nrow=length(proteins), ncol=length(tissues))
# for(p.idx in 1:length(proteins)){
#   predDat <- data.frame(mRNA=mats$mRNAs[p.idx, ])
#   for(t.idx in 1:length(tissues)){
#     elasticFit <- glmnet(
#       x = data.matrix(predDat[-t.idx,]), 
#       y= mats$prots[p.idx, ], 
#       alpha = 0.5, 
#       intercept = F)
#     mats$elastic_NI[p.idx, t.idx] <- predict(elasticFit, newx=data.matrix(predDat[t.idx,]), s=c(0.01))
#   }
# }
# 
# # LASSO MODEL - NO INTERCEPT ------------------------------------------------------------
# library(glmnet)
# mats$lasso_NI <- matrix(NA, nrow=length(proteins), ncol=length(tissues))
# for(p.idx in 1:length(proteins)){
#   predDat <- data.frame(mRNA=mats$mRNAs[p.idx, ])
#   for(t.idx in 1:length(tissues)){
#     x <- predDat[-t.idx,]
#     lassoFit <- glmnet(
#       x = data.matrix(predDat[-t.idx,]), 
#       y= mats$prots[p.idx, ], 
#       alpha = 1,
#       intercept=F)
#     mats$lasso_NI[p.idx, t.idx] <- predict(lassoFit, newx=data.matrix(predDat[t.idx,]), s=c(0.01))
#   }
# }
# 
# # ELASTIC NET MODEL ------------------------------------------------------------
# library(glmnet)
# mats$elastic <- matrix(NA, nrow=length(proteins), ncol=length(tissues))
# for(p.idx in 1:length(proteins)){
#   predDat <- data.frame(prot=mats$prots[p.idx, ],  mRNA=mats$mRNAs[p.idx, ])
#   for(t.idx in 1:length(tissues)){
#     elasticFit <- glmnet(x = data.matrix(predDat[-t.idx,]), y= predDat[-t.idx,]$prot, alpha = 0.5)
#     mats$elastic[p.idx, t.idx] <- predict(
#       elasticFit, 
#       newx=data.matrix(predDat[t.idx,]), 
#       s=c(elasticFit$lambda[10]))
#   }
# }
# 
# # ELASTIC NET MODEL -- > CARET ------------------------------------------------------------
# library(caret)
# library(glmnet)
# mats$elastic <- matrix(NA, nrow=length(proteins), ncol=length(tissues))
# for(p.idx in 1:length(proteins)){
#   predDat <- data.frame(prot=mats$prots[p.idx, ],  mRNA=mats$mRNAs[p.idx, ])
#   for(t.idx in 1:length(tissues)){
#     elasticFit <- train(prot ~ mRNA, data = predDat,
#                         method = "glmnet",
#                         trControl = trainControl(
#                           method="repeatedcv",
#                           repeats=3,
#                           number=3))
#     elasticFit <- glmnet(x = data.matrix(predDat[-t.idx,]), y= predDat[-t.idx,]$prot, alpha = 0.5)
#     mats$elastic[p.idx, t.idx] <- predict(
#       elasticFit, 
#       newx=data.matrix(predDat[t.idx,]), 
#       s=c(elasticFit$lambda[10]))
#   }
# }
# 
# # LASSO MODEL ------------------------------------------------------------
# library(glmnet)
# mats$lasso <- matrix(NA, nrow=length(proteins), ncol=length(tissues))
# for(p.idx in 1:length(proteins)){
#   predDat <- data.frame(prot=mats$prots[p.idx, ],  mRNA=mats$mRNAs[p.idx, ])
#   for(t.idx in 1:length(tissues)){
#     x <- predDat[-t.idx,]
#     lassoFit <- glmnet(x = data.matrix(predDat[-t.idx,]), y= predDat[-t.idx,]$prot, alpha = 1)
#     mats$lasso[p.idx, t.idx] <- predict(lassoFit, newx=data.matrix(predDat[t.idx,]), s=c(0.01))
#   }
# }


save(mats, tissues, proteins, file="1_mats.RData")
