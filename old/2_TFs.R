setwd("~/OtherAnalysis/2015_07_03_mRNAvsProtein/DraftMapOfHumanProteome/")

# this file uses transcription factors to predict protein. 
# Currently only the first 4 TFs are used because using all would have to many features
# at this point this is not very useful yet

load("1_mats.RData")
options(stringsAsFactors=F)

# WHICH GENE TO LOOK AT FOR EACH PROTEIN? ---------------------------------
str(tfs <- read.table("data/TranscriptionFactors.txt", sep = "\t", header =T)$To)
str(tfs.idx <- which(proteins %in% tfs)[1:4])

# ChipDB
(load("data/chipbase.RData"))
str(chipBase)
sum(chipBase$tf.ensg %in% proteins)
str(chipBase2 <- subset(chipBase, target.ensg %in% proteins & tf.ensg %in% proteins))

# LINEAR MODEL ------------------------------------------------------------

mats$lm_TFs <- matrix(NA, nrow=length(proteins), ncol=length(tissues))
for(p.idx in 1:length(proteins)){
  predDat <- data.frame(prot=mats$prots[p.idx, ], t(mats$mRNAs[unique(c(p.idx, tfs.idx)),]))
  for(t.idx in 1:length(tissues)){
    lmFit <- lm(prot ~ 1 + ., data=predDat[-t.idx,])
    mats$lm_TFs[p.idx, t.idx] <- predict(lmFit, predDat[t.idx,])
  }
}

# ELASTIC NET CARET-------------------------------------------------------------------
library(caret)
library(doMC)
registerDoMC(cores = 11)
mats$elastic_TFs <- matrix(NA, nrow=length(proteins), ncol=length(tissues))
for(p.idx in 1:length(proteins)){
  predDat <- data.frame(prot=mats$prots[p.idx, ], t(mats$mRNAs[unique(c(p.idx, tfs.idx)),]))
  for(t.idx in 1:length(tissues)){
    (elasticFit <- train(prot ~ ., data = predDat[-t.idx,],
                        method = "glmnet",
                        trControl = trainControl(method="LOOCV")))
#     as.matrix(coef(elasticFit$finalModel, s=elasticFit$bestTune$lambda))
    mats$elastic_TFs[p.idx, t.idx] <- predict(elasticFit, predDat[t.idx,])
  }
}

# ELASTIC NET WITHOUT INTERCEPT-------------------------------------------------------------------

mats$elastic_TFs_NI <- matrix(NA, nrow=length(proteins), ncol=length(tissues))
for(p.idx in 1:length(proteins)){
  predDat <- data.frame(prot=mats$prots[p.idx, ], t(mats$mRNAs[unique(c(p.idx, tfs.idx)),]))
  for(t.idx in 1:length(tissues)){
    (elasticFit <- train(prot ~ ., data = predDat[-t.idx,],
                         method = "glmnet",
                         trControl = trainControl(method="LOOCV"),
                         intercept=F
                         ))
#   as.matrix(coef(elasticFit$finalModel, s=elasticFit$bestTune$lambda))
    mats$elastic_TFs_NI[p.idx, t.idx] <- predict(elasticFit, predDat[t.idx,])
  }
}

# SAVE --------------------------------------------------------------------

save(mats, tissues, proteins, file="2_mats.RData")
