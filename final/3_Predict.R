source("../0_SETUP.R")
setwd("final/")


# FUNCTIONS FOR PREDICTION ------------------------------------------------

addPredictionsToMats <- function(mats, proteins, tissues){
  
  # RATIO -------------------------------------------------------------------
  mats$Ratio <- matrix(NA, nrow=length(proteins), ncol=length(tissues))
  for(p.idx in 1:length(proteins)){
    ratio <- median(mats$prots[p.idx, ]/mats$mRNAs[p.idx, ], na.rm=T)
    for(t.idx in 1:length(tissues)){
      mats$Ratio[p.idx, t.idx] <- mats$mRNAs[p.idx, t.idx] * ratio
    }
  }
  
  # RATIO with Leave one out strategy-------------------------------------------------------------------
  mats$RatioLOO <- matrix(NA, nrow=length(proteins), ncol=length(tissues))
  for(p.idx in 1:length(proteins)){
    for(t.idx in 1:length(tissues)){
      ratio <- median(mats$prots[p.idx, -t.idx]/mats$mRNAs[p.idx, -t.idx], na.rm=T)
      mats$RatioLOO[p.idx, t.idx] <- mats$mRNAs[p.idx, t.idx] * ratio
    }
  }
  
  # CONTROL-----------------------
  mats$Control <- matrix(NA, nrow=length(proteins), ncol=length(tissues))
  for(t.idx in 1:length(tissues)){
    mats$Control[,t.idx] <- apply(mats$prots, 1, median, na.rm=T)
  }
  
  # CONTROL_LOO-----------------------
  mats$ControlLOO <- matrix(NA, nrow=length(proteins), ncol=length(tissues))
  for(t.idx in 1:length(tissues)){
    mats$ControlLOO[,t.idx] <- apply(mats$prots[, -t.idx], 1, median, na.rm=T)
  }
  
  return(mats)
}


# DO THE PREDICTION -------------------------------------------------------
load("1_Mats.RData")
length(mats)
mats <- addPredictionsToMats(mats, proteins, tissues)
length(mats)
save(mats, tissues, proteins, file="3_Mats.RData")
rm(mats)

load("2_Mats_norm.RData")
length(mats)
mats <- addPredictionsToMats(mats, proteins, tissues)
length(mats)
save(mats, tissues, proteins, file="3_Mats_norm.RData")
rm(mats)
