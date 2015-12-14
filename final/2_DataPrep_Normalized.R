# predicting protein simply (only from the mRNA in question)
source("../0_SETUP.R")
setwd("final/")

load("1_Mats.RData")

library(limma)
boxplot(t(mats$mRNA)[,1:6])
apply(t(mats$mRNA)[,1:6], 2, mean)
apply(t(mats$mRNA)[,1:6], 2, sd)
boxplot(scale(t(mats$mRNA)[,1:6]))
apply(scale(t(mats$mRNA)[,1:6]), 2, mean)
apply(scale(t(mats$mRNA)[,1:6]), 2, sd)

mats <- lapply(mats, function(m){
  return(t(scale(t(m))))
})

apply(scale(t(mats$mRNA)[,1:6]), 2, mean)
apply(scale(t(mats$mRNA)[,1:6]), 2, sd)

apply(scale(t(mats$prots)[,1:6]), 2, mean, na.rm=T)
apply(scale(t(mats$prots)[,1:6]), 2, sd, na.rm=T)

save(mats, tissues, proteins, file="2_Mats_norm.RData")
