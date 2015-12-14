source("../0_SETUP.R")
setwd("final/")

library(reshape)
library(ggplot2)
library(gdata)
library(grid)

n = 5
sd=0.3
set.seed(19199)
pDat <- data.frame(
  Tissue=paste0("Tissue", 1:n),
  Protein = c(rnorm(n,1,sd=sd),rnorm(n,2.5,sd=sd),rnorm(n,4,sd=sd)),
  Gene = rep(c("A", "B", "C"), each = n),
  mRNA = c(rnorm(n, 3, sd=sd), rnorm(n, 4,sd=sd), rnorm(n,1,sd=sd))
  )

ratio <- list()
pDat$PredProtein <- NA
for(i in c("A", "B", "C")){
  d <- subset(pDat, Gene == i)
  ratio[[i]] <- median(d$Protein/d$mRNA)
  pDat$PredProtein[pDat$Gene == i] <- pDat$mRNA[pDat$Gene == i] * ratio[[i]]
}

p <- ggplot(pDat, aes(x=mRNA, y=Protein, color=Gene)) + 
  geom_point(size=3, shape=1) + 
  geom_point(data = subset(pDat, Tissue=="Tissue4"), shape=16, size=3)
for(i in names(ratio)){
  p <- p + geom_abline(slope=ratio[[i]], size=1, colour="grey")
}
p + geom_segment(
  data = subset(pDat, Tissue=="Tissue4"), 
  aes(x=mRNA, xend=mRNA, y=Protein, yend=PredProtein, colour=Gene)) + 
  geom_point(size=3, shape=1) + 
  geom_point(data = subset(pDat, Tissue=="Tissue4"), shape=16, size=3) + 
  ylab("Observed protein") + 
  xlim(0,6) + ylim(0,6) + theme_bw(24) + theme(legend.position="none")
ggsave("6_Example1.pdf", width=5, height=5)

ggplot(subset(pDat, Tissue=="Tissue4"), aes(x=PredProtein, y=Protein, color=Gene)) + 
  geom_point(size = 3) + 
  geom_abline(slope = 1, colour = "grey") + 
  geom_segment(aes(x=PredProtein, xend=PredProtein, y=Protein, yend=PredProtein, colour=Gene)) +
  xlim(0,6) + ylim(0,6) + theme_bw(24) + theme(legend.position="none") + 
  xlab("Predicted protein") + ylab("Observed protein")
ggsave("6_Example2.pdf", width=5, height=5)
