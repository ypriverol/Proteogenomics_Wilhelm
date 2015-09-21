setwd("~/OtherAnalysis/2015_07_03_mRNAvsProtein/DraftMapOfHumanProteome/")

files <- list.files(".", recursive=T)
rFiles <- files[which(regexpr("^3\\_.+?\\.R$", files) >0)]
(rFiles <- rFiles[which(regexpr("0_Pipeline\\.R$", rFiles) < 1)])

for(file in rFiles){
  print(file)
  knitr::stitch_rhtml(file, output=sub("\\.R$", "\\.html", file))  
}