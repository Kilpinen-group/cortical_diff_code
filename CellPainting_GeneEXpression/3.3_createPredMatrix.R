setwd(".")

library(oobLegacy)

gexData<-fastRead("input/gexMetacellsPerCellLine.tsv",as.matrix = T)
gex2cpModel<-fastRead("results/GeXmodelMulti.tsv")

cpData<-fastRead("input/cpMetacellsPerCellLine.tsv",as.matrix = T)
cp2gexModel<-fastRead("results/CPmodelMulti.tsv")

computeTheoMat<-function(mat,model){
  nfeatureModel=nrow(model)-1
  res<-apply(mat,1,function(obs){
    apply(model,2,function(feature){
      feature[1]+sum(feature[2:length(feature)]*obs)
    })
  })
  t(res)
}


#create cp theo
scData<-readRDS("countTableCorrected.rds") |> t()
cpTheoMat<-computeTheoMat(scData, gex2cpModel)

saveRDS(cpTheoMat, "outRobj/cpTheoMat.rds")

rm(scData, cpTheoMat)

cpData<- fastRead("")  |> t()

scTheoMat<-computeTheoMat(scData, cp2gexModel)



