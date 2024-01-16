setwd(".")

library(oobLegacy)

cellAnnot<-fastRead("cellAnnotFiltered.tsv")
featureMat<-fastRead("featureMatFiltered.tsv",as.matrix = TRUE)

scLine<-c('d20_11.4','d20_47.2','d20_61.1','d20_61.2','d20_62.4','d20_82.6','d40_11.4','d40_61.2','d40_62.4','d40_82.6','d70_11.4','d70_61.2','d70_62.4','d70_82.6')
cpLine<-cellAnnot$expPop |> as.factor() |> levels()

commonLine<-intersect(scLine,cpLine)

selectedCells<-cellAnnot$expPop%in%commonLine

cellAnnot<-cellAnnot[selectedCells,]
featureMat<-featureMat[,selectedCells]


###doing metacells
#3 metacell per exp Pop
cellPerExpPop<-factorToVectorList(cellAnnot$expPop,factorNames = rn(cellAnnot))

metaCellPerExpPop<-lapply(names(cellPerExpPop),function(level){
  res<-paste0(level,".m", sample.int(3,length(cellPerExpPop[[level]]),replace = T))
  names(res)<-cellPerExpPop[[level]]
  res
}) ; names(metaCellPerExpPop)<-names(cellPerExpPop)

metaCells<-unlist(metaCellPerExpPop)
names(metaCells)<-unlist(lapply(metaCellPerExpPop,names))
cellAnnot$metaCells<-metaCells[rn(cellAnnot)]

cpMetacellsPerCellLine<-aggregMatPerVector(featureMat,cellAnnot$metaCells,FUN = mean) |> t()

fastWrite(cpMetacellsPerCellLine,"input/cpMetacellsPerCellLine.tsv")



