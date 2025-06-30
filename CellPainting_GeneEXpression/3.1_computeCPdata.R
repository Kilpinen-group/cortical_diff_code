setwd(".")

library(oob)

dir.create("out3.1", showWarnings = FALSE)
inputDir <- "zenodoArchive/"

cellAnnot<-fastRead(paste0(inputDir,"/","CP_cellAnnotFiltered.tsv"))
featureMat<-fastRead(paste0(inputDir,"/","CP_featureMatFiltered.tsv"),as.matrix = TRUE)

scLine<-c('d20L11.4','d20L47.2','d20L61.1','d20L61.2','d20L62.4','d20L82.6','d40L11.4','d40L61.2','d40L62.4','d40L82.6','d70L11.4','d70L61.2','d70L62.4','d70L82.6')
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

fastWrite(cpMetacellsPerCellLine,"out3.1/CPGEX_cpMetacellsPerCellLine.tsv")



