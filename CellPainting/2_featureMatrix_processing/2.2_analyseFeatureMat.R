setwd(".")  #change to your working dir

library(oob)
library(GSDS)
library(stringr)

set.seed(666)
dir.create("out2.2",showWarnings = F)

cellAnnot<- fastRead("zenodoArchive/CP_cellAnnotFiltered.tsv")
cellAnnot<-oob::formatAnnotFromMeta(annotDataFrame = cellAnnot,metaAnnot = fastRead("zenodoArchive/CP_metaAnnot.tsv"))
colorScales<-attr(cellAnnot,"colorScales")
colorScales$cellClusters<-c(
    "d20.endoRetNeg"="#F3E600",
    "d20.nucNeg"="#783329",
    "d20.mito"="#D4686A",
    "d70_20.mito"="#E53C23",
    "d20.mature"="#F59B06",
    "d40.early"="#55BDBD",
    "d40.smallCells"="#404899",
    "d40.endoRet"="#007E80",
    "d40.mature"="#2C2968",
    "d70.bigcells"="#6EB52C",
    "d70.mitoNeg"="#115A2A"
)

featureMat<-fastRead("zenodoArchive/CP_featureMatFiltered.tsv", as.matrix = TRUE)

#umap<-UMAP(featureMat,ret_nn = T)

## create new cluster will overwrite the ones used in this script
#cellAnnot$leidenClust<-leidenFromUMAP(umap,n_neighbors = 20, resolution_parameter = 1)

# cellClusters<-rep("",nrow(cellAnnot));names(cellClusters)<-rn(cellAnnot)
# cellClusters[rn(cellAnnot)[cellAnnot$leidenClust=="k01"]]<-"d40.smallCells"
# cellClusters[rn(cellAnnot)[cellAnnot$leidenClust=="k02"]]<-"d20.mature"
# cellClusters[rn(cellAnnot)[cellAnnot$leidenClust=="k03"]]<-"d40.endoRet"
# cellClusters[rn(cellAnnot)[cellAnnot$leidenClust=="k04"]]<-"d20.endoRetNeg"
# cellClusters[rn(cellAnnot)[cellAnnot$leidenClust=="k05"]]<-"d40.mature"
# cellClusters[rn(cellAnnot)[cellAnnot$leidenClust=="k06"]]<-"d70.mitoNeg"
# cellClusters[rn(cellAnnot)[cellAnnot$leidenClust=="k07"]]<-"d20.mito"
# cellClusters[rn(cellAnnot)[cellAnnot$leidenClust=="k08"]]<-"d70_20.mito"
# cellClusters[rn(cellAnnot)[cellAnnot$leidenClust=="k09"]]<-"d40.early"
# cellClusters[rn(cellAnnot)[cellAnnot$leidenClust=="k10"]]<-"d70.bigcells"
# cellClusters[rn(cellAnnot)[cellAnnot$leidenClust=="k11"]]<-"d20.nucNeg"

# cellAnnot$cellClusters<-factor(cellClusters,levels = names(colorScales$cellClusters))

# colnames(umap$embedding)<-c("V1","V2")
# fastWrite(umap$embedding,"out2.2/umapCoord.tsv")
umap <- fastRead("zenodoArchive/CP_umapCoord.tsv",as.matrix = T,row.names = NULL)

pdf("out2.2/projAnnot.pdf",width=13, height=12)
proj2d(umap, useScatterMore=T,colorBy = cellAnnot$diffDay,colorScale = colorScales$diffDay)
proj2d(umap, useScatterMore=T,colorBy = cellAnnot$cellClusters,plotFactorsCentroids = T, colorScale = colorScales$cellClusters)
proj2d(umap, useScatterMore=T,colorBy = cellAnnot$parentalCellLine,colorScale = colorScales$parentalCellLine)
proj2d(umap, useScatterMore=T,colorBy = cellAnnot$expPop,colorScale = colorScales$expPop)
dev.off()

dir.create("out2.2/projPerFeature",showWarnings = F)
for(f in rn(featureMat)){
    pdf(paste0("out2.2/projPerFeature/",f,".pdf"),width=13, height=12)
    proj2d(umap, useScatterMore=T,colorBy = featureMat[f,],
           colorScale = c("darkblue","blue","cyan","lightgreen","orange","red","darkred"))
    dev.off()
}


## correlation matrix well to well
wellAnnot<-fastRead("zenodoArchive/CP_wellAnnot.tsv")
cellAnnot$wellKey<-paste0(cellAnnot$diffDay,substr(cellAnnot$well,1,1),formatNumber2Character(substr(cellAnnot$well,2,4)))

featureMatPerWell<-aggregMatPerVector(featureMat,cellAnnot$wellKey,FUN = mean)
corWell<-cor(featureMatPerWell)

colData = wellAnnot[rn(corWell),c("dataset","CellLine")]
colnames(colData)<-c("diffDay","cellLine")

colData<-oob::formatAnnotFromMeta(colData,metaAnnot = fastRead("zenodoArchive/CP_metaAnnot.tsv")[cn(colData),])

pdf("out2.2/perWellCor.pdf",width = 12,height = 10)
heatmap.DM(corWell,preSet = "cor",colData = colData, colorAnnot = colorScales[cn(colData)])
dev.off()


### feature Modules

# subsample cells to have the same number of cells per cluster
# to avoid bias in the correlation matrix

# group<-cellAnnot$cellClusters |> as.factor()
# names(group) <- rn(cellAnnot)
# drawSize <- min(table(group))
# drawCells <- unlist(lapply(unique(group), function(lvl) {
#     cells <- names(group)[group == lvl]
#     sample(cells, min(drawSize),length(cells))
# }))
drawCells<-scan("zenodoArchive/CP_drawedCells.txt",what = "character")

####


clustFeature<-hierarchicalClustering(featureMat[, drawCells],transpose = FALSE,method.dist = covDist)
best.cutree(clustFeature,graph = T)  # a good partition is 18 clusters
featureGroup<-cutree(clustFeature,k = 18)  |> as.factor()

featureGroupList<-factorToVectorList(featureGroup)
names(featureGroupList)<-c(
  "intensity.endoRet", #1
  "intensity.CGAN",  #2
  "entropy.nuc",   #3
  "intensity.mito",  #4
  "angMoment.mito",  #5
  "angMom.cytoStruct",   #6
  "angMom.nucleus",  #7
  "integrInt.cytoStruct",  #8
  "shape.nuc",  #9
  "area.nuc",  #10
  "area.cell",  #11
  "misc",  #12
  "shape.cell",  #13
  "infoMeas.CGAN",  #14
  "infoMeas.endoRet",  #15
  "massDisplace.mito",  #16
  "intensity.nuc",  #17
  "infoMeas.mito"   #18
)

featureGroup<-VectorListToFactor(featureGroupList)
featureGroup<-featureGroup[rn(featureMat)]

featureAnnot<-data.frame(row.names = names(featureGroup), module=featureGroup)
fastWrite(featureAnnot,"out2.2/featureAnnot.tsv")
##

pdf("out2.2/completeHt.pdf",width = 25,height = 20)
heatmap.DM(featureMat[, drawCells], colData = cellAnnot[drawCells, c("parentalCellLine","diffDay")],
           column_split = group[drawCells], row_split = featureGroup,colorAnnot = colorScales[c("parentalCellLine","diffDay")],
           cluster_row_slices = TRUE, cluster_column_slices = FALSE, show_column_names = FALSE,
           row_title_rot=0) #,colorAnnot = colorSc)
dev.off()

#module activation score per cell
activScores<-activeScorePCAlist(featureMat,featureGroupList)
activScoreMat<-activScores$activScoreMat|>t()

pdf("out2.2/htModuleActivScore.pdf",width = 25,height = 10)
heatmap.DM(activScoreMat[, drawCells], colData = cellAnnot[drawCells, c("parentalCellLine","diffDay")],
           column_split = group[drawCells],colorAnnot = colorScales[c("parentalCellLine","diffDay")],
           cluster_row_slices = TRUE, cluster_column_slices = FALSE, row_split = rn(activScoreMat), show_column_names = FALSE,preSet = "vanilla")
dev.off()

contribList<-unlist(activScores$contributionList,use.names = F)
names(contribList)<-lapply(activScores$contributionList,names) |> unlist()

featureAnnot$moduleContrib<-NA
featureAnnot[names(contribList),"moduleContrib"]<-contribList
featureAnnot$mean<-apply(featureMat,1,mean)

fastWrite(activScoreMat, "out2.2/activScoreMat.tsv")
fastWrite(featureAnnot, "out2.2/featureAnnot.tsv")

###
corFeatures<-cor(featureMat |> t())

pdf("out2.2/featureCorHeatmap.pdf",width = 14,height = 12)
Heatmap(corFeatures, row_split = featureGroup, column_split = featureGroup,col=c("darkblue", "white", "#FFAA00"),border = TRUE,row_title_rot = 0)
dev.off()

###


