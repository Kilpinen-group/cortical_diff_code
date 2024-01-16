setwd(".")

library(oob)
library(ggplot2)
library(DBI)
library(plotly)
library(igraph)
library(ggbeeswarm)
library(ggrepel)
library(stringr)
library(gridExtra)
library(patchwork)
library(ComplexHeatmap)

cellAnnot<- fastRead("input/cellAnnotFiltered.tsv")
cellAnnot$parentalCellLine<-paste0("L",cellAnnot$parentalCellLine)
cellAnnot$cellLine<-paste0("L",str_replace(cellAnnot$cellLine,pattern = "_",replacement = "R"))
cellAnnot$expPop<-paste0(cellAnnot$diffDay,cellAnnot$parentalCellLine)
cellAnnot$expPopWithReplicate<-paste0(cellAnnot$diffDay,cellAnnot$cellLine)


cellAnnot<-oob::formatAnnotFromMeta(annotDataFrame = cellAnnot,metaAnnot = fastRead("input/metaAnnot.tsv"))

colorScales<-attr(cellAnnot,"colorScales")

featureMat<-fastRead("input/featureMatFiltered.tsv", as.matrix = TRUE)


PCAdat<-fastPCA(
  featureMat,
  nPC = 30,
  maxit = 1000
)

pcrRes<-PCR(PCAdat,cellAnnot[,c("diffDay","leidenClust","expPop","cellLine")],nComponent = 30)


ggBorderedFactors(
  ggplot(pcrRes,aes(x=PC,y=Rsquared,fill=Annotation))+
    geom_beeswarm(pch=21,size=4,cex = .5)+
    xlab("Principal component")+ylab("R²")+
    scale_fill_manual(values = oobColors(nlevels(pcrRes$Annotation)))+
    theme(
      panel.grid.major.y = element_line(colour = "grey75"),
      panel.grid.minor.y = element_line(colour = "grey75"),
      panel.background = element_rect(fill = NA,colour="black")
    )
)


umap<-UMAP(featureMat,ret_nn = T)
cellAnnot$leidenClust<-leidenFromUMAP(umap,n_neighbors = 20, resolution_parameter = 1)

fastWrite(cellAnnot,"outText/cellAnnot.tsv")

colnames(umap$embedding)<-c("V1","V2")
fastWrite(umap$embedding,"outText/umapCoord.tsv")

pdf("outPlot/projAnnot.pdf",width=13, height=12)
proj2d(umap, useScatterMore=T,colorBy = cellAnnot$diffDay,colorScale = colorScales$diffDay)
proj2d(umap, useScatterMore=T,colorBy = cellAnnot$leidenClust,plotFactorsCentroids = T)
proj2d(umap, useScatterMore=T,colorBy = cellAnnot$parentalCellLine,colorScale = colorScales$parentalCellLine)
proj2d(umap, useScatterMore=T,colorBy = cellAnnot$expPop,colorScale = colorScales$expPop)
dev.off()

experimentalCols<-c("diffDay" ,"parentalCellLine", "expPop", "leidenClust")
markerPerAnnot<-list()
for(annot in experimentalCols) markerPerAnnot[[annot]]<-getMarkers(featureMat,cellAnnot[,annot])

colorScalesHt<-colorScales[intersect(names(colorScales),experimentalCols)]
for(feature in experimentalCols){
  pdf(paste0("outPlot/markerHt_",feature,".pdf"),width = 20,height = 12)
  htMarker(featureMat,cellAnnot[,feature],markerData = extractFeatureMarkerData(markerPerAnnot[[feature]]),
           colData = cellAnnot[,experimentalCols,drop=F],maxDrawSize = 500,colorAnnot=colorScalesHt,
           name="feature value\n(positive markers)")
  htMarker(featureMat,cellAnnot[,feature],markerData = - extractFeatureMarkerData(markerPerAnnot[[feature]]),
           colData = cellAnnot[,experimentalCols,drop=F],maxDrawSize = 500,colorAnnot=colorScalesHt,
           name="feature value\n(negative markers)")
  dev.off()
}

for(feature in experimentalCols){
  fastWrite(markerPerAnnot[[feature]],paste0("outText/markers_",feature,".tsv"))
}

bestFeaturePerClust<-rn(markerPerAnnot$leidenClust)[extractFeatureMarkerData(markerPerAnnot$leidenClust) |> apply(2,which.max)] |> unique()

pdf("outPlot/violinMarke4Cluster.pdf",width = 30,height = 6)
print(
  plotExpr(expr = featureMat[bestFeaturePerClust,], group=cellAnnot$leidenClust |> as.factor(), log10Plus1yScale = F,colorScale =  oobColors(unique(cellAnnot$leidenClust) |> len()), returnGraph = T)
)
dev.off()

#bestFeatPerDiffPoint
bestFeatureDiffPoint<-rn(markerPerAnnot$diffDay)[extractFeatureMarkerData(markerPerAnnot$diffDay) |> apply(2,whichTop,top = 5)] |> unique()

pdf("outPlot/violinMarkerDiffPoint.pdf",width = 10,height = 6)
print(
  plotExpr(expr = featureMat[bestFeatureDiffPoint,], group=cellAnnot$diffDay, log10Plus1yScale = F,colorScale =  oobColors(3), returnGraph = T)
)
dev.off()

dir.create("outPlot/projPerFeature",showWarnings = F)
for(f in rn(featureMat)){
  pdf(paste0("outPlot/projPerFeature/",f,".pdf"),width=13, height=12)
  proj2d(umap, useScatterMore=T,colorBy = featureMat[f,],colorScale = c("darkblue","blue","cyan","lightgreen","orange","red","darkred"))
  dev.off()
}


#k08 vs k06
samplesK08<-rn(cellAnnot)[cellAnnot$leidenClust=="k08" & cellAnnot$diffDay=="d70"]
samplesK06<-rn(cellAnnot)[cellAnnot$leidenClust=="k06" & cellAnnot$diffDay=="d70"]

samplesK06<-sample(samplesK06,len(samplesK08))

sampleOfInterest<-c(samplesK06,samplesK08)

resLinearTest<- multiLinearModel(exprData = featureMat[,sampleOfInterest],colData = cellAnnot[sampleOfInterest,],contrast = c("leidenClust","k06","k08"))

cellsOfHeatmap<-sample(sampleOfInterest,500)

bestk06<-rn(resLinearTest)[oob::whichTop(resLinearTest$log2FoldChange,top = 10)]
bestk08<-rn(resLinearTest)[oob::whichTop(resLinearTest$log2FoldChange,top = 10,decreasing = FALSE)]

vectorIsBest<-c( rep("K06_marker",len(bestk06)), rep("K08_marker",len(bestk08))); names(vectorIsBest)<-c(bestk06,bestk08)

pdf("outPlot/diffValk06k08.pdf",width=12,height=8)
heatmap.DM(featureMat[names(vectorIsBest),cellsOfHeatmap],colData = cellAnnot[cellsOfHeatmap,experimentalCols,drop=F],colorAnnot=colorScalesHt,
           column_split = cellsOfHeatmap %in% samplesK06, row_split = vectorIsBest,
           cluster_row_slices = FALSE, cluster_column_slices = FALSE)
dev.off()

#export 2 web app

spatialInfoAnnot<-readRDS("../createFeatureMat/robj/spatialInfoAnnot.rds")[cn(featureMat),]

export2WebApp<-data.frame(x=umap$embedding[,1],y=umap$embedding[,2],diffDay=cellAnnot$diffDay,cluster=cellAnnot$leidenClust,
                          image_filename=paste0("day",cellAnnot$diffDay|>as.character() |> substr(2,3),"/imgsOut/",
                                                substr(cellAnnot$imageFileName,1,nchar(cellAnnot$imageFileName)-5),".png"),
                          bbox_minX=spatialInfoAnnot[,"Cel_Shape_BBoxMin_X"],
                          bbox_maxX=spatialInfoAnnot[,"Cel_Shape_BBoxMax_X"],
                          bbox_minY=spatialInfoAnnot[,"Cel_Shape_BBoxMin_Y"],
                          bbox_maxY=spatialInfoAnnot[,"Cel_Shape_BBoxMax_Y"],
                          centerX = spatialInfoAnnot[,"Nuc_Location_Center_X"],
                          centerY= spatialInfoAnnot[,"Nuc_Location_Center_Y"],
                          t(featureMat))
fastWrite(export2WebApp,"D:/PostdocUnsync/05_cellPainting/webapp/data/webAppFeatures.tsv")

###enrichment of cell line
source("function.R")


enrichResParental<-enrich2vector(cellAnnot$cellClusters,cellAnnot$parentalCellLine)
enrichResDiffDay <-enrich2vector(cellAnnot$cellClusters,cellAnnot$diffDay)

ggOE<-reshape2::melt(cbind(enrichResDiffDay$log2OEmat,enrichResParental$log2OEmat), value.name = "log2OE",varnames=c("cluster","attr"))

pvals<-cbind(enrichResDiffDay$pvalMat,enrichResParental$pvalMat) |>  data.frame() |> unlist()
ggOE$pvals<-pvals
ggOE$padj<-p.adjust(ggOE$pvals,method = "BH")
ggOE$isSignif<-ggOE$padj<0.01

ggplot(ggOE,aes(x=cluster,y=log2OE,label=attr,fill=attr,shape=isSignif))+
  geom_point(size=3,color="black")+
  geom_text_repel()+
  theme_bw()+
  scale_fill_manual(values = colorScales$cellClusters)+
  scale_shape_manual(values = c(21,23))+
  ylim(-10,10)
ggsave("outPlot/OE.pdf",width = 16,height = 9)


binLeiden<-sapply(levels(cellAnnot$cellClusters), function(x) x == cellAnnot$cellClusters); colnames(binLeiden)<-levels(cellAnnot$cellClusters)
binParentalLine<-sapply(levels(cellAnnot$parentalCellLine), function(x) x == cellAnnot$parentalCellLine); colnames(binParentalLine)<-levels(cellAnnot$parentalCellLine)
binDiffDay<- sapply(levels(cellAnnot$diffDay |> as.factor()), function(x) x == cellAnnot$diffDay); colnames(binDiffDay)<-levels(cellAnnot$diffDay |> as.factor())

binsExpVar<-cbind(binDiffDay,binParentalLine)

AMIs<-matrixFromDimnames(colnames(binLeiden),colnames(binsExpVar)) #oob function: create an empty matrix from 2 vector with row/colnames
for(k in rownames(AMIs)){
  for( l in colnames(AMIs)){
    AMIs[k,l]<-aricode::AMI(binLeiden[,k],binsExpVar[,l])
  }
}

ggAMIs<-reshape2::melt(AMIs,value.name = "AMI",varnames=c("cluster","attr"))
ggAMIs$attr<-as.factor(ggAMIs$attr)


ggAll<-ggOE
ggAll$AMI<-ggAMIs$AMI
ggAll$attrType<-"Line"
ggAll[grep("^d",ggAll$attr),"attrType"]<-"diffDay"

ggplot(ggAll,aes(y=AMI,x=log2OE,label=paste0(cluster,"_",attr),fill=cluster))+
  geom_point(shape=21,size=3)+geom_text_repel()+
  theme_bw()+facet_wrap(~attrType)+
  scale_fill_manual(values=colorScales$cellClusters)
ggsave("outPlot/OEvsAMI.pdf",width = 16,height = 7)

ggplot(ggAll[ggAll$attrType=="Line",],aes(y=AMI,x=log2OE,label=paste0(cluster,"_",attr),fill=cluster))+
  geom_point(shape=21,size=3)+geom_text_repel()+
  theme_bw()+facet_wrap(~attrType)+
  scale_fill_manual(values=colorScales$cellClusters)
ggsave("outPlot/OEvsAMI_lineOnly.pdf",width = 8,height = 7)

#simple PCA
pca<-fastPCA(featureMat)
pca2d(pca,useScatterMore=TRUE,colorBy = cellAnnot$diffDay)

###total heatmap

group<-cellAnnot$leidenClust |> as.factor()
names(group) <- rn(cellAnnot)
drawSize <- min(table(group))
drawCells <- unlist(lapply(unique(group), function(lvl) {
  cells <- names(group)[group == lvl]
  sample(cells, min(drawSize),length(cells))
}))
group[drawCells]

clustFeature<-hierarchicalClustering(featureMat[, drawCells],transpose = FALSE,method.dist = covDist)
best.cutree(clustFeature,graph = T)
featureGroup<-cutree(clustFeature,k = 12)  |> as.factor()

pdf("outPlot/completeHt.pdf",width = 25,height = 20)
heatmap.DM(featureMat[, drawCells], colData = cellAnnot[drawCells, c("parentalCellLine","diffDay")],
           column_split = group[drawCells], row_split = featureGroup,colorAnnot = colorScales[c("parentalCellLine","diffDay")],
           cluster_row_slices = TRUE, cluster_column_slices = TRUE, show_column_names = FALSE) #,colorAnnot = colorSc)
dev.off()

#per well stats

wellAnnot<-fastRead("../createFeatureMat/in/wellAnnot.tsv")
cellAnnot$wellKey<-paste0(cellAnnot$diffDay,substr(cellAnnot$well,1,1),formatNumber2Character(substr(cellAnnot$well,2,4)))

featureMatPerWell<-aggregMatPerVector(featureMat,cellAnnot$wellKey,FUN = mean)
corWell<-cor(featureMatPerWell)


old2newCellLineName<-function(x){
  paste0("L",substr(x,1,4),"R",substr(x,6,9))
}

colData = wellAnnot[rn(corWell),c("dataset","CellLine")]
colnames(colData)<-c("diffDay","cellLine")
colData$cellLine<-old2newCellLineName(colData$cellLine)

colData<-oob::formatAnnotFromMeta(colData,metaAnnot = fastRead("input/metaAnnot.tsv")[cn(colData),])

pdf("outPlot/perWellCor.pdf",width = 12,height = 10)
heatmap.DM(corWell,preSet = "cor",colData = colData, colorAnnot = colorScales[cn(colData)])
dev.off()


#feature Modules
corFeatures <- cor(featureMat[,drawCells]|>t())

pdf("outPlot/correlationHeatmap.pdf",width=12,height=11)
heatmap.DM(corFeatures,preSet="cor")
dev.off()

vectCorFeature <- corFeatures[upper.tri(corFeatures)]

corVals<-seq(0.5,0.9,0.01)
nBCluster<-rep(0,len(corVals))

for(i in seq_along(corVals)){
  corVal<-corVals[i]
  graphOfFeature <-graph_from_adjacency_matrix(corFeatures > corVal, mode = "upper", diag = F)


  c1 = cluster_fast_greedy(graphOfFeature)
  nBCluster[i]<-max(c1$membership)
}

ggplot(data.frame(corVals,nBCluster),aes(x=corVals, y=nBCluster))+
  geom_point()

graphOfFeature <-graph_from_adjacency_matrix(corFeatures > 0.82, mode = "upper", diag = F)
c1 = cluster_fast_greedy(graphOfFeature)


pdf("outPlot/moduleOfFeatures.pdf",width = 7,height = 7)
plot(c1,
     graphOfFeature,
     vertex.label = NA,
     vertex.size = 2)
dev.off()

graphStat <-
  data.frame(
    centrality = centr_degree(graphOfFeature)$res,
    mod = c1$membership,
    row.names = names(V(graphOfFeature))
  )
modules <- unique(c1$membership)

featurePerModule<-c1$membership
names(featurePerModule)<-c1$names
featurePerModule<-factorToVectorList(as.factor(featurePerModule))
featurePerModule<-featurePerModule[sapply(featurePerModule,len)>4]


pdf("outPlot/completeHtPerGraphModule.pdf",width = 25,height = 20)
heatmap.DM(featureMat[c1$names, drawCells], colData = cellAnnot[drawCells, c("parentalCellLine","diffDay")],
           column_split = group[drawCells], row_split = c1$membership,colorAnnot = colorScales[c("parentalCellLine","diffDay")],
           cluster_row_slices = TRUE, cluster_column_slices = TRUE, show_column_names = FALSE)
dev.off()


####


clustFeature<-hierarchicalClustering(featureMat[, drawCells],transpose = FALSE,method.dist = covDist)
best.cutree(clustFeature,graph = T)
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
fastWrite(featureAnnot,"outText/featureAnnot.tsv")
##
cellClusters<-rep("",nrow(cellAnnot));names(cellClusters)<-rn(cellAnnot)

cellClusters[rn(cellAnnot)[cellAnnot$leidenClust=="k01"]]<-"d40.smallCells"
cellClusters[rn(cellAnnot)[cellAnnot$leidenClust=="k02"]]<-"d20.mature"
cellClusters[rn(cellAnnot)[cellAnnot$leidenClust=="k03"]]<-"d40.endoRet"
cellClusters[rn(cellAnnot)[cellAnnot$leidenClust=="k04"]]<-"d20.endoRetNeg"
cellClusters[rn(cellAnnot)[cellAnnot$leidenClust=="k05"]]<-"d40.mature"
cellClusters[rn(cellAnnot)[cellAnnot$leidenClust=="k06"]]<-"d70.mitoNeg"
cellClusters[rn(cellAnnot)[cellAnnot$leidenClust=="k07"]]<-"d20.mito"
cellClusters[rn(cellAnnot)[cellAnnot$leidenClust=="k08"]]<-"d70_20.mito"
cellClusters[rn(cellAnnot)[cellAnnot$leidenClust=="k09"]]<-"d40.early"
cellClusters[rn(cellAnnot)[cellAnnot$leidenClust=="k10"]]<-"d70.bigcells"
cellClusters[rn(cellAnnot)[cellAnnot$leidenClust=="k11"]]<-"d20.nucNeg"


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

cellAnnot$cellClusters<-factor(cellClusters,levels = names(colorScales$cellClusters))

pdf("outPlot/projCellCLusters.pdf",width=13, height=12)
proj2d(umap, useScatterMore=T,colorBy = cellAnnot$cellClusters,plotFactorsCentroids = T,colorScale = colorScales$cellClusters)
dev.off()

group<-cellAnnot$cellClusters;names(group)<-rn(cellAnnot)

fastWrite(cellAnnot,"outPlot/cellAnnot.tsv")


##

pdf("outPlot/completeHt.pdf",width = 25,height = 20)
heatmap.DM(featureMat[, drawCells], colData = cellAnnot[drawCells, c("parentalCellLine","diffDay")],
           column_split = group[drawCells], row_split = featureGroup,colorAnnot = colorScales[c("parentalCellLine","diffDay")],
           cluster_row_slices = TRUE, cluster_column_slices = FALSE, show_column_names = FALSE,
           row_title_rot=0) #,colorAnnot = colorSc)
dev.off()


activScores<-activeScorePCAlist(featureMat,featureGroupList)
activScoreMat<-activScores$activScoreMat|>t()

pdf("outPlot/htModuleActivScore.pdf",width = 25,height = 10)
heatmap.DM(activScoreMat[, drawCells], colData = cellAnnot[drawCells, c("parentalCellLine","diffDay")],
           column_split = group[drawCells],colorAnnot = colorScales[c("parentalCellLine","diffDay")],
           cluster_row_slices = TRUE, cluster_column_slices = FALSE, row_split = rn(activScoreMat), show_column_names = FALSE,preSet = "vanilla")
dev.off()

contribList<-unlist(activScores$contributionList,use.names = F)
names(contribList)<-lapply(activScores$contributionList,names) |> unlist()

featureAnnot$moduleContrib<-NA
featureAnnot[names(contribList),"moduleContrib"]<-contribList
featureAnnot$mean<-apply(featureMat,1,mean)

fastWrite(activScoreMat, "outText/activScoreMat.tsv")
fastWrite(featureAnnot, "outText/featureAnnot.tsv")
####

clusterPerImage<-sapply(unique(cellAnnot$imageKey),function(image){
  cellOfImage<-rn(cellAnnot)[cellAnnot$imageKey==image]
  cellAnnot[cellAnnot$imageKey==image,"cellClusters"] |> table()
}) |> t()

clusterPerImage<-clusterPerImage[rowSums(clusterPerImage)>=10,]

clusterPerImage_clustering<-hierarchicalClustering(clusterPerImage,method.dist = "pearson",transpose = FALSE)
plot(clusterPerImage_clustering,hang=-1)

best.cutree(clusterPerImage_clustering,graph = T)
clusters<-cutree(clusterPerImage_clustering,k = 10)
clusters<-paste0("imgk", formatNumber2Character(clusters))

clusterPerImage_prop<-clusterPerImage/rowSums(clusterPerImage)*100
clustersImg<-unique(clusters)|>sort()

dday<-substr(rn(clusterPerImage),1,3)
ddayPerCluster<-table(data.frame(clusters,dday))

linesPerExpPop<-rep("",nrow(clusterPerImage));names(linesPerExpPop)<-rn(clusterPerImage)
for(img in rn(clusterPerImage)) linesPerExpPop[img]<-cellAnnot[cellAnnot$imageKey==img,"expPop"][1] |> as.character()
linesPerExpPopPerCluster<-table(data.frame(clusters,linesPerExpPop))

linesPerParCellLine<-rep("",nrow(clusterPerImage));names(linesPerParCellLine)<-rn(clusterPerImage)
for(img in rn(clusterPerImage)) linesPerParCellLine[img]<-cellAnnot[cellAnnot$imageKey==img,"parentalCellLine"][1] |> as.character()
linesPerParentalLinePerCluster<-table(data.frame(clusters,linesPerParCellLine))

countPerCluster<-aggregMatPerVector(clusterPerImage,clusters,sum)

ggDataKimg<-reshape2::melt(countPerCluster,value.name = "count",varnames=c("imgCluster","CPcluster"))
ggDataDday<-reshape2::melt(ddayPerCluster,value.name = "count",varnames=c("imgCluster","diffDay"))
ggDataExpPop<-reshape2::melt(linesPerExpPopPerCluster,value.name = "count",varnames=c("imgCluster","expPop"))
ggDataParentalLine<-reshape2::melt(linesPerParentalLinePerCluster,value.name = "count",varnames=c("imgCluster","cellLine"))


g1<-ggplot(ggDataKimg,mapping=aes(x=imgCluster,fill=CPcluster,y=count))+
  geom_bar(position="stack", stat="identity",color="black")+
  scale_fill_manual(values=colorScales$cellClusters)
g2<-ggplot(ggDataDday,aes(x=imgCluster,fill=diffDay,y=count))+
  geom_bar(position="stack", stat="identity",color="black")+
  scale_fill_manual(values=colorScales$diffDay)
g3<-ggplot(ggDataExpPop,aes(x=imgCluster,fill=expPop,y=count))+
  geom_bar(position="stack", stat="identity",color="black")+
  scale_fill_manual(values=colorScales$expPop)
g4<-ggplot(ggDataParentalLine,aes(x=imgCluster,fill=cellLine,y=count))+
  geom_bar(position="stack", stat="identity",color="black")+
  scale_fill_manual(values=colorScales$parentalCellLine)

pdf("outPlot/imgClusterDistrib.pdf",width = 7,height = 14)
print(g1 / g2 / g4 /g3)
dev.off()

###

umapRes<-UMAP(activScoreMat[rn(activScoreMat)!="misc",],n_neighbors = 10)

proj2d(umapRes,colorBy = cellAnnot$cellClusters,colorScale = colorScales$cellClusters,
       useScatterMore = T)

trimapRes<-TRIMAP(activScoreMat[rn(activScoreMat)!="misc",])
proj2d(trimapRes,colorBy = cellAnnot$cellClusters,colorScale = colorScales$cellClusters,
       useScatterMore = T)
proj2d(trimapRes,colorBy = cellAnnot$diffDay,colorScale = colorScales$diffDay,
       useScatterMore = T)




###co-occurance matrix

clusterPerImage<-sapply(unique(cellAnnot$imageKey),function(image){
  cellOfImage<-rn(cellAnnot)[cellAnnot$imageKey==image]
  cellAnnot[cellAnnot$imageKey==image,"cellClusters"] |> table()
}) |> t()

CPclusters<-levels(cellAnnot$cellClusters)

cooccMat<-sapply(CPclusters,function(cluster){
  picturesWithCluster<-unique(cellAnnot[cellAnnot$cellClusters==cluster,"imageKey"])
  tempMat<-clusterPerImage[picturesWithCluster,]
  tempMat[tempMat>1]<-1
  colSums(tempMat)/max(colSums(tempMat))
})

library(grid)
library(ComplexHeatmap)

Heatmap(cooccMat,
          rect_gp = gpar(type = "none"),
          cluster_rows = FALSE, cluster_columns = FALSE,
  cell_fun = function(j, i, x, y, w, h, fill) {
    if(i > j) {
     grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
    }
  },
  col=circlize::colorRamp2(c(0,0.1,0.25,0.5,1),colors = c("white","yellow","orange","red","purple"))
)

#make mat by diff day
table(cellAnnot$diffDay,cellAnnot$cellClusters)

#per well analysis
wellAnnot$CellLine<-paste0("L",str_replace(wellAnnot$CellLine,"_","R"))
wellAnnot$parentalCellLine<-substr(wellAnnot$CellLine,1,5)

cellAnnot$well_day<-paste0(cellAnnot$diffDay,substr(cellAnnot$well,1,1),oob::formatNumber2Character(substr(cellAnnot$well,2,3)))
wellClustTable<-sapply(unique(cellAnnot$well_day),function(well){
  cellOfImage<-rn(cellAnnot)[cellAnnot$well_day==well]
  cellAnnot[cellAnnot$well_day==well,"cellClusters"] |> table()
}) |> t()

rn(wellAnnot)[wellAnnot$CellLine=="L11.4R4.2"]
wellAnnot<-wellAnnot[rn(wellClustTable),] #only keep well with data
wellAnnot$diffDay<-wellAnnot$dataset
fastWrite(wellAnnot,"outText/wellAnnot.tsv")


for(day in levels(cellAnnot$diffDay)){
    wells<-rn(wellAnnot)[wellAnnot$diffDay==day]
    wellClustTableOfDay<-wellClustTable[wells,]
    ggWell<-reshape2::melt(wellClustTableOfDay,value.name = "count",varnames=c("well","CPcluster"))
    ggWell$CPcluster = factor(ggWell$CPcluster ,levels = levels(cellAnnot$cellClusters))
    ggWell$CellLine<-wellAnnot[ggWell$well |> as.character(),"CellLine"]
    ggWell$well<-substr(ggWell$well,4,6)
    ggplot(ggWell,mapping=aes(x=well,fill=CPcluster,y=count))+
      geom_bar(position="stack", stat="identity",color="black")+
      scale_fill_manual(values=colorScales$cellClusters)+
      facet_wrap(~CellLine,scales="free_x")
    ggsave(filename = paste0("outPlot/wellplot",day,".pdf"),width = 12,height = 12)

    }


library(EBImage)
dat<-"d70"
wells<-rn(wellAnnot)[wellAnnot$diffDay==day]
wellClustTableOfDay<-wellClustTable[wells,]
ggWell<-reshape2::melt(wellClustTableOfDay,value.name = "count",varnames=c("well","CPcluster"))
ggWell$CPcluster = factor(ggWell$CPcluster ,levels = levels(cellAnnot$cellClusters))
ggWell$CellLine<-wellAnnot[ggWell$well |> as.character(),"CellLine"]
ggWell$well<-substr(ggWell$well,4,6)

findImageName<-function(day,wells){
  wellTrim<-sapply(wells,function(well){
    if(substr(well,2,2)=="0"){
      paste0(substr(well,1,1),substr(well,3,3))
    }else{
      well
    }
  })
  dayFull<-paste0("day",substr(day,2,4))
  return(paste0("D:/PostdocUnsync/05_cellPainting/outBIAS/",dayFull,"/concatWell/",wellTrim,".tiff"))
}

ggWell$imgName<-findImageName(day,ggWell$well)

##
markerPerAnnot[["cellClusters"]]<-getMarkers(featureMat,cellAnnot$cellClusters)
fastWrite(markerPerAnnot[["cellClusters"]],paste0("outText/markers_","cellClusters",".tsv"))


#feature correlation Heatmap
corFeatures<-cor(featureMat |> t())

pdf("outPlot/featureCorHeatmap.pdf",width = 14,height = 12)
Heatmap(corFeatures, row_split = featureGroup, column_split = featureGroup,col=c("darkblue", "white", "#FFAA00"),border = TRUE,row_title_rot = 0)
dev.off()



#test umap
cellsOfd20_70<- rn(cellAnnot)[cellAnnot$cellClusters == "d70_20.mito"]
umapShape<-UMAP(featureMat[grep("Shape",rn(featureMat),value = T),cellsOfd20_70])

proj2d(umapShape, colorBy=cellAnnot[cellsOfd20_70,]$diffDay)

d20_70dat <- featureMat[,cellsOfd20_70] |> t() |> data.frame()
d20_70dat$minMajAxisRatio <- d20_70dat$Cel_Shape_MajorAxisLength - d20_70dat$Cel_Shape_MinorAxisLength

d20_70dat$diffDay<-cellAnnot[cellsOfd20_70,]$diffDay

ggplot(d20_70dat, aes(x=diffDay,y=minMajAxisRatio))+
  geom_violin()

pdf("outPlot/clusterPerDay.pdf",width = 14,height = 8)
print(
  proj2d(umap, colorBy=cellAnnot$cellClusters,returnGraph = TRUE,colorScale = colorScales$cellClusters, useScatterMore = T) + facet_wrap(cellAnnot$diffDay)
)
dev.off()