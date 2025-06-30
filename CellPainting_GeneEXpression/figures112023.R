setwd(".")

library(oob)
library(ggplot2)
library(ggrepel)
library(scattermore)
library(glmnet)

source("functions.R")

gexAnnot<-formatAnnotFromMeta(fastRead("../data/gexAnnot.tsv"),
  fastRead("../data/metaAnnotGex.tsv"))

cpAnnot<-formatAnnotFromMeta(fastRead("../data/cpAnnot.tsv"), 
  fastRead("../data/metaAnnotCP.tsv"))

#gexData<-fastRead("../data/gexData.tsv",as.matrix=TRUE)
cpData<-fastRead("../data/cpData.tsv",as.matrix=TRUE)
cpFromGexTheoMat<-fastRead("../analysePredMatV2/outText/cpFromGexTheoMat.tsv",
  as.matrix=TRUE)[rn(gexAnnot),] 

colorScalesCP <-attr(cpAnnot,"colorScale")

gexAnnot<-gexAnnot[gexAnnot$finalAnnotation!="",]

commonCellLine<-intersect(unique(gexAnnot$parentalCellLine),
  unique(cpAnnot$parentalCellLine))

gexCellForEnrich<-gexAnnot$newAnnotCellType!="" & gexAnnot$parentalCellLine %in% commonCellLine
cpCellForEnrich<-cpAnnot$parentalCellLine %in% commonCellLine


enrichGex<-enrich2vector(gexAnnot$newAnnotCellType[gexCellForEnrich],
  gexAnnot$parentalCellLine[gexCellForEnrich])

enrichCP<-enrich2vector(cpAnnot$cellClusters[cpCellForEnrich],
  cpAnnot$parentalCellLine[cpCellForEnrich])

log2OEmatGex<-enrichGex$log2OEmat
log2OEmatGex[log2OEmatGex==-Inf]<- -7

log2OEmatCp<-enrichCP$log2OEmat

ggDataGex<-reshape2::melt(log2OEmatGex, varnames=c("cellGroup","cellLine"),value.name = "log2OE")
ggDataCp<- reshape2::melt(log2OEmatCp,  varnames=c("cellGroup","cellLine"),value.name = "log2OE")
ggDataGex$modality<-"gex"
ggDataCp$modality<-"cp"

ggData<-rbind(ggDataGex,ggDataCp)
ggData$group<-paste0(ggData$cellLine,ggData$modality)

ggplot(ggData,mapping=aes(x=log2OE,y=group,label=cellGroup))+
  geom_point()+
  geom_text_repel()+
  geom_hline(yintercept = seq(2.5,10.5,2))+
  geom_vline(xintercept = 0)+
  theme_bw()+
  xlim(c(-7,4))

ggsave("./outPlots/enrichCP_Gex_cellLine.pdf",  width=12,height=10)

###direct correlation

#cluster cpdata further
cpMitoNeg<-cpData[cpAnnot$cellClusters=="d70.mitoNeg",]

clusteringMitoNeg<-hierarchicalClustering(cpMitoNeg |> t())
best.cutree(clusteringMitoNeg, graph=TRUE)
clusterMitoNeg<-paste0("d70.mitoNeg.k",cutree(clusteringMitoNeg,k=4))
names(clusterMitoNeg)<-rn(cpMitoNeg)

cpAnnot$clustersPlus<-cpAnnot$cellClusters |> as.character()
cpAnnot[names(clusterMitoNeg),"clustersPlus"]<-clusterMitoNeg


metaCP<-oob::aggregMatPerVector(cpData, cpAnnot$clustersPlus)
metaCPfromGex<-oob::aggregMatPerVector(cpFromGexTheoMat, gexAnnot$finalAnnotation)

featMarkerOfCellType<-getMarkers(cpFromGexTheoMat |> t(),gexAnnot$finalAnnotation)
featMarkerOfCP<-getMarkers(cpData |> t(),cpAnnot$cellClusters)
fastWrite(featMarkerOfCellType,"./outText/featMarkerOfCellType.tsv")

featMarkerOfCellTypeExtracted<-extractFeatureMarkerData(featMarkerOfCellType)
featMarkerOfCPExtracted<-extractFeatureMarkerData(featMarkerOfCP)

featMarkerOfCellTypeExtracted<-featMarkerOfCellTypeExtracted[,cn(featMarkerOfCellTypeExtracted)!="NA."]

directCor<-cor(metaCP |> t(),metaCPfromGex |> t()) |> t()

pdf("./outPlots/heatmapDirectCorrelation.pdf",width=8,height=8)
heatmap.DM(directCor,preSet="cor",
  clustering_distance_rows=covDist, clustering_distance_columns=covDist,
  colorScale = c("darkblue", "lightblue", "white", "#FFAA00","red3"))
dev.off()
#by Feature

featureCPmodule<-fastRead("../data/featurePerModule.tsv")

featMarkerOfCellType<-getMarkers(cpFromGexTheoMat |> t(),gexAnnot$finalAnnotation)
featMarkerOfCP<-getMarkers(cpData |> t(),cpAnnot$clustersPlus)
#fastWrite(featMarkerOfCellType,"./outText/featMarkerOfCellType.tsv")

featMarkerOfCellTypeExtracted<-extractFeatureMarkerData(featMarkerOfCellType)
featMarkerOfCPExtracted<-extractFeatureMarkerData(featMarkerOfCP)

corFeatureScore<-cor(featMarkerOfCellTypeExtracted,featMarkerOfCPExtracted)

pdf("./outPlots/heatmapMarkerScoreCorrelation.pdf",width=8,height=8)
heatmap.DM(corFeatureScore,preSet="cor",
  clustering_distance_rows=covDist, clustering_distance_columns=covDist,
    colorScale = c("darkblue", "lightblue", "white", "#FFAA00","red3"))
dev.off()


htclusteringCP<-hierarchicalClustering(corFeatureScore)
htclusteringGex<-hierarchicalClustering(corFeatureScore, transpose=FALSE)

pdf("./outPlots/compareCorCPclust_CTypesHt.pdf",width=11,height=8)
print(
  heatmap.DM(corFeatureScore,preSet="cor",returnHeatmap=TRUE,name="cor\n(marker scores)",
    cluster_columns=htclusteringCP, cluster_rows=htclusteringGex,
      colorScale = c("darkblue", "lightblue", "white", "#FFAA00","red3"))+
    heatmap.DM(directCor,preSet="cor",returnHeatmap=TRUE,
    cluster_columns=htclusteringCP,name="cor\n(direct)",
    colorScale = c("darkblue", "lightblue", "white", "#FFAA00","red3"))
)
dev.off()

#cell cycle

#gexFromCpTheoMat<-fastRead("../analysePredMatV2/outText/gexFromCpTheoMat.tsv")

cor(cpFromGexTheoMat,gexAnnot$S.Score )
cor(cpFromGexTheoMat,gexAnnot$G2M.Score )

cellCycleScore <- gexAnnot$S.Score - gexAnnot$G2M.Score
corOfFeatureToScore <- cor(cpFromGexTheoMat,cellCycleScore)


correlatedFeatToCCscore <- cn(cpFromGexTheoMat)[abs(corOfFeatureToScore)>0.5]

lmData <- data.frame(ccScore = cellCycleScore, cpFromGexTheoMat[,correlatedFeatToCCscore])

scattermoreplot(cellCycleScore,cpFromGexTheoMat[,"Nuc_Txtr_InfoMeas1_Nuc_3_01_256"])
scattermoreplot(cellCycleScore,cpFromGexTheoMat[,"Cyt_Int_MeanInt_EndoRet"])
scattermoreplot(cellCycleScore,cpFromGexTheoMat[,"Cel_Shape_Eccentricity"])

lmRes<- lm( ccScore~., data=lmData)

summary(lmRes)

makeLassoModel<-function(varToExplain,featureMatrix){
  trainingSamples<-sample(c(TRUE,FALSE),nrow(featureMatrix) ,replace=TRUE,prob=c(2/3,1/3))
  cv_model <- cv.glmnet(featureMatrix[trainingSamples,],varToExplain[trainingSamples], alpha = 1)
  best_lambda <- cv_model$lambda.min
  best_model <- glmnet(featureMatrix[trainingSamples,],varToExplain[trainingSamples], alpha = 1, lambda = best_lambda)

  y_predicted <- predict(best_model, s = best_lambda, newx = featureMatrix[!trainingSamples,])
  y = varToExplain[!trainingSamples]
  #find SST and SSE
  sst <- sum((y - mean(y))^2)
  sse <- sum((y_predicted - y)^2)

  #find R-Squared
  rsq <- 1 - sse/sst
  list(best_model = best_model, best_lambda=best_lambda, rsq=rsq)
}

modelCCscore <- makeLassoModel(cellCycleScore,cpFromGexTheoMat)
modelG2M<- makeLassoModel(gexAnnot$G2M.Score,cpFromGexTheoMat)
modelS <- makeLassoModel(gexAnnot$S.Score,cpFromGexTheoMat)

cp_G2Mscore <- predict(modelG2M$best_model, s = modelG2M$best_lambda, newx = cpData)[,1]
cp_Sscore <- predict(modelS$best_model, s = modelS$best_lambda, newx = cpData)[,1]
cp_CCscore <- cp_Sscore - cp_G2Mscore

dtCCCP <- data.frame(row.names=rn(cpAnnot), G2M.Score = cp_G2Mscore, S.Score=cp_Sscore, delta = cp_CCscore)
dtCCCP$Phase <- ""

dtCCCP[dtCCCP$S.Score < 0 & dtCCCP$G2M.Score < 0 , "Phase"] <- "G1"
dtCCCP[dtCCCP$delta < 0 & dtCCCP$Phase !=  "G1" , "Phase"] <- "G2M"
dtCCCP[dtCCCP$delta >= 0 & dtCCCP$Phase !=  "G1" , "Phase"] <- "S"

table(dtCCCP$Phase)
table(gexAnnot$Phase)

scattermoreplot(cp_CCscore, predict(modelCCscore$best_model, s = modelCCscore$best_lambda, newx = cpData)[,1])
qplotDensity(cp_G2Mscore)
qplotDensity(gexAnnot$G2M.Score)


ccdf<-data.frame(ccScore = cellCycleScore,  phase=gexAnnot$Phase,g2m=gexAnnot$G2M.Score, s=gexAnnot$S.Score)
ggplot( ccdf, aes(x= phase, y = g2m))+
  geom_scattermore()

ggplot( ccdf[ccdf$phase!="G1",], aes(x= phase, y = s-g2m))+
  geom_scattermore()

table(ccdf$phase, ccdf$g2m<0 & ccdf$s>0)

cpProj<-fastRead("../data/cpProj.tsv",as.matrix=T)
proj2d(cpProj, colorBy = dtCCCP$Phase,useScatterMore=TRUE, returnGraph=TRUE) + facet_wrap(~dtCCCP$Phase)
proj2d(cpProj, colorBy = dtCCCP$G2M.Score,useScatterMore=TRUE,  colorScale=c("darkblue","blue","grey75","orange","red"))
proj2d(cpProj, colorBy = dtCCCP$S.Score,useScatterMore=TRUE,  colorScale=c("darkblue","blue","grey75","orange","red"))
proj2d(cpProj, colorBy = dtCCCP$delta,useScatterMore=TRUE,  colorScale=c("darkblue","blue","grey75","orange","red"))
proj2d(cpProj, colorBy = dtCCCP$G2M.Score/dtCCCP$S.Score,useScatterMore=TRUE,  colorScale=c("darkblue","blue","grey75","orange","red"))

dtCCCP$cellClusters <- cpAnnot$cellClusters

ggplot(dtCCCP, aes(x=cellClusters,y=G2M.Score,fill=cellClusters))+
  geom_violin(scale="width")+
  geom_boxplot(fill="black",color="white",width=.2,outlier.alpha=0)+
  scale_fill_manual(values=colorScalesCP$cellClusters)+
  theme_bw()

ggsave("./outPlots/G2MscorePerCPcluster.pdf",width=8,height=6)

ggplot(dtCCCP, aes(x=cellClusters,y=S.Score,fill=cellClusters))+
  geom_violin(scale="width")+
  geom_boxplot(fill="black",color="white",width=.2,outlier.alpha=0)+
  scale_fill_manual(values=colorScalesCP$cellClusters)+
  theme_bw()

ggsave("./outPlots/SscorePerCPcluster.pdf",width=8,height=6)

ggplot(dtCCCP[dtCCCP$cellClusters=="d20.nucNeg",], aes(x=G2M.Score,y=S.Score))+
  geom_scattermore()

####
gexFromCpTheoMat<-fastRead("../analysePredMatV2/outText/gexFromCpTheoMat.tsv")

panRGmarkerDt<-fastRead("./inputs/markerOfpanRG.tsv")
panRGmarkerDt<-panRGmarkerDt[panRGmarkerDt$CellType=="panRG-O",]
panRGmarkerDt<-panRGmarkerDt[order(panRGmarkerDt$avg_log2FC,decreasing=T),]

qplotDensity(panRGmarkerDt$avg_log2FC)

bestmarker <- panRGmarkerDt$feature[panRGmarkerDt$avg_log2FC>.5] |> unique()

bestmarker<-intersect(bestmarker, cn(gexFromCpTheoMat))


panRGmark_AS<-activScorePCA(t(gexFromCpTheoMat),bestmarker,returnContribution=FALSE)

pdf("./outPlots/projPanRGmarker.pdf",width=10,height=7)
proj2d(cpProj,colorBy=panRGmark_AS, useScatterMore=TRUE, 
  colorScale=c("darkblue","blue","cyan","lightgreen","orange","red","darkred"))
dev.off()



proj2d(cpProj,colorBy=gexFromCpTheoMat[,"SFRP2"], useScatterMore=TRUE, 
colorScale=c("darkblue","blue","cyan","lightgreen","orange","red","darkred"))

