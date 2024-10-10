setwd(".") #change to your working dir

library(oob)
library(ggplot2)
library(DBI)
library(plotly)
library(igraph)

set.seed(666)

inputDir <- "zenodoArchive/"
dir.create("out2.1", showWarnings = FALSE)

wellAnnot <- fastRead(paste0(inputDir,"CP_wellAnnot.tsv"))

conDay20 <-
  dbConnect(RSQLite::SQLite(), dbname = paste0(inputDir, "CP_rawOutCProfiler_d20.db"))
conDay40 <-
  dbConnect(RSQLite::SQLite(), dbname = paste0(inputDir, "CP_rawOutCProfiler_d40.db"))
conDay70 <-
  dbConnect(RSQLite::SQLite(), dbname = paste0(inputDir, "CP_rawOutCProfiler_d70.db"))

resPerImage <-
  list(
    d20 = dbReadTable(conDay20, "Per_Image"),
    d40 = dbReadTable(conDay40, "Per_Image"),
    d70 = dbReadTable(conDay70, "Per_Image")
  )

for (el in names(resPerImage))
  resPerImage[[el]]$dataset <- el
resPerImage <- do.call("rbind", resPerImage)

resPerImage$imageKey <-
  paste0(resPerImage$dataset, resPerImage$ImageNumber)
resPerImage$wellKey <-
  paste0(
    resPerImage$dataset,
    substr(resPerImage$Image_Metadata_Well, 1, 1),
    formatNumber2Character(substr(resPerImage$Image_Metadata_Well, 2, 3))
  )

rownames(resPerImage) <- resPerImage$imageKey

resPerImage$Image_Metadata_cellLine <-
  wellAnnot[resPerImage$wellKey, "CellLine"]

#powerLogSlope --> Blur QC of image
powerLogSlopeCols <-
  grep("Image_ImageQuality_PowerLogLogSlope",
       cn(resPerImage),
       value = TRUE)

qplotDensity(apply(resPerImage[, powerLogSlopeCols], 1, min), returnGraph = T) +
  geom_vline(xintercept = -2.5) + xlab("Minumum of PowerLogLogSlope of the 4 channels")
ggsave("out2.1/thresPowerLogLogSlope.pdf",
       width = 5,
       height = 4)

resPerImage$passBlurQC <- apply(resPerImage[, powerLogSlopeCols], 1, min) > -2.5
summary(resPerImage$passBlurQC)

resPerCell <-
  list(
    d20 = dbReadTable(conDay20, "Per_Object"),
    d40 = dbReadTable(conDay40, "Per_Object"),
    d70 = dbReadTable(conDay70, "Per_Object")
  )

for (el in names(resPerCell))
  resPerCell[[el]]$dataset <- el
resPerCell <- do.call("rbind", resPerCell)

#attribute image metadata to cell
resPerCell$imageKey <-
  paste0(resPerCell$dataset, resPerCell$ImageNumber)
for (f in c(grep("Image_Metadata", cn(resPerImage), value = TRUE), "passBlurQC"))
  resPerCell[, f] <- resPerImage[resPerCell$imageKey, f]



resPerCell[, "Image_FileName"] <- resPerImage[resPerCell$imageKey, "Image_FileName_OrigNucleusHoechst"] #get filename in resPerCell

#Only keep cells corresponding to images that passed the blur QC
resPerCell <- resPerCell[resPerCell$passBlurQC, ]
resPerCell <- resPerCell[resPerCell$Image_Metadata_QCFlag_isBlurry != 1 &
                           resPerCell$Image_Metadata_QCFlag_isSaturated != 1 , ]
rm(resPerImage)
gc()

resPerCell$Cytoplasm_AreaShape_Area <-
  resPerCell$Cells_AreaShape_Area - resPerCell$Nuclei_AreaShape_Area

###Split in metadata/feature matrix, rename features

quantitativeCols <- c()
for (col in cn(resPerCell)) {
  if (is.numeric(resPerCell[1:2, col])) {
    if (len(unique(resPerCell[1:10000, col])) > 2) { #not logical
      if (strsplitNth(col, split = "_", n = 1) %in% c("Nuclei", "Cytoplasm", "Cells")) {
        if (strsplitNth(col, split = "_", n = 2) != "Number") {
          quantitativeCols <- c(quantitativeCols, col)
        }
      }
    }
  }
}
otherCol <- setdiff(cn(resPerCell), quantitativeCols)

cellAnnot<-resPerCell[,otherCol]
cellAnnot<-cellAnnot[,c('ImageNumber','ObjectNumber','dataset','imageKey','Image_Metadata_Field',
                        'Image_Metadata_Well',
                        'Image_Metadata_cellLine','Image_FileName')]

colnames(cellAnnot)<-c('imageNumber','objectNumber','diffDay','imageKey','imageField',
                       'well',
                       'cellLine','imageFileName')

cellAnnot$ParentalCellLine<-substr(cellAnnot$cellLine,1,4)

#rename features
featureMat<-resPerCell[,quantitativeCols] |> t()
rm(resPerCell);gc()

featureMat<-featureMat[!grepl("Parent",x= rn(featureMat)),] #unwanted column

hashFamily<- c(
  AreaShape = "Shape",
  Intensity = "Int" ,
  Neighbors = "Neighbors" ,
  Texture = "Txtr",
  Location = "Location"
)

hashType <- c(
  Area = "Area",
  Eccentricity = "Eccentricity",
  EquivalentDiameter = "EquivDiameter",
  EulerNumber = "EulerNumber",
  Extent = "Extent",
  FormFactor = "FormFactor",
  MaximumRadius = "MaxRadius",
  MeanRadius = "MeanRadius",
  MedianRadius = "MedRadius",
  Orientation = "Orientation",
  Solidity = "Solidity",
  IntegratedIntensityEdge = "IntegrIntEdge",
  LowerQuartileIntensity = "LowerQInt",
  MaxIntensityEdge = "MaxIntEdge",
  MaxIntensity = "MaxInt",
  MeanIntensityEdge = "MeanIntEdge",
  MeanIntensity = "MeanInt",
  MedianIntensity = "MedInt",
  MinIntensityEdge = "MinIntEdge",
  MinIntensity = "MinInt",
  NumberOfNeighbors = "NNeighbors",
  PercentTouching = "PctTouching",
  SumAverage = "SumAvg",
  AngularSecondMoment = "AngularScndMoment",
  Correlation = "Corr",
  DifferenceEntropy = "DiffEntro",
  DifferenceVariance = "DiffVar",
  Entropy = "Entropy",
  SumEntropy = "SumEntropy",
  InfoMeas1 = "InfoMeas1",
  InfoMeas2 = "InfoMeas2",
  InverseDifferenceMoment = "InvDiffMnt",
  InverseDifferenceMomnt = "InvDiffMnt",
  MajorAxisLength = "MajorAxisLength",
  MinFeretDiameter = "MinFeretDiameter",
  Perimeter = "Perimeter",
  MADIntensity = "MADInt",
  StdIntensity = "StdInt",
  UpperQuartileIntensity = "UpperQInt",
  Zernike = "Zernike",
  BoundingBoxArea = "BBoxArea",
  BoundingBoxMaximum = "BBoxMax",
  BoundingBoxMinimum = "BBoxMin",
  Center = "Center",
  Compactness = "Compactness",
  BoundingBoxArea = "BBoxArea",
  ConvexArea = "ConvexArea",
  SumVariance = "SumVar",
  Variance = "Var",
  Contrast = "Contrast",
  MaxFeretDiameter = "MaxFeretDiameter",
  IntegratedIntensity = "IntegratedInt",
  MinorAxisLength="MinorAxisLength",
  IntegratedIntensity="IntegratedInt",
  MassDisplacement="MassDisplacement",
  StdIntensityEdge="StdIntEdge",
  CenterMassIntensity="CenterMassInty",
  AngleBetweenNeighbors="AngleBtNghbors",
  FirstClosestDistance="FrstClosestDist",
  FirstClosestObjectNumber="FrstClosestObjNb",
  SecondClosestDistance="ScndClosestDist",
  SecondClosestObjectNumber="ScndClosestObjNb"
  
)

hashChannel<-c(
  OrigCGAN568="CGAN",
  OrigEndoReticulum488="EndoRet",
  OrigMitochondria647="Mito",
  Adjacent="Adjacent",
  OrigEndoRetclm488="EndoRet",
  OrigMitochondr647="Mito",
  OrigEndoReticulm488="EndoRet",
  OrigEndoReticlm488="EndoRet",
  OrigMitochondri647="Mito",
  OrgEndRtclm488="EndoRet",
  OrgMtchndr647="Mito",
  OrigNucleusHoechst="Nuc",
  OrigNucleusHchst="Nuc"
)

renamedFeatureList<-strsplit(rn(featureMat),"_",fixed = T)
renamedFeatureList<-lapply(renamedFeatureList,function(x){
  x[1]<-substr(x[1],1,3)
  if(is.na(hashFamily[x[2]])) print(paste0("family:",x[2]))
  x[2] <- hashFamily[x[2]]
  if(is.na(hashType[x[3]])) print(paste0("type:",x[3]))
  x[3]<- hashType[x[3]]
  if(!is.na(x[4]) & !is.na(hashChannel[x[4]])) x[4]<- hashChannel[x[4]]
  x
})

renamedFeatureList<-sapply(renamedFeatureList, paste0, collapse="_")
names(renamedFeatureList)<-rownames(featureMat)
rownames(featureMat)<-renamedFeatureList

qplotDensity(featureMat["Cyt_Shape_Area",], returnGraph = T) +
           scale_x_log10()
qplotDensity(featureMat["Nuc_Shape_Area",], returnGraph = T) +
  scale_x_log10()

#remove too small / too big cells
featureMat <-featureMat[,
  featureMat["Cyt_Shape_Area",] > 1000 &
  featureMat["Cyt_Shape_Area",] < 10 ^ 5 &
  featureMat["Nuc_Shape_Area",] < 10000
]

featureMat<-na.omit(t(featureMat)) |> t()

#remove too bright / too dark cells
qplotDensity(
  featureMat["Cyt_Int_MeanInt_CGAN",],
  returnGraph = T
) + scale_x_log10()
qplotDensity(
  featureMat["Cyt_Int_MeanInt_Mito",],
  returnGraph = T
) + scale_x_log10()
qplotDensity(
  featureMat["Cyt_Int_MeanInt_EndoRet",],
  returnGraph = T
) + scale_x_log10()
qplotDensity(
  featureMat["Nuc_Int_MeanInt_Nuc",],
  returnGraph = T
) + scale_x_log10()


meanIntColumns<-c("Cyt_Int_MeanInt_CGAN","Cyt_Int_MeanInt_Mito","Cyt_Int_MeanInt_EndoRet","Nuc_Int_MeanInt_Nuc")
cytFluo<-featureMat[c("Cyt_Int_MeanInt_CGAN","Cyt_Int_MeanInt_Mito","Cyt_Int_MeanInt_EndoRet"),] |> apply(2, sum)

QCintDat<-data.frame(Nuc_Int_MeanInt_Nuc=featureMat[c("Nuc_Int_MeanInt_Nuc"),],cytoInt=cytFluo)
QCintDat$dens2d<-oob::pointdensity.nrd(log10(QCintDat)) 

ggplot(data.frame(Nuc_Int_MeanInt_Nuc=featureMat[c("Nuc_Int_MeanInt_Nuc"),],cytoInt=cytFluo,filt=
                    featureMat["Nuc_Int_MeanInt_Nuc",]>0.003 & QCintDat$dens2d>0.0001),
       aes(x=Nuc_Int_MeanInt_Nuc,y=cytoInt,color=filt))+
  scattermore::geom_scattermore()+scale_x_log10()+scale_y_log10()+
  xlab("Nucleus Mean Intensity")+ylab("Sum of Cytoplasm Mean Intensity")


featureMat <-featureMat[, featureMat["Nuc_Int_MeanInt_Nuc",]>0.003 & QCintDat$dens2d>0.0001]
cellAnnot<-cellAnnot[cn(featureMat),]

nFeatPerPlot<-16
pdf("out2.1/allFeatureDistrib.pdf",width = 32,height = 18)
for(i in seq(1,nrow(featureMat),nFeatPerPlot)){
  features<-rn(featureMat)[i:min(i+nFeatPerPlot,nrow(featureMat))]
  plots<-vector(mode ="list", nFeatPerPlot*2)
  j=1;for(f in features){
    plots[[j]] <-  qplotDensity(
      featureMat[f,],
      returnGraph = T
    ) + ggtitle(f)
    transfeat<-featureMat[f,]-min(featureMat[f,],na.rm=T)
    transfeat<-log2(transfeat/max(transfeat,na.rm=T)*1000+1)
    plots[[j+1]] <-  qplotDensity(
      transfeat,
      returnGraph = T
    ) + ggtitle(paste0(f," trans"))
    j=j+2
  }
  multiplot(plotlist = plots,cols = 8)
}
dev.off()


## Regularization of CP features
#extract spatial infos
regularizationDat<-fastRead(paste0(inputDir,"CP_FeatureRegularization.tsv"))
regu<-regularizationDat$transformation;names(regu)<-rn(regularizationDat)
rm(regularizationDat)

spatialAnnotFeatures<-names(regu)[regu=="move2spatialMat"]
spatialInfoAnnot<-t(featureMat[spatialAnnotFeatures,])
featureMat<-featureMat[!rn(featureMat) %in% spatialAnnotFeatures,]
saveRDS(spatialInfoAnnot,"out2.1/spatialInfoAnnot.rds")
rm(spatialInfoAnnot)

saveRDS(featureMat,"out2.1/rawFeatureMat.rds")

#Put everything on the range 0-1000
featureMat<-apply(featureMat,1,function(x){
  x<-x-min(x) #min = 0
  x<-x/max(x) * 1000
  x
}) |> t()

#Per feature type transformation

for(f in names(regu)[regu=="log2(x+1)"]){
  featureMat[f,]<-log2(featureMat[f,]+1)
  featureMat[f,]<-featureMat[f,]/max(featureMat[f,])*1000 #put again on range 0-1000
}

for(f in names(regu)[regu=="log2(x+0.01)"]){
  x<-featureMat[f,]
  x<-log2(x+0.01)
  x<-x - min(x)
  featureMat[f,]<-x/max(x)*1000
  #qplotDensity(featureMat[f,])
}

for(f in names(regu)[regu=="log2(max(x)-x+0.1)"]){
  x<- featureMat[f,]
  x<- log2(max(x)-x+0.1) 
  x<- max(x)-x
  x<- x - min(x)
  featureMat[f,]<-x/max(x)*1000
  #qplotDensity(featureMat[f,])
}

for(f in names(regu)[regu=="log2(max(x)-x+0.5)"]){
  x<- featureMat[f,]
  x<- log2(max(x)-x+0.5) 
  x<- max(x)-x
  x<- x - min(x)
  featureMat[f,]<-x/max(x)*1000
  #plot(qplotDensity(x/max(x)*1000,returnGraph = T)+ggtitle(f))
}

for(f in names(regu)[regu=="log2(max(x)-x+0.01)"]){
  x<- featureMat[f,]
  x<- log2(max(x)-x+0.01) 
  x<- max(x)-x
  x<- x - min(x)
  featureMat[f,]<-x/max(x)*1000
  #plot(qplotDensity(x/max(x)*1000,returnGraph = T)+ggtitle(f))
}

for(f in names(regu)[regu=="2{ log2(max(x)-x+1)}"]){
  x<- featureMat[f,]
  x<- log2(max(x)-x+1)
  x<- log2(max(x)-x+.5)
  featureMat[f,]<-x
  #qplotDensity(x)
}

nFeatPerPlot<-20
pdf("out2.1/allFeatureDistribPostRegularization.pdf",width = 32,height = 18)
for(i in seq(1,nrow(featureMat),nFeatPerPlot)){
  features<-rn(featureMat)[i:min(i+nFeatPerPlot,nrow(featureMat))]
  plots<-vector(mode ="list", len(features));names(plots)<-features
  for(f in features){
    plots[[f]] <-  qplotDensity(
      featureMat[f,],
      returnGraph = T
    ) + ggtitle(f)
  }
  multiplot(plotlist = plots,cols = 5)
}
dev.off()

saveRDS(featureMat,"out2.1/regFeatureMat.rds")
fastWrite(cellAnnot,"out2.1/cellAnnot.tsv")

## Dataset Subsampling, trying to have equal number of cells from each experimental population
cellAnnot$expPopWithReplicate<-paste0(cellAnnot$diffDay, "_", cellAnnot$cellLine)
cellAnnot$expPop<-substr(cellAnnot$expPopWithReplicate,1,8)

tablePops <- table(as.factor(cellAnnot$expPopWithReplicate))
tablePops<-tablePops[names(tablePops)!="d40_82.6_3.1"] #only 3 cells

samplePerExpPop<- lapply(names(tablePops), function(lvl) {
  subsampling <- rn(cellAnnot)[which(cellAnnot$expPopWithReplicate == lvl)]
  sample(subsampling,size = min(length(subsampling),2000))
});names(samplePerExpPop)<-names(tablePops)
subsampling<-unlist(samplePerExpPop)
subsampling<-sample(subsampling,len(subsampling)) #shuffle the result

featureMat<-featureMat[,subsampling]
cellAnnot<-cellAnnot[subsampling,]


## FEATURE SELECTION

# keep feature among those that have same values 
corFeatures <- cor(t(featureMat))
vectCorFeature <- corFeatures[upper.tri(corFeatures)]


graphOfFeature <-
  graph_from_adjacency_matrix(corFeatures > .99, mode = "upper", diag = F)


c1 = cluster_fast_greedy(graphOfFeature)

pdf("out2.1/moduleOfFeatures.pdf",width = 7,height = 7)
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

#select manually in each module the feature to keep
featurePerModule[[1]]<-"Cyt_Int_StdInt_EndoRet"
featurePerModule[[2]]<-"Cel_Int_StdInt_CGAN"
featurePerModule[[3]]<-"Nuc_Txtr_Var_Nuc_3_00_256"
featurePerModule[[4]]<-"Cyt_Int_MeanInt_EndoRet"
featurePerModule[[5]]<-"Cyt_Int_StdInt_Mito"
featurePerModule[[6]]<-"Cyt_Int_MeanInt_CGAN"
featurePerModule[[7]]<-"Cyt_Int_MeanInt_Mito"
featurePerModule[[8]]<-"Cyt_Txtr_AngularScndMoment_Mito_3_00_256"
featurePerModule[[9]]<-"Cyt_Txtr_Entropy_Mito_3_00_256"
featurePerModule[[10]]<-"Cyt_Txtr_AngularScndMoment_EndoRet_3_00_256"
featurePerModule[[11]]<-"Cyt_Txtr_AngularScndMoment_CGAN_3_00_256"
featurePerModule[[12]]<-"Cyt_Txtr_SumAvg_Mito_3_00_256"
featurePerModule[[13]]<-"Cyt_Txtr_Entropy_EndoRet_3_00_256"
featurePerModule[[14]]<-"Nuc_Txtr_Entropy_Nuc_3_00_256"
featurePerModule[[15]]<-"Cyt_Txtr_Entropy_CGAN_3_00_256"
featurePerModule[[16]]<-"Cyt_Int_MinInt_Mito"
featurePerModule[[17]]<-"Nuc_Txtr_AngularScndMoment_Nuc_3_00_256"
featurePerModule[[18]]<-"Cyt_Txtr_DiffVar_EndoRet_3_00_256"
featurePerModule[[19]]<-"Cyt_Int_MinInt_EndoRet"
featurePerModule[[20]]<-"Nuc_Txtr_SumAvg_Nuc_3_00_256"
featurePerModule[[21]]<-"Cyt_Txtr_DiffVar_Mito_3_00_256"
featurePerModule[[22]]<-"Cyt_Int_MinInt_CGAN"
featurePerModule[[23]]<-"Cyt_Txtr_DiffVar_CGAN_3_00_256"
featurePerModule[[24]]<-"Cyt_Int_MaxInt_EndoRet"
featurePerModule[[25]]<-"Cyt_Int_MaxInt_CGAN"
featurePerModule[[26]]<-"Cyt_Int_IntegrIntEdge_CGAN"
featurePerModule[[27]]<-"Nuc_Shape_MinorAxisLength"
featurePerModule[[28]]<-"Nuc_Shape_MeanRadius"
featurePerModule[[29]]<-"Nuc_Shape_Area"
featurePerModule[[30]]<-"Cyt_Int_MeanIntEdge_Mito"
featurePerModule[[31]]<-"Cyt_Int_MeanIntEdge_EndoRet"
featurePerModule[[32]]<-"Cyt_Txtr_Contrast_EndoRet_3_00_256"
featurePerModule[[33]]<-"Cyt_Int_MeanIntEdge_CGAN"
featurePerModule[[34]]<-"Cyt_Int_MaxInt_Mito"
featurePerModule[[35]]<-"Cyt_Int_IntegrIntEdge_EndoRet"
featurePerModule[[36]]<-"Cel_Shape_Area"

selectedFeatures <- unlist(featurePerModule)

#Heatmap of all feature on a random subset of cells
pdf("out2.1/completeHeatmap.pdf",width = 12,height = 30)
heatmap.DM(featureMat[selectedFeatures,sample.int(ncol(featureMat),  1000)],
          name = "regularized\nfeature\nvalues")
dev.off()

#zoom on Zernike feature.. they are bad
zernikeFeat<-grep("Zernike",selectedFeatures,value = T)

pdf("out2.1/HeatmapZernike.pdf",width = 20,height = 30)
heatmap.DM(featureMat[zernikeFeat,sample.int(ncol(featureMat),  1000)],
           name = "regularized\nfeature\nvalues")
dev.off()

#UMAP, unsupervised clustering, export to web app for removing bad clusters by inspecting images
umap<-UMAP(featureMat[setdiff(selectedFeatures,zernikeFeat),],ret_nn = T)
cellAnnot$leidenClust<-leidenFromUMAP(umap,n_neighbors = 20, resolution_parameter = 1)

proj2d(umap, useScatterMore=T,colorBy = cellAnnot$diffDay)
proj2d(umap, useScatterMore=T,colorBy = cellAnnot$leidenClust,plotFactorsCentroids = T)

selectedFeature<-setdiff(selectedFeatures,c(zernikeFeat,"Cel_Neighbors_AngleBtNghbors_Adjacent"))
selectedSamples<-rn(cellAnnot)[!cellAnnot$leidenClust %in% paste0("k",10:14)] 
# Originally, bad clusters are k10,11,12,13,14 based, of the pictures of cells
# belonging to these clusters

cellAnnot<-cellAnnot[selectedSamples,]
featureMat<-featureMat[selectedFeature,selectedSamples]

fastWrite(cellAnnot,"out2.1/cellAnnotFiltered.tsv")
fastWrite(featureMat,"out2.1/featureMatFiltered.tsv")
