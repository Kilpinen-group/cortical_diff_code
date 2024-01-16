setwd(".")

library(oobLegacy)
library(ComplexHeatmap)

linearModels<-fastRead("inputs/GeXmodelLinear.tsv")
featurePerModule<-fastRead("inputs/featureAnnot.tsv")
# pca2d(PCA(linearModels),plotText = T)

intercepts<-linearModels[1,rn(featurePerModule)]
linearModels<-linearModels[-1,rn(featurePerModule)]

qplotBarplot(colSums(abs(linearModels)))

#fig 1A

feature2Show<-c("Cyt_Int_MeanInt_Mito", "Cyt_Int_MeanInt_CGAN", "Cyt_Int_MeanInt_EndoRet","Nuc_Int_MeanInt_Nuc", "Cel_Shape_Area", "Nuc_Shape_Area", "Cel_Neighbors_PctTouching_Adjacent")

coefMatrixOfFeature2Show<-linearModels[,feature2Show]
maxCoefPerGene<-apply(abs(coefMatrixOfFeature2Show),1,max)

qplotDensity(maxCoefPerGene)
sum(maxCoefPerGene > 2.5)
#bestGenesOfSelFeature <- sort(maxCoefPerGene,decreasing = T)[1:20] |> names()

bestGenesOfSelFeature<- apply(coefMatrixOfFeature2Show,2,function(x) names(x)[order(abs(x),decreasing = T)] )[1:2,] |> as.vector() |> unique()

pdf("outPlot/heatmapFig1A.pdf",width = 6,height = 8.6)
heatmap.DM(coefMatrixOfFeature2Show[bestGenesOfSelFeature,],preSet = "vanilla")
dev.off()

modules<-unique(featurePerModule$module)
medCoefPerModel<-aggregMatPerVector(linearModels,featurePerModule$module,FUN = median)

#coefOfMod<-list()
#for(mod in modules){
#    featureOfMod<-rn(featurePerModule)[featurePerModule$module == mod]
#    linearModels[,featureOfMod]
#}

head(medCoefPerModel)


dbs<-oob::getDBterms(geneSym = rn(medCoefPerModel),database = c("kegg","goBP","goCC","goMF"))

resEnrichPerModule<-list()
for(mod in modules){
    resEnrichPerModule[[mod]]<-enrich.fcs(medCoefPerModel[,mod],db_terms = dbs)
}
for(mod in modules) fastWrite(resEnrichPerModule[[mod]],paste0("outText/",mod,".tsv"))

for(mod in modules){
    pdf(paste0("outPlot/volcano_",mod,".pdf"),width = 9,height = 8)
    volcanoPlot(resEnrichPerModule[[mod]],effectSizeCol = "NES", adjPvalCol = "padj",labelCol = "pathway")
    dev.off()
}

cn(linearModels)


