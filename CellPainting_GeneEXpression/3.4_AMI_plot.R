setwd(".")
setwd("C:/Users/dmeister/OneDriveUniv/KilpinenPostdoc/tryCPcode/Code/")

library(oob)
library(FNN)
library(ggplot2)
library(glmnet)
library(stringr)
library(ComplexHeatmap)

gexData<-t(readRDS("../createPredMatV2/inputs/gexMat.rds"))
#randomizeCell
randomDraw<-sample.int(nrow(gexData), nrow(gexData) )
gexData<-gexData[randomDraw,]

gexAnnot<-fastRead("inputs/cellAnnotGex.tsv")[randomDraw,]

#only keep cell line in common?

cpData<-fastRead("../createPredMatV2/inputs/CPfeatureMat.tsv",as.matrix = TRUE) |> t()

cpAnnot<-formatAnnotFromMeta(fastRead("inputs/cellAnnotCP.tsv"),fastRead("inputs/metaAnnotCP.tsv"))
gexAnnot<-formatAnnotFromMeta(annotDataFrame = gexAnnot,metaAnnot = fastRead("inputs/metaAnnotGex.tsv"))

colorScales<-attr(cpAnnot,"colorScales")

colorScaleCTfile<-fastRead("./inputs/colorCellType.tsv",row.names = NULL,header = TRUE)
colorScales$newAnnotCellType<-colorScaleCTfile$Color
colorScales$MapRefPolioudakisRharmonyPredCtype<-colorScaleCTfile$Color[1:15]
names(colorScales$newAnnotCellType)<-colorScaleCTfile$Lvl
names(colorScales$MapRefPolioudakisRharmonyPredCtype)<-colorScaleCTfile$Lvl[1:15]
gexAnnot$newAnnotCellType<-factor(gexAnnot$newAnnotCellType, levels = colorScaleCTfile$Lvl)
gexAnnot$MapRefPolioudakisRharmonyPredCtype<-factor(gexAnnot$MapRefPolioudakisRharmonyPredCtype, levels = colorScaleCTfile$Lvl[1:15])

gexFromCpTheoMat<-readRDS("../createPredMatV2/results/gexFromCpTheoMat.rds")
cpFromGexTheoMat<-readRDS("../createPredMatV2/results/cpFromGexTheoMat.rds")[randomDraw,]

cpProj<-fastRead("inputs/umapCP.tsv",as.matrix = T)
gexProj<-fastRead("inputs/umapGex.tsv",as.matrix = T)[randomDraw,]

CPknn<-get.knnx(cpFromGexTheoMat,cpData, k=1)

predictedAnnotCP2GEX<-data.frame(row.names = rn(cpAnnot))
for(annot in cn(gexAnnot)){
    predictedAnnotCP2GEX[,annot]<-gexAnnot[CPknn$nn.index,annot]
}
fastWrite(predictedAnnotCP2GEX,"outText/predictedAnnotCP2GEX.tsv")

pdf("outPlot/projGexAnnotOnCPproj.pdf", width=14,height=12)
proj2d(cpProj,  colorBy = predictedAnnotCP2GEX$MapRefPolioudakisRharmonyPredCtype, colorScale = colorScales$MapRefPolioudakisRharmonyPredCtype, 
 useScatterMore = T,legendTitle = "MapRefPolioudakisRharmonyPredCtype")
proj2d(cpProj,  colorBy = predictedAnnotCP2GEX$newAnnotCellType, colorScale = colorScales$newAnnotCellType,   
useScatterMore = T, legendTitle = "newAnnotCellType")
proj2d(cpProj,  colorBy = predictedAnnotCP2GEX$diffDay |> as.factor() ,   
    useScatterMore = T, legendTitle = "differentiationDay\n(NearestNeighbour)",colorScale = colorScales$diffDay)
proj2d(cpProj,  colorBy = cpAnnot$diffDay |> as.factor() ,   useScatterMore = T, 
    legendTitle = "differentiationDay(CP)",,colorScale = colorScales$diffDay)
dev.off()

progCellLines<-c("PgS","PgG2M","IP")
cellOrder<-c(rn(predictedAnnotCP2GEX)[predictedAnnotCP2GEX$MapRefPolioudakisRharmonyPredCtype%in%progCellLines],
rn(predictedAnnotCP2GEX)[!predictedAnnotCP2GEX$MapRefPolioudakisRharmonyPredCtype%in%progCellLines])

pdf("outPlot/projGexAnnotOnCPprojCheat.pdf", width=14,height=12)
proj2d(cpProj[cellOrder,],  colorBy = predictedAnnotCP2GEX[cellOrder,"MapRefPolioudakisRharmonyPredCtype"], 
    colorScale = colorScales$MapRefPolioudakisRharmonyPredCtype, 
    useScatterMore = T,legendTitle = "MapRefPolioudakisRharmonyPredCtype")
dev.off()

#facetWrap by newAnnotCellType
pdf("outPlot/projGexAnnotOnCPprojFacetCellType.pdf", width=12,height=12)
print(
    proj2d(cpProj,  colorBy = predictedAnnotCP2GEX$MapRefPolioudakisRharmonyPredCtype,   useScatterMore = T, legendTitle = "MapRefPolioudakisRharmonyPredCtype",returnGraph = TRUE,
    colorScale = colorScales$MapRefPolioudakisRharmonyPredCtype)+
        ggplot2::facet_wrap(~predictedAnnotCP2GEX$newAnnotCellType,ncol = 3)
)
print(
    proj2d(cpProj,  colorBy = predictedAnnotCP2GEX$newAnnotCellType,   useScatterMore = T, legendTitle = "newAnnotCellType",returnGraph = TRUE,
    colorScale = colorScales$newAnnotCellType)+
        ggplot2::facet_wrap(~predictedAnnotCP2GEX$newAnnotCellType,ncol = 3)
)
print(
    proj2d(cpProj,  colorBy = predictedAnnotCP2GEX$differentiationDay |> as.factor() ,   useScatterMore = T, legendTitle = "differentiationDay\n(NearestNeighbour)",returnGraph = TRUE)+
        ggplot2::facet_wrap(~predictedAnnotCP2GEX$newAnnotCellType,ncol = 3)
)
print(
    proj2d(cpProj,  colorBy = cpAnnot$diffDay |> as.factor(),   useScatterMore = T, legendTitle = "differentiationDay(CP)",returnGraph = TRUE)+
        ggplot2::facet_wrap(~predictedAnnotCP2GEX$newAnnotCellType,ncol = 3)
)
dev.off()



#facetWrap by diff day
pdf("outPlot/projGexAnnotOnCPprojFacetDiffDayCP.pdf", width=12,height=12)
print(
    proj2d(cpProj,  colorBy = predictedAnnotCP2GEX$MapRefPolioudakisRharmonyPredCtype,   useScatterMore = T, legendTitle = "MapRefPolioudakisRharmonyPredCtype",returnGraph = TRUE,
    colorScale = colorScales$MapRefPolioudakisRharmonyPredCtype)+
        ggplot2::facet_wrap(~cpAnnot$diffDay |> as.factor(),ncol = 3)
)
print(
    proj2d(cpProj,  colorBy = predictedAnnotCP2GEX$newAnnotCellType,   useScatterMore = T, legendTitle = "newAnnotCellType",returnGraph = TRUE,
    colorScale = colorScales$newAnnotCellType)+
        ggplot2::facet_wrap(~cpAnnot$diffDay |> as.factor(),ncol = 3)
)
print(
    proj2d(cpProj,  colorBy = predictedAnnotCP2GEX$differentiationDay |> as.factor() ,   useScatterMore = T, legendTitle = "differentiationDay\n(NearestNeighbour)",returnGraph = TRUE)+
        ggplot2::facet_wrap(~cpAnnot$diffDay |> as.factor(),ncol = 3)
)
print(
    proj2d(cpProj,  colorBy = cpAnnot$diffDay |> as.factor(),   useScatterMore = T, legendTitle = "differentiationDay(CP)",returnGraph = TRUE)+
        ggplot2::facet_wrap(~cpAnnot$diffDay |> as.factor(),ncol = 3)
)
dev.off()

#what is the good prediction rate?

isGoodPred<-predictedAnnotCP2GEX$diffDay == cpAnnot$diffDay
table(isGoodPred)/nrow(cpAnnot)
#random prediction rate
randomPred<-predictedAnnotCP2GEX$diffDay == sample(c("d20","d40","d70"),nrow(predictedAnnotCP2GEX),replace = T)
table(randomPred)/nrow(cpAnnot)

isPartition1<-sample(c(TRUE,FALSE),nrow(gexData),replace = T)

gExp1<-gexData[isPartition1,]
gExp2<-gexData[!isPartition1,]


gexPartitionknn<-get.knnx(gExp1,gExp2, k=1)

predDiffDayGexPartition<-gexAnnot[isPartition1,"diffDay"][gexPartitionknn$nn.index]
gexxRandomCorrespondance<-predDiffDayGexPartition==gexAnnot[!isPartition1,"diffDay"]
table(gexxRandomCorrespondance)/len(gexxRandomCorrespondance)

rm(gExp1,gExp2);gc()

#examine what variable predict wrong prediction outcome for diff day
predGoodDf<-data.frame(
    isGoodPred = ifelse(isGoodPred, 1, 0),
    cpAnnot[, c("leidenClust","parentalCellLine")],
    cpData,knnDist=CPknn$nn.dist[,1]
)

options(expressions = 5e5)
model_matrix <- model.matrix(isGoodPred ~ ., data = predGoodDf)

fit <- glmnet(model_matrix, predGoodDf$isGoodPred, family = "binomial",alpha = .5,
    lambda=cv.glmnet(model_matrix, predGoodDf$isGoodPred, family = "binomial")$lambda.1se)

write.table(as.matrix(coef(fit)),file = "./outText/contribGoodPred.tsv",
    row.names = TRUE, quote = FALSE,sep="\t")

## genes

fit <- glmnet(gexFromCpTheoMat, ifelse(isGoodPred, 1, 0), family = "binomial",alpha = .5,
    lambda=cv.glmnet(gexFromCpTheoMat, ifelse(isGoodPred, 1, 0), family = "binomial")$lambda.1se)

coefs<-as.matrix(coef(fit))
coefs<-coefs[coefs!=0,]

gene2test<-names(coefs)[coefs>0]

write(gene2test,"./outText/gene2test.txt")

##CPdata

fit <- glmnet(cpData, ifelse(isGoodPred, 1, 0), family = "binomial",alpha = .5,
    lambda=cv.glmnet(cpData, ifelse(isGoodPred, 1, 0), family = "binomial")$lambda.1se)

coefs<-as.matrix(coef(fit))
coefs<-coefs[coefs!=0,]

#gex 2 cp analysis
GEXknn<-get.knnx(gexFromCpTheoMat,gexData, k=1)

predictedAnnotGEX2CP<-data.frame(row.names = rn(gexAnnot))
for(annot in cn(cpAnnot)) predictedAnnotGEX2CP[,annot]<-cpAnnot[GEXknn$nn.index,annot]

fastWrite(predictedAnnotGEX2CP,"outText/predictedAnnotGEX2CP.tsv")

pdf("outPlot/projCPAnnotOnGex.pdf", width=14,height=12)
proj2d(gexProj,  colorBy = predictedAnnotGEX2CP$leidenClust,   useScatterMore = T,legendTitle = "CP leiden clust")
proj2d(gexProj,  colorBy = predictedAnnotGEX2CP$diffDay ,   useScatterMore = T, colorScale = colorScales$diffDay,
    legendTitle = "differentiationDay\n(NearestNeighbour)")
proj2d(gexProj,  colorBy = gexAnnot$diffDay |> as.factor() , colorScale = colorScales$diffDay, 
    useScatterMore = T, legendTitle = "differentiationDay(Gex)")
dev.off()

##
pdf("./outPlot/barplotFreqCell.pdf",  width=10,height=9)
ggData<-data.frame(predictedAnnotCP2GEX,cpAnnot)
ggData$expPop<-factor(str_replace(ggData$expPop.1,"_","L"), levels = names(colorScales$expPop))

print(ggplot(ggData,mapping=aes(x=leidenClust,fill=newAnnotCellType))+
	geom_bar(position="fill", color="black")+theme_bw()+
	scale_fill_manual(values = colorScales$newAnnotCellType))

print(ggplot(ggData,mapping=aes(x=leidenClust,fill=MapRefPolioudakisRharmonyPredCtype))+
	geom_bar(position="fill", color="black")+theme_bw()+
	scale_fill_manual(values = colorScales$MapRefPolioudakisRharmonyPredCtype))

print(ggplot(ggData,mapping=aes(x=newAnnotCellType,fill=leidenClust))+
	geom_bar(position="fill", color="black")+theme_bw()+
	scale_fill_manual(values = oobColors(11)))

print(ggplot(ggData,mapping=aes(x=MapRefPolioudakisRharmonyPredCtype,fill=leidenClust))+
	geom_bar(position="fill", color="black")+theme_bw()+
	scale_fill_manual(values = oobColors(11)))

print(ggplot(ggData,mapping=aes(x=MapRefPolioudakisRharmonyPredCtype,fill=expPop))+
	geom_bar(position="fill", color="black")+theme_bw()+
	scale_fill_manual(values = colorScales$expPop))

dev.off()

pdf("./outPlot/barplotFreqCellPerDonorPerDiffDay.pdf",  width=10,height=9)
print(
    ggplot(ggData,mapping=aes(x=parentalCellLine,fill=newAnnotCellType))+
    geom_bar(position="fill", color="black")+theme_bw()+
    scale_fill_manual(values = colorScales$newAnnotCellType)+
    facet_wrap(~diffDay, nrow=3)
)
print(
    ggplot(ggData,mapping=aes(x=parentalCellLine,fill=MapRefPolioudakisRharmonyPredCtype))+
    geom_bar(position="fill", color="black")+theme_bw()+
    scale_fill_manual(values = colorScales$MapRefPolioudakisRharmonyPredCtype)+
    facet_wrap(~diffDay, nrow=3)
)
dev.off()

#compute closest cell matrix
cpFromGexClosestNeighbor<-cpData[GEXknn$nn.index,]
fastWrite(cpFromGexClosestNeighbor,"./outText/cpFromGexClosestNeighbor.tsv")

fastWrite(gexData,"./outText/gexData.tsv")
fastWrite(cpData,"./outText/cpData.tsv")
fastWrite(gexProj,"./outText/gexProj.tsv")
fastWrite(cpProj,"./outText/cpProj.tsv")
fastWrite(gexFromCpTheoMat,"./outText/gexFromCpTheoMat.tsv")
fastWrite(cpFromGexTheoMat,"./outText/cpFromGexTheoMat.tsv")
fastWrite(gexAnnot,"./outText/gexAnnot.tsv")
fastWrite(cpAnnot,"./outText/cpAnnot.tsv")

#### correlation heatmap
#meta cell method


pdf("outPlot/projMitoFeatures.pdf",  width=10,height=9)
proj2d(gexProj,colorBy = cpFromGexTheoMat[,"Cyt_Int_MeanInt_Mito"],
    useScatterMore=TRUE, main = "Cyt_Int_MeanInt_Mito",
    colorScale = c("darkblue","blue","cyan","lightgreen","orange","red","darkred"))
proj2d(gexProj,colorBy = cpFromGexTheoMat[,"Cyt_Txtr_AngularScndMoment_Mito_3_00_256"],
    useScatterMore=TRUE, main = "Cyt_Txtr_AngularScndMoment_Mito_3_00_256",
    colorScale = c("darkblue","blue","cyan","lightgreen","orange","red","darkred"))
dev.off()
#compute more resolved UMAP

proj2d(cpProj,cpAnnot$leidenClust,useScatterMore = T)
proj2d(gexProj,gexAnnot$newAnnotCellType,useScatterMore = T)
save.image("./outRobj/workspace.RData")

###
markerPerCPcluster<-getMarkers(gexFromCpTheoMat |> t(),cpAnnot$leidenClust)
markerPerCPcluster<-oob::extractFeatureMarkerData(markerPerCPcluster)
markerPerCPcluster<-as.matrix(markerPerCPcluster)

oob::qplotDensity(markerPerCPcluster |> as.vector())

dbTerms<-oob::getDBterms(rn(markerPerCPcluster),database = c("kegg","goCC"))

x<-isOfInterest
corrIdGenes = NULL; database = c("kegg", "reactom", 
    "goBP", "goCC", "goMF"); minSize = 2; maxSize = 500; returnGenes = FALSE;
    keggDisease = FALSE; species = "Human"; customAnnot = NULL; 
    db_terms = dbTerms; speciesData = NULL

enrich.ora<-function (x, corrIdGenes = NULL, database = c("kegg", "reactom", 
    "goBP", "goCC", "goMF"), minSize = 2, maxSize = 500, returnGenes = FALSE, 
    keggDisease = FALSE, species = "Human", customAnnot = NULL, 
    db_terms = NULL, speciesData = NULL) 
{
    validDBs <- c("kegg", "reactom", "goBP", "goCC", "goMF", 
        "custom")
    if (sum(database %in% validDBs) == 0) 
        stop(paste0("Error, valid values for database are: ", 
            paste0(validDBs, collapse = ", ")))
    if (is.null(customAnnot) & "custom" %in% database) 
        stop("You must give a value a list in customAnnot if database=custom")
    if (is.data.frame(x) | is.matrix(x)) {
        tempx <- x
        x <- tempx[, 1]
        names(x) <- rownames(tempx)
    }
    if (class(x) != "logical") 
        stop("Values must be logical (TRUE or FALSE)")
    if (class(names(x)) != "character") 
        stop("Values must be named with genes symbol")
    if (is.null(db_terms)) 
        db_terms <- getDBterms(geneSym = names(x), corrIdGenes = corrIdGenes, 
            database = database, customAnnot = customAnnot, keggDisease = keggDisease, 
            species = species)
    nInterest <- length(which(x))
    nuniverse <- length(x)
    results <- list()
    for (db in names(db_terms)) {
        len_term <- sapply(db_terms[[db]], length)
        db_terms[[db]] <- db_terms[[db]][len_term >= minSize & 
            len_term <= maxSize]
        nGeneByterm <- sapply(db_terms[[db]], length)
        nGeneOfInterestByterm <- sapply(db_terms[[db]], function(term) {
            return(length(which(x[term])))
        })
        results[[db]] <- data.frame(row.names = names(db_terms[[db]]))
        results[[db]]$term <- names(db_terms[[db]])
        results[[db]]$database <- db
        parameterList4Enrich <- vector(mode = "list", length = length(results[[db]]$term))
        for (i in seq_along(results[[db]]$term)) {
            parameterList4Enrich[[i]] <- list(intersectionSize = nGeneOfInterestByterm[i], 
                setSizes = c(nInterest, nGeneByterm[i]), universeSize = nuniverse)
        }
        resEnrich <- lapply(parameterList4Enrich, function(params) do.call("enrichSetIntersection", 
            params))
        results[[db]]$nGene <- nGeneByterm
        results[[db]]$obsOverlap <- nGeneOfInterestByterm
        results[[db]]$expectOverlap <- sapply(resEnrich, function(x) x$expected)
        results[[db]]$log2OE <- sapply(resEnrich, function(x) x$log2OE)
        results[[db]]$pval <- sapply(resEnrich, function(x) x$pval)
        results[[db]]$padj <- 0
        if (returnGenes) {
            results[[db]]$genes <- db_terms[[db]]
        }
    }
    results <- do.call("rbind", results)
    results$padj <- p.adjust(results$pval, method = "BH")
    return(results)
}


enrichPerCPcluster<-list()
for(cluster in cn(markerPerCPcluster)){
    isOfInterest<-markerPerCPcluster[,cluster]>3
    names(isOfInterest)<-rownames(markerPerCPcluster)
    enrichPerCPcluster[[cluster]]<-enrich.ora(isOfInterest,db_terms = dbTerms)
}

for(cluster in cn(markerPerCPcluster)) fastWrite(enrichPerCPcluster[[cluster]],paste0("./outText/enrichORA",cluster,".tsv"))

###

dir.create("./outPlot/projFeaturesOnUMAP",showWarnings = FALSE)

for(feature in cn(cpFromGexTheoMat)){
    pdf(paste0("outPlot/projFeaturesOnUMAP/",feature,".pdf"),  width=20,height=20)
    print(proj2d(gexProj,colorBy = cpFromGexTheoMat[,feature],
        useScatterMore=TRUE, main = feature,returnGraph = TRUE,
        colorScale = c("darkblue","blue","cyan","lightgreen","orange","red","darkred"))+
        ggplot2::facet_wrap(~gexAnnot$newAnnotCellType |> as.factor(),ncol = 4))
    dev.off()
}

plotExpr(cpFromGexTheoMat[,feature]|>t(),group=cellAnnot)

###BEST cp marker of gene expression
gexAnnot<-fastRead("../data/gexAnnot.tsv")
featureCPmodule<-fastRead("inputs/featureAnnotCP.tsv")

featMarkerOfCellType<-getMarkers(cpFromGexTheoMat |> t(),gexAnnot$finalAnnotation)
featMarkerOfCP<-getMarkers(cpData |> t(),cpAnnot$cellClusters)
fastWrite(featMarkerOfCellType,"./outText/featMarkerOfCellType.tsv")

featMarkerOfCellTypeExtracted<-extractFeatureMarkerData(featMarkerOfCellType)
featMarkerOfCPExtracted<-extractFeatureMarkerData(featMarkerOfCP)

featMarkerOfCellTypeExtracted<-featMarkerOfCellTypeExtracted[,cn(featMarkerOfCellTypeExtracted)!="NA."]

pdf("outPlot/featureMarkerScoreCPclust_cellTypes.pdf",width=14, height=12)
heatmap.DM(featMarkerOfCPExtracted, preSet = NULL,   center=FALSE,returnHeatmap = TRUE,
    row_split = featureCPmodule, row_title_rot=0)+
    heatmap.DM(featMarkerOfCellTypeExtracted, preSet = NULL,   center=FALSE,returnHeatmap = TRUE)
dev.off()

pdf("outPlot/corCTypeCPclusterOnspecifictyScore.pdf",  width=7.5,height=9)
heatmap.DM(cor(featMarkerOfCellTypeExtracted,featMarkerOfCPExtracted),preSet = "cor")
dev.off()

drawGex<-sample(rn(gexAnnot), size=nrow(gexAnnot))

pdf("outPlot/projForPoster.pdf",  width=10,height=9)
proj2d(gexProj[drawGex,],colorBy=gexAnnot[drawGex,]$diffDay,useScatterMore = T,colorScale=colorScales$diffDay)
proj2d(gexProj[drawGex,],,colorBy=gexAnnot[drawGex,]$newAnnotCellType,useScatterMore = T,colorScale=colorScales$newAnnotCellType)
dev.off()
