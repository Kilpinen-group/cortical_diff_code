setwd("C:/Users/dmeister/OneDriveUniv/KilpinenPostdoc/tryCPcode/Code/")
setwd(".")

inputDir <- "zenodoArchive/"
dir.create("out2.3", showWarnings = FALSE)

library(oob)
library(harmony)
library(plotly)

colorScaleCluster <-c(
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

print_app <- function(widget) {
    
    # Generate random file name
    temp <- paste(tempfile('plotly'), 'html', sep = '.')
    
    # Save. Note, leaving selfcontained=TRUE created files that froze my browser
    htmlwidgets::saveWidget(widget, temp, selfcontained = FALSE)
    
    # Launch with desired application
    system(sprintf("chromium-browser -app=file://%s", temp))
    
    # Return file name if it's needed for any other purpose
    temp
}

diffAbundance <- function(clusters,samples, sampleAnnot){
    require(edgeR)
    abundance <- table(clusters,samples)
    print(abundance)
    sampleCol <- colnames(sampleAnnot)[1]
    sampleAnnot<- sampleAnnot[colnames(abundance),,drop=F]
    y.ab <- DGEList(abundance, samples=sampleAnnot)
    keep <- filterByExpr(y.ab, group=y.ab$samples[,sampleCol])
    y.ab <- y.ab[keep,T]
    print(summary(keep))
    design <- model.matrix(as.formula(paste0("~factor(",sampleCol,")")), y.ab$samples)
    y.ab <- estimateDisp(y.ab, design, trend="none")
    fit.ab <- glmQLFit(y.ab, design, robust=TRUE, abundance.trend=FALSE)
    res <- glmQLFTest(fit.ab, coef=ncol(design))
    res$abundance <- abundance
    print(summary(decideTests(res)))
    res
}



tmp<-fastRead(paste0(inputDir,"CP_featureMatFiltered.tsv"))
selectedFeat <- rn(tmp)
rm(tmp);gc()


cellAnnotRun1<-fastRead(paste0(inputDir,"CPmerge_cellAnnotRun1.tsv"))
cellAnnotRun2<-fastRead(paste0(inputDir,"CPmerge_cellAnnotRun2.tsv"))
rawFeatMatRun1<-fastRead(paste0(inputDir,"CPmerge_rawFeatMatRun1.tsv"),as.matrix = TRUE)
rawFeatMatRun2<-fastRead(paste0(inputDir,"CPmerge_rawFeatMatRun2.tsv"),as.matrix = TRUE)

commonCols <- intersect(colnames(cellAnnotRun1), colnames(cellAnnotRun2))
commonColsFeat <- intersect(colnames(rawFeatMatRun1), colnames(rawFeatMatRun2))
rawFeatMat <- rbind(rawFeatMatRun1[,commonColsFeat], rawFeatMatRun2[,commonColsFeat])

cellAnnot <- rbind(cellAnnotRun1[,commonCols], cellAnnotRun2[,commonCols])
cellAnnot<- cellAnnot[cellAnnot$parentalCellLine != "pool",]
rawFeatMat <- rawFeatMat[rn(cellAnnot),]

cellAnnot$oldClusters <- cellAnnot$cellClusters

cellAnnot$isKabuki = FALSE
cellAnnot[cellAnnot$parentalCellLine %in% c('aask4','ierp4','oadp4','qeti2'),"isKabuki"] = TRUE

######
regularizationDat<-fastRead(paste0(inputDir,"CP_FeatureRegularization.tsv"))
regu<-regularizationDat$transformation;names(regu)<-rn(regularizationDat)
rm(regularizationDat)

badcells <- scan(paste0(inputDir,"CPmerge_cell2remove.txt"), what = "character") #bad quality cells to remove
cellAnnot <- cellAnnot[!rownames(cellAnnot) %in% badcells,]
rawFeatMat <- rawFeatMat[!rownames(rawFeatMat) %in% badcells,]


featureMat <- rawFeatMat

#Put everything on the range 0-1000
featureMat<-apply(featureMat,1,function(x){
    x<-x-min(x) #min = 0
    x<-x/max(x) * 1000
    x
}) |> t()

rm(rawFeatMatRun1, rawFeatMatRun2, cellAnnotRun1, cellAnnotRun2, rawFeatMat)


#Per feature type transformation

for(f in names(regu)[regu=="log2(x+1)"]){
    featureMat[,f]<-log2(featureMat[,f]+1)
    featureMat[,f]<-featureMat[,f]/max(featureMat[,f])*1000 #put again on range 0-1000
}

for(f in names(regu)[regu=="log2(x+0.01)"]){
    x<-featureMat[,f]
    x<-log2(x+0.01)
    x<-x - min(x)
    featureMat[,f]<-x/max(x)*1000
    #qplotDensity(featureMat[,f])
}

for(f in names(regu)[regu=="log2(max(x)-x+0.1)"]){
    x<- featureMat[,f]
    x<- log2(max(x)-x+0.1) 
    x<- max(x)-x
    x<- x - min(x)
    featureMat[,f]<-x/max(x)*1000
    #qplotDensity(featureMat[,f])
}

for(f in names(regu)[regu=="log2(max(x)-x+0.5)"]){
    x<- featureMat[,f]
    x<- log2(max(x)-x+0.5) 
    x<- max(x)-x
    x<- x - min(x)
    featureMat[,f]<-x/max(x)*1000
    #plot(qplotDensity(x/max(x)*1000,returnGraph = T)+ggtitle(f))
}

for(f in names(regu)[regu=="log2(max(x)-x+0.01)"]){
    x<- featureMat[,f]
    x<- log2(max(x)-x+0.01) 
    x<- max(x)-x
    x<- x - min(x)
    featureMat[,f]<-x/max(x)*1000
    #plot(qplotDensity(x/max(x)*1000,returnGraph = T)+ggtitle(f))
}

for(f in names(regu)[regu=="2{ log2(max(x)-x+1)}"]){
    x<- featureMat[,f]
    x<- log2(max(x)-x+1)
    x<- log2(max(x)-x+.5)
    featureMat[,f]<-x
    #qplotDensity(x)
}

featureMat <- featureMat[,selectedFeat]

resHarmony <- harmony::RunHarmony(featureMat |> t(), cellAnnot$run, verbose = TRUE,lambda=1)

## 
pca<-fastPCA(resHarmony[cellAnnot$run =="run1",],transpose = F ,nPC = 50)
pcaAll <- oob::pcaAddSamples(pca, resHarmony[cellAnnot$run =="run2",], transpose = F)
umapRes <- UMAP(pcaAll$x[,3:50],transpose = FALSE, n_neighbors = 20, ret_nn = T)


cellAnnotD20 <- cellAnnot[cellAnnot$diffDay == "d20",]
sampleDrawD20<-rn(cellAnnot)[cellAnnot$diffDay == "d20"]
sampleDrawD20 <- sample(sampleDrawD20, length(sampleDrawD20))
umapd20 <- UMAP(pcaAll$x[cellAnnot$diffDay == "d20",1:50],transpose = FALSE, n_neighbors = 20)


pdf("out2.3/umapRun1_2_day20.pdf",width = 10,height = 10)
proj2d(umapRes$embedding[sampleDrawD20,], colorBy = cellAnnot[sampleDrawD20,]$run, useScatterMore = T, legendTitle = "run")
proj2d(umapRes$embedding[sampleDrawD20,], colorBy = cellAnnot[sampleDrawD20,]$isKabuki, useScatterMore = T, legendTitle = "isKabuki")
proj2d(umapRes$embedding[sampleDrawD20,], colorBy = cellAnnot[sampleDrawD20,]$diffDay, useScatterMore = T, legendTitle= "diffDay")
proj2d(umapRes$embedding[sampleDrawD20,], colorBy = cellAnnot[sampleDrawD20,]$parentalCellLine, useScatterMore = T, legendTitle= "parentalCellLine")
proj2d(umapRes$embedding[sampleDrawD20,], colorBy = cellAnnot[sampleDrawD20,]$oldClusters, useScatterMore = T, legendTitle= "oldClusters",
           colorScale = colorScaleCluster)
dev.off()

cellLines <- cellAnnot$parentalCellLine |> unique()
expGroups <- ifelse(cellLines %in% c('aask4','ierp4','oadp4','qeti2'),"Kabuki","WildType")

sampleAnnot<-data.frame(expGroups=expGroups,row.names = cellLines)


#####
library(FNN)

umap_R1 <- UMAP(resHarmony[cellAnnot$run=="run1",],transpose = FALSE, n_neighbors = 20, n_components = 10)
umap_R2 <- UMAP(resHarmony[cellAnnot$run=="run2",],transpose = FALSE, n_neighbors = 20, n_components = 10)

knn_result <- knn(resHarmony[cellAnnot$run=="run1",], resHarmony[cellAnnot$run=="run2",], cl = cellAnnot$oldClusters[!is.na(cellAnnot$oldClusters)], k = 5)
knn_distances <- attr(knn_result, "nn.dist")


g<-ggplotly(qplotDensity(apply(knn_distances,1,min),returnGraph = T))
print_app(g)


cellAnnot$reannotatedClusters <- cellAnnot$oldClusters
distance_threshold <- 400

dataset2_clusters <- sapply(1:sum(cellAnnot$run=="run2"), function(i) {
    min_distance <- min(knn_distances[i, ])  
    if (min_distance < distance_threshold) {
        return(as.character(knn_result[i]) )
    } else {
        return(NA)
    }
})

# Update cluster_annotation for dataset 2
cellAnnot$reannotatedClusters[is.na(cellAnnot$oldClusters)] <- dataset2_clusters


###diff abundance reannotated clusters
cellAnnotD20$reannotatedClusters <- cellAnnot$reannotatedClusters[cellAnnot$diffDay == "d20"]
table_reannot_run <- table(cellAnnotD20$reannotatedClusters,cellAnnotD20$run,useNA = "ifany")
heatmap.DM(table_reannot_run, preSet = "vanilla", showValues = T, name = "Cluster abundance")

ggData <- data.frame(table_reannot_run)
colnames(ggData) <- c("cluster","run","count")

ggplot(ggData,aes(x=cluster,fill=run,y=count))+
    geom_bar(stat="identity")+
    scale_fill_manual(values = oob::oobColors(2))+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave("out2.3/clusterAbundanceReannotatedS20.pdf",width = 7,height = 7)

diffAbundances <- list()
for(day in c("20","40","70")){
    dday <- paste0("d",day)
    diffAb<-diffAbundance(cellAnnot$reannotatedClusters[cellAnnot$diffDay == dday],
                                    cellAnnot$parentalCellLine[cellAnnot$diffDay ==dday],
                                    sampleAnnot)
    fastWrite(data.frame(diffAb$table),paste0("out2.3/diffAbundanceD",day,".tsv"))
    
    ggData <- data.frame(topTags(diffAb))
    ggData$cluster <- rn(ggData)
    
    pdf(paste0("out2.3/diffAbundanceD",day,".pdf"),width = 7,height = 7)
    print(
        volcanoPlot(ggData,effectSizeCol = "logFC", adjPvalCol = "FDR",labelCol = "cluster",returnGraph = T)+
        ggtitle(paste0("Differential abundance between Kabuki and WildType cell lines, day",day))+
        geom_hline(yintercept = 0.05,color = "red")
    )
    dev.off()
    
    diffAbundances[[dday]] <- diffAb
    
    pdf(paste0("out2.3/diffAbundanceD",day,"_heatmap.pdf"),width = 11,height = 9)
    heatmap.DM(diffAb$abundance,name = "abundance", preSet = "cor",showValues = T,useProb = F)
    heatmap.DM(diffAb$fitted.values,name = "fitted.values", preSet = "cor",showValues = T,useProb = F)
    dev.off()
}
