library(Seurat)
library(ggpubr)
library(gridExtra)
library(grid)
library(gridtext)
library(scater)
library(scales)
library(readxl)
library(ggplot2)
library(dplyr)
library(harmony)
library(GenomicRanges)
library(patchwork)

args = commandArgs(trailingOnly=TRUE)
queryRed <- args[1]
referenceRed <- args[2]
dataset <- args[3]


if (referenceRed=="pca"){
  redUmap="umap.pca"
  print(redUmap)
} else if (referenceRed=="harmony"){
  redUmap="umap.harmony"
  print(redUmap)
}

query <- readRDS("seuratObjects/mapRef/regOutCellPlexAndProtocol.RDS")

if (dataset=="Polioudakis"){
  print(dataset)
  reference <- readRDS("seuratObjects/mapRef/regOutPolioudakis.RDS")
  celltypeMapping <- "Cluster"
} else if (dataset=="BhaduriOrganoids") {
  print(dataset)
  reference <- readRDS("seuratObjects/mapRef/regOutBhaduriOrganoids.RDS")
  celltypeMapping <- "typeSubType"
}

query <- RunPCA(query, features = VariableFeatures(object = query))
reference <- RunPCA(reference, features = VariableFeatures(object = reference))


#############
### Query ###
#############

dirProcessing <- "seuratObjects/mapRef/harmonyResults/"

#Pre-compute batch-corrected PCA for query even though might not be used later
PCAbeforHarmony <- DimPlot(query, reduction = "pca", pt.size = .1, group.by = 'sampleId')+ggtitle("Before Harmony")
set.seed(123)
query <- query %>%
  RunHarmony(c("sampleId"), plot_convergence = TRUE)
harmony_embeddings <- Embeddings(query, 'harmony')
harmony_embeddings[1:5, 1:5]
PCAafterHarmony <- DimPlot(query, reduction = "harmony", pt.size = .1, group.by = 'sampleId')+ggtitle("After Harmony")
comparisonPCA <- PCAbeforHarmony + PCAafterHarmony +
  plot_layout(guides="collect")
pdf(file=paste0(dirProcessing,"queryPCAcomparison_",dataset,".pdf"), width=6, height = 4)
plot(comparisonPCA)
dev.off()


if (queryRed=="pca"){
  print(queryRed)
  query <- query %>%
    FindNeighbors(reduction = "pca", dims = 1:15) %>%
    FindClusters(resolution = 0.8) %>%
    identity()
  query <- RunUMAP(query,  dims = 1:15, reduction = queryRed, reduction.name = redUmap, reduction.key = "UMAP_", return.model = T)
} else if (queryRed=="harmony") {
  print(queryRed)
  query <- query %>%
    FindNeighbors(reduction = "harmony", dims = 1:15) %>%
    FindClusters(resolution = 0.8) %>%
    identity()
  query <- RunUMAP(query, dims = 1:15, reduction = queryRed, reduction.name = redUmap, reduction.key = "UMAPhar_", return.model = T)
}


#################
### Reference ###
#################


#### group by might change in Bhaduri otherwise error !!

if (dataset=="Polioudakis"){
  
  print(dataset)
  PCAbeforHarmony <- DimPlot(reference, reduction = "pca", pt.size = .1, group.by = 'Library')+ggtitle("Before Harmony")
  set.seed(123)
  reference <- reference %>%
    RunHarmony(c("Library"), plot_convergence = TRUE)
  harmony_embeddings <- Embeddings(reference, 'harmony')
  PCAafterHarmony <- DimPlot(reference, reduction = "harmony", pt.size = .1, group.by = 'Library')+ggtitle("After Harmony")
  comparisonPCA <- PCAbeforHarmony + PCAafterHarmony +
    plot_layout(guides="collect")
  pdf(file=paste0(dirProcessing,"referencePCAcomparison_",dataset,".pdf"), width=6, height = 4)
  plot(comparisonPCA)
  dev.off()

  
} else if (dataset=="BhaduriOrganoids"){
  
  print(dataset)
  PCAbeforHarmony <- DimPlot(reference, reduction = "pca", pt.size = .1, group.by = 'Sample')+ggtitle("Before Harmony")
  set.seed(123)
  reference <- reference %>%
    RunHarmony(c("Sample"), plot_convergence = TRUE)
  harmony_embeddings <- Embeddings(reference, 'harmony')
  PCAafterHarmony <- DimPlot(reference, reduction = "harmony", pt.size = .1, group.by = 'Sample')+ggtitle("After Harmony")
  comparisonPCA <- PCAbeforHarmony + PCAafterHarmony +
    plot_layout(guides="collect")
  pdf(file=paste0(dirProcessing,"referencePCAcomparison_",dataset,".pdf"), width=6, height = 4)
  plot(comparisonPCA)
  dev.off()
  
}



if (referenceRed=="pca"){
  
  print(referenceRed)
  reference <- reference %>%
    FindNeighbors(reduction = "pca", dims = 1:15) %>%
    FindClusters(resolution = 0.8) %>%
    identity()
  reference <- RunUMAP(reference,  dims = 1:15, reduction = referenceRed, n.neighbors = 10, reduction.name = redUmap, reduction.key = "UMAP_", return.model = T)
  
} else if (referenceRed=="harmony") {
  
  print(referenceRed)
  reference <- reference %>%
    FindNeighbors(reduction = "harmony", dims = 1:15) %>%
    FindClusters(resolution = 0.8) %>%
    identity()
  reference <- RunUMAP(reference, dims = 1:15, reduction = referenceRed, n.neighbors = 10, reduction.name = redUmap, reduction.key = "UMAPhar_", return.model = T)
  
}


#############################
### Mapping the reference ###
#############################


if (referenceRed=="harmony"){
  
  print(paste0("DimRedObjectModified:", referenceRed))
  reference[['harmony2']] <- CreateDimReducObject(embeddings = reference[['harmony']]@cell.embeddings,
                                                  key = "harmony2_",
                                                  loadings = reference[['pca']]@feature.loadings,
                                                  assay = "RNA")
}


if (queryRed=="harmony"){
  
  print(paste0("DimRedObjectModified:", queryRed))
  query[['pca']] <- CreateDimReducObject(embeddings = query[['harmony']]@cell.embeddings,
                                         key = "PC_",
                                         loadings = query[['pca']]@feature.loadings,
                                         assay = "RNA")
  
}


if (referenceRed=="harmony"){
  
  print(paste0("FindTransferAnchors:", referenceRed))
  anchor  <- FindTransferAnchors(reference = reference,
                                 reference.reduction = "harmony2",
                                 query = query ,
                                 dims = 1:30)
  
  query <- MapQuery(anchorset = anchor,
                    query = query,
                    reference = reference,
                    refdata = list(celltype = celltypeMapping),
                    reduction.model = redUmap)
  
} else if (referenceRed=="pca"){
  
  print(paste0("FindTransferAnchors:", referenceRed))
  anchor  <- FindTransferAnchors(reference = reference,
                                 reference.reduction = "pca",
                                 query = query,
                                 dims = 1:30)
  
  query <- MapQuery(anchorset = anchor,
                    query = query,
                    reference = reference,
                    refdata = list(celltype = celltypeMapping),
                    reduction.model = redUmap)
  
}


dfResults <- data.frame(robustID=paste0(gsub("-.+","",colnames(query)), "_", query$sampleId),
                        predicted.celltype=query$predicted.celltype,
                        predicted.celltype.score=query$predicted.celltype.score,
                        mapRef=paste0(dataset,"_Q",queryRed,"_R",referenceRed))

print("Percentage of cells with predicted.cell.type.score > 0.5")
print(round(table(dfResults$predicted.celltype.score>0.5)["TRUE"]*100/sum(table(dfResults$predicted.celltype.score>0.5)),2))

saveRDS(dfResults, file=paste0("seuratObjects/mapRef/tabs/",dataset,"_Q",queryRed,"_R",referenceRed,".RDS"))





























