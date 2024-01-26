## d40_70

library(Seurat)
library(readxl)
library(ggplot2)
library(dplyr)
library(harmony)

### change PATH ###
rootDir <- "d40_70/"
###################

summfilesFile <- paste0(rootDir,"d40d70_summaryFiles.txt")
summfiles <- read.table(summfilesFile, header=F)

qcTab <- sapply(summfiles$V1, function(x){
  tt <- read.csv(x)
  tt$analysis <- "LiveseyModified"
  tt
}, simplify=F)

qcTab <- do.call("rbind", qcTab)
rownames(qcTab) <- NULL
qcTab$file <- summfiles$V1
qcTab$sample <- sapply(strsplit(qcTab$file,"/"), function(x) x[10])


metadata <- read_excel(path = paste0(rootDir, "d40d70_infoDesign.xlsx"),
                       sheet = "Design")
metadata <- as.data.frame(metadata)
metadata <- subset(metadata, QC=="Passed")
metadata$QC <- NULL
qcTab <- merge(metadata, qcTab, by.x="10x Samples", by.y="sample")


## Read CellRanger output
dir10xfiles <- paste0(rootDir, "data/count_211123_A00464_0418_BHNVTFDSX2/")
list_10xfiles <- paste0(dir10xfiles,dir(dir10xfiles)[grepl("Pool",dir(dir10xfiles))],"/outs/filtered_feature_bc_matrix/")

perc <- sapply(list_10xfiles, function(x){
  
  ## Load Demuxlet output (it's generated in "count_211123_A00464_0418_BHNVTFDSX2/*/outs/" directory)
  demuxlet_file <- dir(gsub("/filtered_feature_bc_matrix/","", x), full.names=T)[grepl(".best",dir(gsub("/filtered_feature_bc_matrix/","", x)))]
  demuxlet <- read.table(demuxlet_file, header=T, sep="\t")
  singletons_perc <- unname(round(table(grepl("SNG-", demuxlet$BEST))*100/sum(table(grepl("SNG-", demuxlet$BEST))),2)["TRUE"])
  
  tmp <- data.frame(sample=sapply(strsplit(demuxlet_file, "/"), function(x) x[10]),
                    percSingletons=singletons_perc)
  print(x)
  print(table(demuxlet[grepl("SNG", demuxlet$BEST),]$BEST))
  
  return(tmp)
  
}, simplify=F)

perc <- do.call("rbind", perc)
rownames(perc) <- NULL

qcTab <- merge(qcTab, perc, by.x="10x Samples", by.y="sample")


processing <- function(return.obj=T, mito=F){
  
  allobj <- sapply(list_10xfiles, function(x){
    
    tenXrun <- Read10X(data.dir = x)
    tenXrun <- CreateSeuratObject(counts = tenXrun, project = "d40_70")
    sample <- sapply(strsplit(x,"/"), function(x) x[10])
    
    tenXrun@meta.data$sampleId <-  sample
    tenXrun@meta.data$differentiationDay <- qcTab[match(sample, qcTab$`10x Samples`),]$differentiationDay
    
    demuxlet_file <- dir(gsub("/filtered_feature_bc_matrix/","", x), full.names=T)[grepl(".best",dir(gsub("/filtered_feature_bc_matrix/","", x)))]
    demuxlet <- read.table(demuxlet_file, header=T, sep="\t")
    
    tenXrun <- tenXrun[,which(!is.na(match(rownames(tenXrun@meta.data), demuxlet$BARCODE)))]
    tenXrun <- tenXrun[,which(grepl("SNG-", demuxlet$BEST))]
    demuxlet <- demuxlet[grepl("SNG-", demuxlet$BEST),]
    
    stopifnot(dim(demuxlet)[1]==dim(tenXrun)[2])
    
    tenXrun@meta.data$donorId <- gsub("SNG-","",demuxlet[match(rownames(tenXrun@meta.data), demuxlet$BARCODE),]$BEST)

    metadata <- read_excel(path = paste0(rootDir, "d40d70_infoDesign.xlsx"),
                           sheet = "Pool")
    metadata <- as.data.frame(metadata)
    metadata <- subset(metadata, Sample==sample)
    
    tenXrun@meta.data$Condition <- metadata[match(tenXrun@meta.data$donorId, metadata$iPSC_line),]$Condition
    
    limsNFeatures <- c(2000,8000)
    limsMito <- 15
    
    idDef <- unname(sapply(strsplit(x,"/"), function(y) y[11]))
    tenXrun[["percent.mt"]] <- PercentageFeatureSet(tenXrun, pattern = "^MT-")
    
    plot1 <- FeatureScatter(tenXrun, feature1 = "nCount_RNA", feature2 = "percent.mt")+
      ggtitle(idDef)+theme(plot.title=element_text(hjust=0.5, size=11),
                           legend.position="none")
    plot2 <- FeatureScatter(tenXrun, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")+
      ggtitle(idDef)+
      theme(plot.title=element_text(hjust=0.5, size=11),
            legend.position="none")
    qcplot1 <- plot1 + plot2
    
    tenXrun[["percent.ribo"]]<- PercentageFeatureSet(tenXrun, pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA")
    ##QC-ribosomal content
    stopifnot(all(mean(tenXrun@meta.data$percent.ribo)>10 | mean(tenXrun@meta.data$percent.ribo)<20))
    
    dirQCplots <- paste0(rootDir,"QC/")
    pdf(file=paste0(dirQCplots,"QCplot1_",idDef, ".pdf"))
    plot(qcplot1)
    dev.off()
    
    tmpTabQCFeature <- data.frame(nFeature_RNA=unname(sort(tenXrun$nFeature_RNA)),
                                  index=1:length(sort(tenXrun$nFeature_RNA)),
                                  sample=unname(sample))
    
    tmpTabQCMito <- data.frame(percent.mt=unname(sort(tenXrun$percent.mt)),
                               index=1:length(sort(tenXrun$percent.mt)),
                               sample=sample)

    
    tenXrun <- subset(tenXrun, subset = nFeature_RNA > limsNFeatures[1] & nFeature_RNA < limsNFeatures[2] & percent.mt < limsMito)
    
    if (return.obj==TRUE){
      return(tenXrun)
      
    } else {
      
      if (mito==TRUE){
        
        return(tmpTabQCMito)
        
      } else {
        
        return(tmpTabQCFeature)
        
      }
    }
    
    
  }, simplify=F)
  
  
  return(allobj)
  
  
}

allobj <- processing(return.obj=T, mito=F)
nCountQC <- processing(return.obj=F, mito=F)
mitoQC <- processing(return.obj=F, mito=T)

nCountQC <- do.call("rbind", nCountQC)
rownames(nCountQC) <- NULL
mitoQC <- do.call("rbind", mitoQC)
rownames(mitoQC) <- NULL


saveRDS(allobj, file="seuratObjects/CorticalDiff_d40_70.RDS")