

## Cellplex_d20

library(Seurat)
library(readxl)
library(ggplot2)
library(dplyr)
library(harmony)
library(GenomicRanges)


### change PATH ###
rootDir <- "cellplex/"
###################

summfilesFile <- paste0(rootDir,"d20_summaryFiles.txt")
summfiles <- read.table(summfilesFile, header=F)

qcTab <- sapply(summfiles$V1, function(x){
  tt <- read.csv(x)
  
  tmp <- data.frame("Estimated.number.of.cells"=subset(tt, Metric.Name=="Estimated number of cells" & Library.Type=="Gene Expression" & Group.Name=="GEX_1")$Metric.Value,
                    "Mean.Reads.per.Cell"=subset(tt, Metric.Name=="Mean reads per cell" & Library.Type=="Gene Expression" & Group.Name=="GEX_1")$Metric.Value,
                    "Median.Genes.per.Cell"=subset(tt, Metric.Name=="Median genes per cell" & Library.Type=="Gene Expression")$Metric.Value,
                    "Number.of.Reads"=subset(tt, Metric.Name=="Number of reads" & Library.Type=="Gene Expression" & Group.Name=="GEX_1")$Metric.Value,
                    "Valid.Barcodes"=subset(tt, Metric.Name=="Valid barcodes" & Library.Type=="Gene Expression" & Group.Name=="GEX_1")$Metric.Value,
                    "Sequencing.Saturation"=subset(tt, Metric.Name=="Sequencing saturation" & Library.Type=="Gene Expression" & Group.Name=="GEX_1")$Metric.Value,
                    "Q30.Bases.in.Barcode"=subset(tt, Metric.Name=="Q30 barcodes" & Library.Type=="Gene Expression")$Metric.Value,
                    "Q30.Bases.in.RNA.Read"=subset(tt, Metric.Name=="Q30 RNA read" & Library.Type=="Gene Expression")$Metric.Value,
                    "Q30.Bases.in.UMI"=subset(tt, Metric.Name=="Q30 UMI" & Library.Type=="Gene Expression")$Metric.Value,
                    "Reads.Mapped.to.Genome"= subset(tt, Metric.Name=="Mapped to genome" & Library.Type=="Gene Expression"& Group.Name=="GEX_1")$Metric.Value,
                    "Reads.Mapped.Confidently.to.Genome"= subset(tt, Metric.Name=="Confidently mapped to genome" & Library.Type=="Gene Expression"& Group.Name=="GEX_1")$Metric.Value,
                    "Reads.Mapped.Confidently.to.Intergenic.Regions"=subset(tt, Metric.Name=="Confidently mapped to intergenic regions" & Library.Type=="Gene Expression"& Group.Name=="GEX_1")$Metric.Value,
                    "Reads.Mapped.Confidently.to.Intronic.Regions"=subset(tt, Metric.Name=="Confidently mapped to intronic regions" & Library.Type=="Gene Expression"& Group.Name=="GEX_1")$Metric.Value,
                    "Reads.Mapped.Confidently.to.Exonic.Regions"=subset(tt, Metric.Name=="Confidently mapped to exonic regions" & Library.Type=="Gene Expression"& Group.Name=="GEX_1")$Metric.Value,
                    "Reads.Mapped.Confidently.to.Transcriptome"=subset(tt, Metric.Name=="Confidently mapped to transcriptome" & Library.Type=="Gene Expression"& Group.Name=="GEX_1")$Metric.Value,
                    "Reads.Mapped.Antisense.to.Gene"=subset(tt, Metric.Name=="Confidently mapped antisense" & Library.Type=="Gene Expression"& Group.Name=="GEX_1")$Metric.Value,
                    "Fraction.Reads.in.Cells"=subset(tt, Metric.Name=="Confidently mapped antisense" & Library.Type=="Gene Expression"& Group.Name=="GEX_1")$Metric.Value,
                    "Total.Genes.Detected"=subset(tt, Metric.Name=="Total genes detected" & Library.Type=="Gene Expression")$Metric.Value,
                    "Median.UMI.Counts.per.Cell"=subset(tt, Metric.Name=="Median UMI counts per cell" & Library.Type=="Gene Expression")$Metric.Value)

  tmp$analysis <- "Cellplex"
  tmp
}, simplify=F)

qcTab <- do.call("rbind", qcTab)
rownames(qcTab) <- NULL

qcTab$file <- gsub("/outs/.+","",summfiles$V1)
qcTab$sample <- paste0(sapply(strsplit(qcTab$file,"/"), function(x) x[10]))


tmp <- sapply(unique(qcTab$sample), function(y){
  
  temporal <- subset(qcTab, sample==y)
  
  temporal2 <- sapply(colnames(temporal), function(z){
    tmp3 <- unique(temporal[,z])
    if (length(tmp3)>1){
      mean(as.numeric(gsub("%","",gsub(",","",tmp3))))
    } else {
      tmp3
    }
    }, simplify=F)
  
  tmp2 <- as.data.frame(do.call("cbind", temporal2))
  tmp2

}, simplify=F)

tmp <- do.call("rbind", tmp)
rownames(tmp) <- NULL

qcTab <- tmp


### This requires the infoDesign.xlsx file !!!!!

metadata <- read_excel(path = paste0(rootDir, "d20_infoDesign.xlsx"),
                       sheet = "Design")
metadata <- as.data.frame(metadata)
metadata <- subset(metadata, QC=="Passed")
metadata$QC <- NULL
qcTab <- merge(metadata, qcTab, by.x="10x Samples", by.y="sample")


## this is in the summary!!!
##  Cells assigned to a sample 11516 (50.89%)
## Cells assigned to this sample    843 (3.73%)


returnSingletons <- function(summfiles, type="all"){
  
  if (type=="all"){
    
    percTab <- sapply(summfiles$V1, function(x){
      tt <- read.csv(x)
      perc <- tt[tt$Metric.Name=="Cells assigned to a sample" & tt$Library.Type=="Gene Expression" & tt$Group.Name=="GEX_1",]$Metric.Value
      perc <- as.numeric(gsub("%","",(unlist(regmatches(perc, gregexpr( "(?<=\\().+?(?=\\))", perc, perl = T))))))
      
      tmp <- data.frame(sample=sapply(strsplit(x, "/"), function(x) x[10]),
                        percSingletons=perc)
      
    }, simplify=F)
    
    percTab <- do.call("rbind",percTab[!duplicated(percTab)])
    rownames(percTab) <- NULL
    
  } else {
    

    percTab <- sapply(summfiles$V1, function(x){
      tt <- read.csv(x)
      totalCells <- tt[tt$Metric.Name=="Cells assigned to a sample" & tt$Library.Type=="Gene Expression" & tt$Group.Name=="GEX_1",]$Metric.Value
      totalCells <- as.numeric(gsub(",","",gsub(" \\(.+","",totalCells)))
      numCells <- as.numeric(gsub(",","",tt[tt$Metric.Name=="Cells",]$Metric.Value))
      tmp <- data.frame(sample=sapply(strsplit(x, "/"), function(x) x[13]),
                        numCells=numCells,
                        perc=signif(numCells*100/totalCells,2))
      tmp
    }, simplify=F)

    percTab <- do.call("rbind",percTab[!duplicated(percTab)])
    rownames(percTab) <- NULL

    
  }
  
  return(percTab)

}


perc <- returnSingletons(summfiles, type="all")
qcTab <- merge(qcTab, perc, by.x="10x Samples", by.y="sample")
individualSampleperc <- returnSingletons(summfiles, type="single")



## start processing!!!
list_10xfiles <- paste0(rootDir,"d20_list_10xfiles.txt")
list_10xfiles <- read.table(list_10xfiles)$V1

processing <- function(list_10xfiles, return.obj=T, mito=F){
  
  allobj <- sapply(list_10xfiles, function(x){
    print(paste0("Processing sample ",x))
    tenXrun <- Read10X(data.dir = x)[[1]]
    tenXrun <- CreateSeuratObject(counts = tenXrun, project = "Cellplex_Adithi")
    sample <- sapply(strsplit(x,"/"), function(x) x[10])
    
    tenXrun@meta.data$sampleId <-  sample
    tenXrun@meta.data$differentiationDay <- qcTab[match(sample, qcTab$`10x Samples`),]$differentiationDay
    
    metadata <- read_excel(path = paste0(rootDir, "d20_infoDesign.xlsx"),
                           sheet = "Design")
    metadata <- as.data.frame(metadata)
    metadata <- subset(metadata, `10x Samples`==sample)
    
    tenXrun@meta.data$donorId <- sapply(strsplit(x,"/"), function(x) x[13])
    
    metadata <- read_excel(path = paste0(rootDir, "d20_infoDesign.xlsx"),
                           sheet = "Pool")
    metadata <- as.data.frame(metadata)
    metadata <- subset(metadata, Sample==sample)
    
    tenXrun@meta.data$Condition <- metadata[match(tenXrun@meta.data$donorId, metadata$iPSC_line),]$Condition
    
    limsNFeatures <- c(2000,8000)
    limsMito <- 15
    
    idDef <- unname(sapply(strsplit(x,"/"), function(y) y[13]))
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
    pdf(file=paste0(dirQCplots,"QCplot_",idDef, ".pdf"))
    plot(qcplot1)
    dev.off()
    
    tmpTabQCFeature <- data.frame(nFeature_RNA=unname(sort(tenXrun$nFeature_RNA)),
                                  index=1:length(sort(tenXrun$nFeature_RNA)),
                                  sample=sapply(strsplit(x,"/"), function(x) paste0(x[10],"/",x[13])))
    
    tmpTabQCMito <- data.frame(percent.mt=unname(sort(tenXrun$percent.mt)),
                               index=1:length(sort(tenXrun$percent.mt)),
                               sample=sapply(strsplit(x,"/"), function(x) paste0(x[10],"/",x[13])))
    
    
    tenXrun <- subset(tenXrun, subset = nFeature_RNA > limsNFeatures[1] & nFeature_RNA < limsNFeatures[2] & percent.mt < limsMito)
    
    if (return.obj==TRUE){
      print(dim(tenXrun))
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

allobj <- processing(list_10xfiles, return.obj=T, mito=F)
nCountQC <- processing(list_10xfiles, return.obj=F, mito=F)
mitoQC <- processing(list_10xfiles, return.obj=F, mito=T)

nCountQC <- do.call("rbind", nCountQC)
rownames(nCountQC) <- NULL
mitoQC <- do.call("rbind", mitoQC)
rownames(mitoQC) <- NULL


saveRDS(allobj, file="seuratObjects/CellPlexAdithi.RDS")

## Probably do not need those RDS objects below, they were intended just for QC
#saveRDS(nCountQC, file="cellplex/aggr/nFeatures/CellPlexAdith_nFeat.RDS")
#saveRDS(mitoQC, file="cellplex/aggr/mitoCounts/CellPlexAdith_mito.RDS")
#saveRDS(qcTab, file="cellplex/aggr/qcTab/CellPlexAdith_qctab.RDS")


