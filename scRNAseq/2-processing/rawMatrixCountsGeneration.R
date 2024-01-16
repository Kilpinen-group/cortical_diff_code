
library(Seurat)
library(wesanderson)
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

#######################################
## Query: CorticalDiff (d20-d40-d70) ##
#######################################


seuratObjects <- "seuratObjects/"
seuratObjects <- paste0(seuratObjects,dir(seuratObjects))
seuratObjects <- c(seuratObjects[grepl("ModLiveseyRosa.RDS", seuratObjects)],
                   seuratObjects[grepl("CellPlexAdithi.RDS", seuratObjects)])

allobj <- sapply(seuratObjects, function(x){
  
  tmp <- readRDS(x)
  
}, simplify=F)

allobj <- unlist(allobj)
names(allobj) <- gsub(".+RDS.","", names(allobj))

maskAdithi <- grepl("BHVYHGDSX2", names(allobj))
names(allobj)[maskAdithi] <- 
  sapply(strsplit(names(allobj)[maskAdithi],"/"), function(x) paste0(x[10],"/",x[13]))

maskRosa <- grepl("rosaDataset", names(allobj))
names(allobj)[maskRosa] <- sapply(strsplit(names(allobj)[maskRosa], "/"), function(x) x[10])

mlo.big <- merge(allobj[[1]], y=do.call("c",allobj[2:length(allobj)]),
                 project="OTCorticalNeurons")
rm(allobj)

mlo.big@meta.data$donorId_simplified <- mlo.big@meta.data$donorId
mask_notHipsci <- !grepl("HPSI", mlo.big@meta.data$donorId_simplified)
mlo.big@meta.data$donorId_simplified[mask_notHipsci] <- gsub(".+KIL.", "", mlo.big@meta.data$donorId_simplified[mask_notHipsci])
mlo.big@meta.data$donorId_simplified[mask_notHipsci] <- gsub("-",".", sapply(strsplit(mlo.big@meta.data$donorId_simplified[mask_notHipsci],"_"), function(x) x[1]))


## Add biological and technical replicates info

mlo.big@meta.data[mlo.big@meta.data$orig.ident=="RosaLiveseyMod",]$orig.ident <- "AtlasProject"

mlo.big@meta.data$bioRep <- NA
mlo.big@meta.data$techRep <- NA
mlo.big@meta.data$Protocol <- NA

metadata <- as.data.frame(read_excel(path="/lustre/scratch126/cellgen/kilpinen/pp9/scRNA-seq/OpenTargetsNPC/OT_pilotAndCellplex.xlsx",
                                     sheet="rosaProtocol"))
maskRosa <- mlo.big@meta.data$orig.ident=="AtlasProject"
mlo.big@meta.data$Protocol[maskRosa] <- metadata[match(paste0(mlo.big@meta.data[maskRosa,]$sampleId,"-",mlo.big@meta.data[maskRosa,]$donorId_simplified),
                                                       paste0(metadata$library10x,"-",metadata$line)),]$protocol
mlo.big@meta.data$techRep[maskRosa] <- metadata[match(paste0(mlo.big@meta.data[maskRosa,]$sampleId,"-",mlo.big@meta.data[maskRosa,]$donorId_simplified),
                                                      paste0(metadata$library10x,"-",metadata$line)),]$techRep



metadata <- as.data.frame(read_excel(path="/lustre/scratch126/cellgen/kilpinen/pp9/scRNA-seq/OpenTargetsNPC/OT_pilotAndCellplex.xlsx",
                                     sheet="Adithi"))
maskAdithi <- mlo.big@meta.data$orig.ident=="Cellplex_Adithi"
mlo.big@meta.data$Protocol[maskAdithi] <- metadata[match(paste0(mlo.big@meta.data[maskAdithi,]$sampleId,"-",mlo.big@meta.data[maskAdithi,]$donorId),
                                                         paste0(metadata$library10x,"-",metadata$line)),]$protocol
mlo.big@meta.data$techRep[maskAdithi] <- metadata[match(paste0(mlo.big@meta.data[maskAdithi,]$sampleId,"-",mlo.big@meta.data[maskAdithi,]$donorId),
                                                        paste0(metadata$library10x,"-",metadata$line)),]$techRep
mlo.big@meta.data$bioRep[maskAdithi] <- metadata[match(paste0(mlo.big@meta.data[maskAdithi,]$sampleId,"-",mlo.big@meta.data[maskAdithi,]$donorId),
                                                       paste0(metadata$library10x,"-",metadata$line)),]$bioRep

## genes expressed in less than 0.1% total cells are removed
selected_f <- names(which(rowSums(mlo.big@assays$RNA[]>0)>0.001*dim(mlo.big)[1]))
mlo.big <- subset(mlo.big, features = selected_f)
mlo.big

saveRDS(mlo.big, "seuratObjects/rawCountsCorticalDiff70D.RDS")


###########################
## Polioudakis reference ##
###########################

## Data download: http://solo.bmap.ucla.edu/shiny/webapp/

pathToFolder <- "reference/codex/"
load(paste0(pathToFolder,"raw_counts_mat.rdata"))
meta <- read.csv(paste0(pathToFolder,"cell_metadata.csv"), header=T, row.names=1)
reference <- CreateSeuratObject(counts = raw_counts_mat, project = "Polioudakis", meta.data=meta)
reference <- reference[,!is.na(reference$Number_UMI)]
selected_f <- names(which(rowSums(reference@assays$RNA[]>0)>0.001*dim(reference)[1]))
reference <- subset(reference, features = selected_f)
saveRDS(reference, "seuratObjects/rawCountsPolioudakis2019.RDS")


################################
## Bhaduri organoid reference ##
################################

## Data download: https://cells.ucsc.edu/?ds=organoidreportcard

library(data.table)
library(Seurat)
library(readxl)

pathToFolder <- "reference/BhaduriOrganoids/"
mat <- fread(paste0(pathToFolder,"organoidsmatrix_nonorm.txt.gz"))
meta <- read.table(paste0(pathToFolder,"meta.tsv"), header=T, sep="\t", as.is=T, row.names=1)
meta$rnamesMeta <- rownames(meta)
genes = mat[,1][[1]]
mat <- as.matrix(mat[,-1])
rownames(mat) <- genes

## Bhaduri sent by mail the keys matching for reproducibility
matchIds <- as.data.frame(read_excel("reference/BhaduriOrganoids/Keys_Bhadurietal_Nature_2020_Organoid_Cells_Matrix_to_Metadata.xlsx"))
colnames(matchIds) <- c("rnamesMeta","cnamesMat")

dim(matchIds)[1]
# [1] 242349

dim(mat)[2]
# [1] 341492

dim(meta)[1]
# [1] 235121

## Not all cell ids mappings are present in the matrix of expression, but not with metadata table.
table(is.na(match(matchIds$cnamesMat, colnames(mat))))
# FALSE   TRUE 
# 238545   3804

## Not all cell ids from metadata are present in the cell ids mapping table.
table(is.na(match(rownames(meta), matchIds$rnamesMeta)))
# FALSE   TRUE 
# 194642  40479 


## Only 190,838 cells from the raw counts matrix have a corresponding annotation with the metadata table.
meta <- merge(meta, matchIds, by="rnamesMeta")
table(meta$cnamesMat %in% colnames(mat))
# FALSE   TRUE 
# 3804 190838

meta <- subset(meta, iPSCorhESC=="iPSC")
meta$typeSubType <- paste0(meta$Type,"-", meta$Subtype)
meta$typeSubType <- gsub("Astrocyte-Astrocyte","Astro", meta$typeSubType)
meta$typeSubType <- gsub("Unknown-Unknown","Unk", meta$typeSubType)
meta$typeSubType <- gsub("ExcitatoryNeuron-DeepLayer","ExDp", meta$typeSubType)
meta$typeSubType <- gsub("ExcitatoryNeuron-UpperLayer","ExU", meta$typeSubType)
meta$typeSubType <- gsub("IPC-MatureIPC","IPC-Mature", meta$typeSubType)
meta$typeSubType <- gsub("Astrocyte-hindbrainAstrocyte","AstroHindb", meta$typeSubType)
meta$typeSubType <- gsub("ExcitatoryNeuron-glycolyticneurons","ExNeuroGlyco", meta$typeSubType)
meta$typeSubType <- gsub("ExcitatoryNeuron-Newborn","ExNeuroNew", meta$typeSubType)
meta$typeSubType <- gsub("ExcitatoryNeuron-panNeuron","ExPanNeuro", meta$typeSubType)
meta$typeSubType <- gsub("InhibitoryNeuron-Interneuron","InN", meta$typeSubType)
meta$typeSubType <- gsub("Outlier-Outlier","Outlier", meta$typeSubType)
meta$typeSubType <- gsub("RadialGlia-earlyRG","earlyRG", meta$typeSubType)
meta$typeSubType <- gsub("RadialGlia-hindbrainRG","hRG", meta$typeSubType)
meta$typeSubType <- gsub("RadialGlia-lowquality","LQ-RG", meta$typeSubType)
meta$typeSubType <- gsub("RadialGlia-panRG","panRG", meta$typeSubType)
meta$typeSubType <- gsub("RadialGlia-glycolyticRG","glycoRG", meta$typeSubType)


mat <- mat[,colnames(mat) %in% rownames(meta)]
gc()


reference <- CreateSeuratObject(counts = mat, project = "BhaduriOrganoids", meta.data=meta)
rm(mat)
rm(meta)
gc()

reference <- reference[, !(reference$Type=="Outlier")]
selected_f <- names(which(rowSums(reference@assays$RNA[]>0)>0.001*dim(reference)[1]))
reference <- subset(reference, features = selected_f)


saveRDS(reference, "seuratObjects/emilia/rawCountsBhaduriOrganoids.RDS")


