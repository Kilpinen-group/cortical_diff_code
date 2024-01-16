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

reference <- readRDS("seuratObjects/rawCountsBhaduriOrganoids.RDS")

## following article processing steps
reference[["percent.mt"]] <- PercentageFeatureSet(reference, pattern = "^MT-")
reference <- subset(reference, subset = nFeature_RNA > 500 & percent.mt < 10)
reference <- NormalizeData(reference, normalization.method = "LogNormalize", scale.factor = 10000)
reference <- FindVariableFeatures(reference, selection.method = "vst", nfeatures = 3000)

# Read in the expression matrix The first row is a header row, the first column is rownames
pathToFile="extdata/cell_cycle_vignette_files/nestorawa_forcellcycle_expressionMatrix.txt"
exp.mat <- read.table(file = pathToFile,
                      header = TRUE,
                      as.is = TRUE, row.names = 1)

# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

reference <- CellCycleScoring(reference, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
reference@meta.data$old.ident <- NULL

reference$CC.Difference <- reference$S.Score - reference$G2M.Score

reference <- ScaleData(reference,
                       vars.to.regress="CC.Difference",
                       features = rownames(reference))

saveRDS(reference, "seuratObjects/mapRef/regOutBhaduriOrganoids.RDS")



