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

query <- readRDS("seuratObjects/rawCountsCorticalDiff70D.RDS")

## normalise
query <- NormalizeData(query, normalization.method = "LogNormalize", scale.factor = 10000)

## Find highly variable features
query <- FindVariableFeatures(query, selection.method = "vst", nfeatures = 3000)

# Read in the expression matrix The first row is a header row, the first column is rownames
pathToFile="extdata/cell_cycle_vignette_files/nestorawa_forcellcycle_expressionMatrix.txt"
exp.mat <- read.table(file = pathToFile,
                      header = TRUE,
                      as.is = TRUE, row.names = 1)

# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

query <- CellCycleScoring(query, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
query@meta.data$old.ident <- NULL

query$CC.Difference <- query$S.Score - query$G2M.Score
query <- ScaleData(query,
                      vars.to.regress = "CC.Difference",
                      features = rownames(query))


saveRDS(query, "seuratObjects/mapRef/regOutCellPlexAndProtocol.RDS")



