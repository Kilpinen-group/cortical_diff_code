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

reference <- readRDS("seuratObjects/rawCountsPolioudakis2019.RDS")

## normalise
reference <- NormalizeData(reference, normalization.method = "LogNormalize", scale.factor = 10000)

## Find highly variable features
reference <- FindVariableFeatures(reference, selection.method = "vst", nfeatures = 3000)

# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

reference <- CellCycleScoring(reference, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)

reference$CC.Difference <- reference$S.Score - reference$G2M.Score
reference <- ScaleData(reference,
                       vars.to.regress = c('Number_UMI','Library','Donor',"CC.Difference"), # as per original article
                       features = rownames(reference))


saveRDS(reference, "seuratObjects/mapRef/regOutPolioudakis.RDS")