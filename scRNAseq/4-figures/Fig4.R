## Fig.4a
all <- SetIdent(all, value = all$mid_unmapped_bucket)
unmapped <- WhichCells(all, idents = "Unmapped")
DimPlot(all, label=F, cells.highlight=unmapped, cols.highlight = "red", cols= "grey", shuffle=T, sizes.highlight = 0.1, raster=T) + coord_fixed() + theme(legend.position = "none")

## Fig.4b
# pseudobulk matrix based on HQ & LQ cell types
pseudobulk <- data.frame()
for(i in unique(all$mid_unmapped)){
  if(i=="End-LQ" | i=="ExDp2-LQ") next #removing ctypes only present in LQ
  mat <- GetAssayData(object = subset(all, mid_unmapped==i), slot = "data")
  geneAvgs <- rowMeans(mat)
  pseudobulk <- rbind(pseudobulk, geneAvgs)
}

pseudobulk <- t(pseudobulk)
colnames(pseudobulk) <- unique(all$mid_unmapped[which(all$mid_unmapped!="ExDp2-LQ" & all$mid_unmapped!="End-LQ")])
rownames(pseudobulk) <- rownames(all)

# order cols for HQ-LQ comparisons
order <- names(col_vector)[which(names(col_vector) %in% unique(all$mid_gest))]
qual <- rep(c("-HQ","-LQ"),length(order))
order <- paste0(rep(order, each=2), qual)
order <- order[which(order %in% colnames(pseudobulk))]
pseudobulk <- pseudobulk[,order]

fgsea_sets = msigdbr(species = "human", category = "C2", subcategory = "CP:KEGG")
fgsea_sets <- fgsea_sets %>% split(x = .$gene_symbol, f = .$gs_name)
fgsea_sets <- lapply(fgsea_sets, unique)
sel_pathways <- paste0("KEGG_",c("WNT_SIGNALING_PATHWAY","HEDGEHOG_SIGNALING_PATHWAY",
                                 "GLYCOLYSIS_GLUCONEOGENESIS","OXIDATIVE_PHOSPHORYLATION", "CELL_CYCLE",
                                 "RIBOSOME", "NOTCH_SIGNALING_PATHWAY", "LONG_TERM_DEPRESSION", "LONG_TERM_POTENTIATION",
                                 "STEROID_BIOSYNTHESIS", "AXON_GUIDANCE", "GAP_JUNCTION", "MAPK_SIGNALING_PATHWAY",
                                 "SPLICEOSOME", "DNA_REPLICATION", "REGULATION_OF_ACTIN_CYTOSKELETON", "FOCAL_ADHESION", "P53_SIGNALING_PATHWAY"))
fgsea_sets <- fgsea_sets[sel_pathways]
res <- gsva(pseudobulk, fgsea_sets, min.sz = 10, max.sz = 500)

h1 <- pheatmap::pheatmap(res,fontsize = 12, cellwidth=20, cellheight=15, cluster_cols = FALSE, gaps_col=seq(from=2,by=2,to=(ncol(pseudobulk)-2)))

## Fig. 4c
library(RColorBrewer)
unmapped <- subset(all, mid_unmapped_bucket=='Unmapped')

alltp <- sapply(sort(unique(unmapped$differentiationDay)), function(x){
  
  tmp <- unmapped[,unmapped$differentiationDay==x]
  tt <- t(table(tmp$org_unmapped_bucket, tmp$donorId_simplified))
  tt <- as.data.frame(round(tt*100/rowSums(tt),2))
  tt$tp <- x
  colnames(tt) <- c("Donor","CellType","Percentage","diffDay")
  tt
  
}, simplify=F)

alltp <- do.call("rbind", alltp)
alltp$diffDay <- paste0("Day ", alltp$diffDay)

cols <- col_vector[names(col_vector) %in% alltp$CellType]
order <- names(cols)
alltp$CellType <- factor(alltp$CellType, levels = order)

barplotStckd <- ggplot(alltp, aes(fill=CellType, x=diffDay, y=Percentage))+
  #facet_wrap(~tp, scales = "free_x",)+
  geom_bar(position="fill", stat="identity")+
  ggtitle("Cell type composition within unmapped cells")+
  theme_bw()+
  theme(plot.title=element_text(size=8, face="bold"),
        axis.text.x=element_text(size=6,angle=90,hjust=1, vjust=0.5),
        legend.title = element_text(size = 8, face="bold", hjust=0.5), 
        legend.text  = element_text(size = 8),
        legend.key.size = unit(0.5, "lines"))+
  guides(shape = guide_legend(override.aes = list(size = 1)),
         fill = guide_legend(override.aes = list(size = 1), ncol=1))+
  scale_fill_manual(name="Cell types",
                    values = cols)+
  scale_y_continuous(labels = scales::percent)+
  xlab("")+
  ylab("Percentage")

barplotStckd

## Fig.4d
all <- SetIdent(all, value=all$final)
ex <- c("ExDp1","ExN","ExM","ExM-U")
inh <- c("InCGE","InMGE")

# finding DEGs between pan-neuronal cells and fetal-mapped excitatory or inhibitory cell types
panNeuroExMarkers <- FindMarkers(all, ident.1 = "ExPanNeu-O", ident.2 = ex, logfc.threshold = 0, min.pct = 0, min.cells.feature = 0, min.cells.group = 0)
panNeuroInhMarkers <- FindMarkers(all, ident.1 = "ExPanNeu-O", ident.2 = inh, logfc.threshold = 0, min.pct = 0, min.cells.feature = 0, min.cells.group = 0)

fgsea_sets = msigdbr(species = "human", category = "C2", subcategory = "CP:KEGG")
fgsea_sets <- fgsea_sets %>% split(x = .$gene_symbol, f = .$gs_name)
fgsea_sets <- sapply(fgsea_sets, unique)

marker_list <- list(panVsInh=panNeuroInhMarkers, panVsEx=panNeuroExMarkers)
fgseaResTidy_all <- data.frame()
for(i in 1:length(marker_list)){
  markers <- tibble::rownames_to_column(marker_list[[i]], "feature")
  markers %>%
    arrange(desc(avg_log2FC)) %>%
    head(n = 10)
  de_genes <- markers %>%
    arrange(desc(avg_log2FC)) %>% 
    dplyr::select(feature, avg_log2FC)
  ranks <- tibble::deframe(de_genes)
  head(ranks)
  fgseaRes <- fgseaMultilevel(fgsea_sets, stats = ranks, minSize=10, maxSize=500, nPermSimple = 10000, eps=0)
  #tidying data
  fgseaResTidy <- fgseaRes %>%
    as_tibble() %>%
    filter(padj < 0.05) %>% 
    arrange(desc(NES))
  fgseaResTidy %>% 
    dplyr::select(-leadingEdge, -ES) %>% 
    arrange(padj) %>% 
    head()
  fgseaResTidy$pathway <- gsub("KEGG_","",fgseaResTidy$pathway)
  fgseaResTidy$group <- rep(names(marker_list)[i],nrow(fgseaResTidy))
  fgseaResTidy_all <- rbind(fgseaResTidy_all,fgseaResTidy)
}

# module scores of driving genes from top diff. activated pathways

# from pan vs exc
ribo <- fgseaResTidy_all[fgseaResTidy_all$pathway=="RIBOSOME",]$leadingEdge
all <- AddModuleScore(all, features=ribo, assay='RNA', name="pan.ribo")

# from pan vs inh
oxPhos <- fgseaResTidy_all[fgseaResTidy_all$pathway=="OXIDATIVE_PHOSPHORYLATION",]$leadingEdge
all <- AddModuleScore(all, features=oxPhos, assay='RNA', name="pan.OxPhos")
steroid <- fgseaResTidy_all[fgseaResTidy_all$pathway=="STEROID_BIOSYNTHESIS" & fgseaResTidy_all$group=="panVsInh",]$leadingEdge
all <- AddModuleScore(all, features=steroid, assay='RNA', name="pan.Steroid")
ppar <- fgseaResTidy_all[fgseaResTidy_all$pathway=="PPAR_SIGNALING_PATHWAY",]$leadingEdge
all <- AddModuleScore(all, features=ppar, assay='RNA', name="pan.PPAR")
neuroLig <- fgseaResTidy_all[fgseaResTidy_all$pathway=="NEUROACTIVE_LIGAND_RECEPTOR_INTERACTION",]$leadingEdge
all <- AddModuleScore(all, features=neuroLig, assay='RNA', name="pan.neuroLig")

# plotting
feats <- c("pan.OxPhos1","pan.Steroid1","pan.PPAR1","pan.neuroLig1","pan.ribo1")
RidgePlot(subset(all, subset=final %in% c("ExN","ExPanNeu-O","InCGE")),features=feats, ncol=2)

## Fig.4e
# set up technical replicates
all@meta.data[["Samples"]] <- paste0(all$sampleId,"_",all$donorId)
all.sce <- as.SingleCellExperiment(subset(all, subset = differentiationDay == '40')) # donor x protocol representation only at diff day 40
all_milo <- Milo(all.sce)
all_milo <- buildGraph(all_milo, k = 30, d = 15, reduced.dim = "HARMONY")
set.seed(123)
all_milo <- makeNhoods(all_milo, prop = 0.2, k = 30, d=15, refined = TRUE, reduced_dims = "HARMONY")
plotNhoodSizeHist(all_milo) # peak should be between 50-100
all_milo <- countCells(all_milo, meta.data = as.data.frame(colData(all_milo)), sample="Samples")
head(nhoodCounts(all_milo))
all_design <- data.frame(colData(all_milo))[,c("donorId_simplified", "Samples", "Protocol")]
all_design$Samples <- as.factor(all_design$Samples)
all_design$Protocol <- as.factor(all_design$Protocol)
all_design$donorId_simplified <- as.factor(all_design$donorId_simplified)
all_design <- distinct(all_design)
rownames(all_design) <- all_design$Samples
all_design
all_milo <- calcNhoodDistance(all_milo, d=15, reduced.dim = "HARMONY")
da_results <- testNhoods(all_milo, design = ~ donorId_simplified + Protocol, design.df = all_design)
table(da_results$SpatialFDR<0.05)
head(da_results)
da_results %>%
  arrange(SpatialFDR) %>%
  head()
ggplot(da_results, aes(PValue)) + geom_histogram(bins=50)
all_milo <- buildNhoodGraph(all_milo)
umap_pl <- plotReducedDim(all_milo, dimred = "UMAP", colour_by="Protocol", point_size=0.5) +
  guides(fill="none")

# Plot neighbourhood graph to determine colors for each group
nh_graph_pl <- plotNhoodGraphDA(all_milo, da_results, layout="UMAP",alpha=0.1)
library("patchwork")
umap_pl + nh_graph_pl +
  plot_layout(guides="collect") + coord_fixed()

# beeswarm plot
da_results <- annotateNhoods(all_milo, da_results, coldata_col = "final")
ggplot(da_results, aes(final_fraction)) + geom_histogram(bins=50)
da_results$final <- ifelse(da_results$final_fraction < 0.7, "Mixed", da_results$final)
plotDAbeeswarm(da_results, group.by = "final", alpha=0.1) #red (left) is the original protocol