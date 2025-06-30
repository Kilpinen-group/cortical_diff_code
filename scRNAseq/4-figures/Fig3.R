## Fig.3a
library(monocle3)
p1 <- plot_cells(cds,
                 color_cells_by = "pseudotime",
                 label_cell_groups=FALSE,
                 label_leaves=FALSE,
                 label_branch_points=FALSE,
                 label_roots = TRUE,
                 trajectory_graph_color = "black",
                 rasterize = TRUE,
                 graph_label_size=1.5) + coord_fixed()
p1

## Fig.3d
library(GSVA)
# library(clusterProfiler)

fgsea_sets = msigdbr(species = "human", category = "C2", subcategory = "CP:KEGG")
fgsea_sets <- fgsea_sets %>% split(x = .$gene_symbol, f = .$gs_name)
fgsea_sets <- lapply(fgsea_sets, unique)
# kegg_df <- read.gmt("c2.cp.kegg_legacy.v2023.2.Hs.symbols.gmt")
# fgsea_sets <- split(kegg_df$gene, kegg_df$term)
sel_pathways <- paste0("KEGG_",c("WNT_SIGNALING_PATHWAY","HEDGEHOG_SIGNALING_PATHWAY",
                       "GLYCOLYSIS_GLUCONEOGENESIS","OXIDATIVE_PHOSPHORYLATION", "CELL_CYCLE",
                       "RIBOSOME", "NOTCH_SIGNALING_PATHWAY"))
fgsea_sets <- fgsea_sets[sel_pathways]

# pseudobulked object from processing.R
res <- gsva(pseudobulk, fgsea_sets, min.sz = 10, max.sz = 500)
res <- as.data.frame(res)
res$Pathway <- row.names(res)
row.names(res) <- NULL
res_long <- pivot_longer(res,c(1:(ncol(res)-1)),names_to="CellType", values_to="NES")

i <- "KEGG_GLYCOLYSIS_GLUCONEOGENESIS"
j <- "KEGG_OXIDATIVE_PHOSPHORYLATION"
data <- res_long[res_long$Pathway == i | res_long$Pathway==j,]
data$CellType <- factor(data$CellType, levels = data$CellType[order(data[data$Pathway==j,]$NES, decreasing=T)])
p2 <- ggplot(data, aes(x=NES, y=CellType, fill=Pathway)) + geom_bar(width=.7, position="dodge", stat="identity") + theme_bw() + 
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=2,byrow=TRUE)) + scale_fill_discrete(labels = c("Glycolysis", "OxPhos"))#geom_line() ,group=Pathway, color=Pathway
p2

## Fig.3e
fgsea_sets_mito <- msigdbr(species = "human", category = "C5", subcategory = "GO:CC")
fgsea_sets_mito <- fgsea_sets_mito %>% split(x = .$gene_symbol, f = .$gs_name)
fgsea_sets_mito <- lapply(fgsea_sets_mito, unique)
fgsea_sets_mito <- fgsea_sets_mito[grep('MITOCHONDRIAL', names(fgsea_sets_mito))]
res <- gsva(pseudobulk, fgsea_sets_mito, min.sz = 10, max.sz = 500)
h1 <- pheatmap::pheatmap(t(res),fontsize = 4, cellwidth=20, angle_col = "45", cellheight=15, cluster_cols = FALSE)