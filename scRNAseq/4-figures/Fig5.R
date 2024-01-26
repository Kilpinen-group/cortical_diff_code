## Fig.5c
library(GSVA)
all$exptPop <- paste0(all$donorId_simplified,"_d",all$differentiationDay)

#pseudobulk the data
pseudobulk <- data.frame()
for(i in unique(all$exptPop)){
  mat <- GetAssayData(object = subset(all, exptPop==i), slot = "data")
  geneAvgs <- rowMeans(mat)
  pseudobulk <- rbind(pseudobulk, geneAvgs)
}

pseudobulk <- t(pseudobulk)
colnames(pseudobulk) <- unique(all$exptPop)
rownames(pseudobulk) <- rownames(all)

fgsea_sets = msigdbr(species = "human", category = "C2", subcategory = "CP:KEGG")
fgsea_sets <- fgsea_sets %>% split(x = .$gene_symbol, f = .$gs_name)
fgsea_sets <- lapply(fgsea_sets, unique)
sel_pathways <- paste0("KEGG_",c("WNT_SIGNALING_PATHWAY","HEDGEHOG_SIGNALING_PATHWAY",
                                 "GLYCOLYSIS_GLUCONEOGENESIS","OXIDATIVE_PHOSPHORYLATION", "CELL_CYCLE",
                                 "RIBOSOME", "NOTCH_SIGNALING_PATHWAY"))
fgsea_sets <- fgsea_sets[sel_pathways]

res <- gsva(pseudobulk, fgsea_sets, min.sz = 10, max.sz = 500)

tt <- table(all$exptPop,all$final)
tt <- round(tt*100/rowSums(tt),2)
tt <- tt[,c("InCGE","InMGE","InhN-O")]
tt <- rowSums(tt)
tt <- tt[colnames(pseudobulk)]
library(ComplexHeatmap)
library(circlize)
row_ha = rowAnnotation("% inhibitory neurons" = anno_barplot(tt), annotation_name_rot=90)
h1 <- Heatmap(t(res),  
              name = "GSVA enrichment score",
              right_annotation = row_ha,
              cluster_columns = TRUE,
              show_column_dend = TRUE,
              cluster_column_slices = TRUE,
              column_names_gp = gpar(fontsize = 6),
              column_gap = unit(0.5, "mm"),
              cluster_rows = TRUE,
              show_row_dend = FALSE,
              show_row_names = TRUE,
              row_names_gp = gpar(fontsize = 8),
              use_raster = TRUE,
              show_column_names = TRUE)
h1

# determine segments from trajectory plot
plot_cells(cds,
           color_cells_by = "mid_unmapped_bucket",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_roots = FALSE,
           label_branch_points=FALSE,
           label_principal_points = TRUE, #TRUE allows visualization to choose_graph_segments()
           trajectory_graph_color = "black",
           label_cell_groups=FALSE,
           cell_size = 0.35,
           rasterize = TRUE,
           group_label_size = 6) + scale_color_manual(values=col_vector[unique(all$mid_unmapped_bucket)]) + coord_fixed()

# IP-excitatory segment - Y_448-Y_27, IP-inhibitory - Y_324-Y_198 (plus IP-IP)
segment1 <- choose_graph_segments(cds,
                                  reduction_method = "UMAP",
                                  starting_pr_node = "Y_27",
                                  ending_pr_nodes = "Y_448",
                                  return_list = TRUE,
                                  clear_cds = TRUE)
segment2 <- choose_graph_segments(cds,
                                  reduction_method = "UMAP",
                                  starting_pr_node = "Y_198",
                                  ending_pr_nodes = "Y_324",
                                  return_list = TRUE,
                                  clear_cds = TRUE)
segment3 <- choose_graph_segments(cds,
                                  reduction_method = "UMAP",
                                  starting_pr_node = "Y_324",
                                  ending_pr_nodes = "Y_448",
                                  return_list = TRUE,
                                  clear_cds = TRUE)
segment <- cds[,unique(c(segment1$cells,segment2$cells,segment3$cells))]
rowData(segment)$gene_short_name <- row.names(rowData(segment))

# detect genes which vary over the trajectory within selected segment
set.seed(123)
subset_pr_test_res <- graph_test(segment, neighbor_graph="principal_graph", cores=4)
subset_pr_test_res <- subset_pr_test_res[order(subset_pr_test_res$morans_I, decreasing = TRUE),]
pr_deg_ids <- row.names(subset(subset_pr_test_res, q_value < 0.05))
head(pr_deg_ids)

# building modules of genes differentially expressed across the trajectory
subset_module_df <- find_gene_modules(segment[pr_deg_ids,], resolution=0.001, random_seed=123)
head(subset_module_df)

agg_mat <- aggregate_gene_expression(segment, subset_module_df) #creates an aggregate expression for each identified module per cell in the subset
module_dendro <- hclust(dist(agg_mat)) #hierarchical clustering on a set of dissimilarities (distance matrix between rows of agg_mat)
subset_module_df$module <- factor(subset_module_df$module, levels = row.names(agg_mat)[module_dendro$order]) #order the plot so we can see which modules come on first

# visualize module specificity within segment
plot_cells(segment,
           genes=subset_module_df,
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)

# module 12 is specific to the inhibitory trajectory - find top genes
genes <- subset_module_df[subset_module_df$module==12,]$id
ord <- pr_deg_ids[which(pr_deg_ids %in% genes)] # ordered by moran's I value
genes <- genes[match(ord, genes)]

## Fig.5d
library(ggpubr)
data <- as.data.frame(t(GetAssayData(subset(all, mid_unmapped_bucket %in% c("PgS","PgG2M","oRG","vRG")  &
                                              differentiationDay=='40'), slot='data')[genes[c(1,2,3,5)],]))
data$donor <- unname(all$donorId_simplified[rownames(data)])
data$rep <- unname(all$sampleId[rownames(data)])
data <- pivot_longer(data,c(1:4),names_to = "Gene")
head(data)
ggplot(data,aes(x=donor,y=value, fill=rep))+geom_boxplot()+
  ggtitle("Relative expression of interneuron-associated TFs in d40 progenitors")+theme_bw()+facet_wrap(~Gene)

## Fig. 5e
library(scattermore)
plot_list <- list()
for(i in c('11.4','61.2')){ #unique(all$donorId_simplified) to test all donors
  cds_subset <- cds[
    c(genes[c(5,2)],"RPL19","RPL39"),
    colData(cds) %>%
      subset(
        donorId_simplified == i
      ) %>%
      row.names]
  plot_list[[i]] <- plot_genes_in_pseudotime(cds_subset,
                                             color_cells_by="mid_unmapped_bucket",
                                             nrow = 4, ncol = 1,
                                             horizontal_jitter = 0.5,
                                             min_expr=0.1) + ggtitle(i) + scale_color_manual(values=col_vector[unique(all$mid_unmapped_bucket)]) + geom_scattermore()
}
combined_plot <- Reduce(`+`, plot_list[!sapply(plot_list, is.null)])
combined_plot