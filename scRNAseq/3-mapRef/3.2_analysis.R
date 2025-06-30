# assign colors for cell types
col_vector <- c('PgS'='#e52a25','PgG2M'='#7b170e','vRG'='#486426','oRG'='#213400','panRG-O'='#8bc19b','glycoRG-O'='#01794e','IP'='#bc4419',
                'IPC-Mature-O'='#e6734a','ExDp-O'='#c4e6f9','ExDp1'='#15182f','ExDp2'='#342984','ExNeuNew-O'='#5265ac','ExN'='#275d62','ExM'='#5bf33f',
                'ExM-U'='#103940','ExU-O'='#809bd0','ExPanNeu-O'='#8ccfdd','InhN-O'='#8478b6','InCGE'='#b42585','InMGE'='#612161','OPC'='#bdd028',
                'Astro-O'='#736362','hRG-O'='#90c256','AstroHindb-O'='#c49470','Per'='#fcc800','End'='#4e441e','Unmapped'='#8b8a84')

# pseudobulking
pseudobulk <- data.frame()
for(i in unique(all$mid_unmapped_bucket)){
  if(i=="Unmapped" | i=="Unk") next
  mat <- GetAssayData(object = subset(all, mid_unmapped_bucket==i), slot = "data")
  geneAvgs <- rowMeans(mat)
  pseudobulk <- rbind(pseudobulk, geneAvgs)
}

pseudobulk <- t(pseudobulk)
colnames(pseudobulk) <- unique(all$mid_unmapped_bucket[which(all$mid_unmapped_bucket!="Unmapped",)])
rownames(pseudobulk) <- rownames(all)

# generation of pseudotime trajectory
library(monocle3)

cds <- as.cell_data_set(all)
cds <- estimate_size_factors(cds)
cds <- preprocess_cds(cds, num_dim=100)
cds <- align_cds(cds, alignment_group = "sampleId")
set.seed(1000)
cds <- cluster_cells(cds, resolution=1e-4)

rowData(cds)$gene_short_name <- row.names(rowData(cds))

#trajectory graph
set.seed(123)
cds <- learn_graph(cds, use_partition = TRUE, verbose = FALSE, close_loop=TRUE)

get_earliest_principal_node <- function(cds, time_bin="20"){
  cell_ids <- which(colData(cds)[, "differentiationDay"] == time_bin)
  
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}
cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))