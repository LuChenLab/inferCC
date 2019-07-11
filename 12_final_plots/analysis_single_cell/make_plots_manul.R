#!/usr/bin/env Rscript
# since 2019.07.09

source("/mnt/raid62/Lung_cancer_10x/project_scripts/12_final_plots/analysis_single_cell/basic.R")


### Alveolar II
#### monocle
#### slingshot

setwd("/mnt/raid62/Lung_cancer_10x/final_plots/03_each_cells/Alveolar_II/")

selected_obj <- readRDS("seurat.rds")
temp = str_detect(selected_obj@meta.data$Disease, "Normal")
selected_obj@meta.data$Stage[temp] = selected_obj@meta.data$Disease[temp]

n.pcs = 16
num_cores = 10
low_thresh = 1 ### the threshold of low-expressed
reduced_element = 1 ### the con-founding factors need to be reduced (batch effect)
num_PCs =  16
relative_expr = TRUE ###TRUE when using counts; FALSE when treating PSI
norm_method = 'log' ### 'log' when using counts; 'none' when treating PSI
reduction_method = 'both' ### the method of dimension reduction: tSNE and UMAP, or both
RGE_method = 'DDRTree' ### method to learn cell trajectory (recommoned)
cluster_method = 'densityPeak'


temp = selected_obj@dr$umap@cell.embeddings
### find out which cells needs to kept
choosen = temp$UMAP1 < 10

feature_info = data.frame(row.names=rownames(selected_obj@raw.data[,choosen]), gene_short_name=rownames(selected_obj@raw.data[,choosen]))

### create monocle object
cds <- new_cell_data_set(
    as.matrix(selected_obj@raw.data[,choosen]),
    cell_metadata = selected_obj@meta.data[choosen,],
    gene_metadata = feature_info
)


cds <- preprocess_cds(
    cds, 
    num_dim = n.pcs,
    relative_expr = relative_expr,
    norm_method = 'log', 
    verbse=T,
    residualModelFormulaStr = paste0('~',reduced_element)
)


### replace pca with harmony
# harmony_mat <- 
cds@reducedDims$PCA <- selected_obj@dr$harmony@cell.embeddings[rownames(cds@reducedDims$PCA), 1:n.pcs]

### dimesion reduction
cds <- reduce_dimension(cds, reduction_method = 'UMAP')


### replace umap with Seurat umap results
# reducedDimA(cds)  <- t(selected_obj@dr$umap@cell.embeddings[choosen, ]) -> reducedDims(cds)
print("set umap")
cds@reducedDims$UMAP <- selected_obj@dr$umap@cell.embeddings[rownames(cds@reducedDims$UMAP), ]


print("cluster")
cds <- cluster_cells(cds)
# cds <- partitionCells(cds)
cds <- learn_graph(cds)


# cds <- clusterCells(cds,verbose = F,
#                     method = cluster_method,
#                     res = 1,
#                     louvain_iter = 1,
#                     cores = num_cores)


saveRDS(cds, "monocle3.rds")


pData(cds)[,'Seurat_Cluster'] <- selected_obj@ident[colnames(cds)]
p1 = plot_cells(cds,alpha=0.5, color_cells_by='PatientID')
p2 = plot_cells(cds,alpha=0.5, color_cells_by='batch')

p <- cowplot::plot_grid(p1, p2)
ggsave(filename = "monocle3_umap_paitent_batch.pdf", plot = p, width = 12, height = 6, dpi = 600, units = "in")


### using Benign or Stage I as root
stages = sort(unique(selected_obj@meta.data$Stage[choosen]))

if ("Normal" %in% stages) {
    root = "Normal"
    stages = stages[stages != "Normal"]
} else {
    root = stages[1]
}

get_earliest_principal_node <- function(cds, time_bin="Normal"){
    cell_ids <- which(colData(cds)[, "Stage"] == time_bin)
    
    closest_vertex <-
        cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
    closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
    root_pr_nodes <-
        igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                                  (which.max(table(closest_vertex[cell_ids,]))))]
    
    root_pr_nodes
}

# node_ids = get_correct_root_state(cds,cell_phenotype='Stage', root)
cds <- order_cells(cds, root_pr_nodes = get_earliest_principal_node(cds, "I"))
p1 = plot_cells(
    cds,
    alpha=0.5,
    label_cell_groups=FALSE,
    label_leaves=TRUE,
    label_branch_points=TRUE,
    graph_label_size=1.5
) + labs(title = root[1]) + theme(legend.position = "none")
p2 = plot_cells(
    cds,
    alpha=0.5,
    color_cells_by='Stage',
    label_cell_groups=FALSE,
    label_leaves=TRUE,
    label_branch_points=TRUE,
    graph_label_size=1.5
) + scale_colour_manual(values = stage_colors)

p <- cowplot::plot_grid(p1, p2)

ggsave(filename = "monocle3_trajectory.pdf", plot = p, width = 12, height = 6, dpi = 600, units = "in")


p <- plot_cells(
    cds,
    alpha=0.5,
    color_cells_by='Stage',
    label_cell_groups=FALSE,
    label_leaves=TRUE,
    label_branch_points=TRUE,
    graph_label_size=1.5
) + facet_grid(.~Stage) +
    scale_colour_manual(values = stage_colors) 

ggsave(filename = "monocle3_trajectory_by_stage.pdf", plot = p, width = 6 * length(unique(cds$Stage)), height = 6, dpi = 600, units = "in")

## find pseudotime temporal genes
pr_graph_test <- principalGraphTest(cds, k=3, cores=10)

pr_graph_signif = setDT(rownames_to_column(pr_graph_test, var = 'gene_symbol'))
### choose the significant junctions
pr_graph_signif = pr_graph_signif[morans_I>0.5&qval<0.05][order(qval,-morans_test_statistic)]
cluster_marker_res <- find_cluster_markers(cds, pr_graph_test, group_by = 'Stage', morans_I_threshold = 0.5)

write.xlsx(cluster_marker_res, "monocle3_marker_genes.xlsx")

### need the expressed percentage > 0.2
cluster_marker_res <- setDT(cluster_marker_res)[percentage > 0.25]

top_marker <- cluster_marker_res[order(-specificity,-percentage,-morans_I)][,head(.SD,10),by=Group][order(Group)]
top_marker <- unique(as.character(rev(top_marker$gene_short_name)))

pdf("heatmap_stage_top_marker_monocle.pdf", width = 8, height = 12)
DoHeatmap(object = selected_obj, genes.use = top_marker, group.by = "Stage", slim.col.label = TRUE, remove.key = TRUE, cells.use = colnames(selected_obj@raw.data[, choosen]))
dev.off()

if (length(top_marker) > 1) {
    p <- plot_markers_by_group(cds, rev(top_marker), group_by = 'Stage', ordering_type='none')
    ggsave("monocle3_dotplot_marker_genes.pdf", plot = p, width = 8, height = 12, dpi = 600, units = "in")
    
    
    p <- plot_cell_clusters(cds, markers = top_marker, ncol=5)
    
    width = 15
    if (length(top_marker) < 5) {
        width = 3 * length(top_marker)
    }
    
    ggsave("monocle3_umap_marker_genes.pdf", plot = p, width = width, height = 3 * ceiling(length(top_marker) / 5), dpi = 600, units = "in")
}

