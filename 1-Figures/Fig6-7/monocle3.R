library(monocle3)
library(destiny)

setwd("LungCancer10x/")


stage_colors = c(
    "I"="#65A9A3", 
    "II"="#4A933E", 
    "III"="#EC7A21", 
    "IV"="#D73F47", 
    "LUAD_Normal"="#FECC1B", 
    "LUSC_Normal"="#778793"
)


### function to convert seurat obj to monocle3 obj
convert_seurat_monocle <- function(
    obj, 
    choosen=NULL, 
    n.pcs = 10, 
    norm_method = 'log', 
    relative_expr = TRUE,
    reduced_element = 1
) {
    
    if (is.null(choosen)) {
        choosen = colnames(obj@raw.data)
    }
    
    feature_info = data.frame(
        row.names=rownames(obj@raw.data[,choosen]), 
        gene_short_name=rownames(obj@raw.data[,choosen])
    )
    
    ### create monocle object
    cds <- new_cell_data_set(
        as.matrix(obj@raw.data[,choosen]),
        cell_metadata = obj@meta.data[choosen,],
        gene_metadata = feature_info
    )
    
    print("normalize")
    # cds <- estimateSizeFactors(cds)
    # cds <- estimateDispersions(cds)
    
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
    cds@reducedDims$PCA <- obj@dr$harmony@cell.embeddings[rownames(cds@reducedDims$PCA), 1:n.pcs]
    
    ### dimesion reduction
    cds <- reduce_dimension(cds, reduction_method = 'UMAP')
    
    
    ### replace umap with Seurat umap results
    # reducedDimA(cds)  <- t(obj@dr$umap@cell.embeddings[choosen, ]) -> reducedDims(cds)
    print("set umap")
    cds@reducedDims$UMAP <- obj@dr$umap@cell.embeddings[rownames(cds@reducedDims$UMAP), ]
    
    
    print("cluster")
    cds <- cluster_cells(cds)
    # cds <- partitionCells(cds)
    cds <- learn_graph(cds)
    
    
    # cds <- clusterCells(cds,verbose = F,
    #                     method = cluster_method,
    #                     res = 1,
    #                     louvain_iter = 1,
    #                     cores = num_cores)
    
    
    cds
}



obj <- readRDS("03_each_cells/ATII_Basal/Basal.rds")


### Diffusion map
make_dc_plot <- function(dpt, obj, dims = 6, output = "Basal") {

    mtx = as.data.frame(dpt@dm@eigenvectors)
    rownames(mtx) <- rownames(dpt@dm@data_env$data)
    mtx$Stage <- obj@meta.data[rownames(mtx), "Stage"]
    
    for(i in seq(1, dims, 2)) {
        temp <- mtx[, c(paste0("DC", i), paste0("DC", i + 1), "Stage")]
        colnames(temp) <- c("DC1", "DC2", "Stage")
        
        p <- ggplot(temp, aes(x=DC1, y = DC2, color = Stage)) +
            geom_point() +
            scale_color_manual(values = stage_colors) +
            theme_bw() +
            theme(
                aspect.ratio = 1,
                legend.position = c(0.6, 0.05),
                axis.text = element_text(size = 15),
                axis.title = element_text(size = 18),
                legend.text = element_text(size = 15),
                legend.background = element_blank()
            ) +
            guides(color = guide_legend(override.aes = list(size = 5), nrow = 1)) +
            labs(color = "", x = paste0("DC", i), y = paste0("DC", i + 1))
        
        ggsave(
            filename = paste0("03_each_cells/ATII_Basal/", output, "_DM_", i, ".pdf"),
            plot = p,
            width = 6,
            height = 6
        )
    }
}


dmp <- DiffusionMap(obj@dr$pca@cell.embeddings[, 1:10])
dpt <- DPT(dmp)

saveRDS(dpt, "03_each_cells/ATII_Basal/Basal_DPT.rds")
dpt <- readRDS("03_each_cells/ATII_Basal/Basal_DPT.rds")

make_dc_plot(dpt, obj)


mtx = as.data.frame(dpt@dm@eigenvectors)
rownames(mtx) <- rownames(dpt@dm@data_env$data)
mtx$ident <- obj@meta.data[rownames(mtx), "ident"]


p <- ggplot(mtx, aes(x=DC1, y = DC2, color = as.factor(ident))) +
    geom_point() +
    theme_bw() +
    theme(
        aspect.ratio = 1,
        legend.position = c(0.6, 0.1),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 15),
        legend.background = element_blank()
    ) +
    guides(color = guide_legend(override.aes = list(size = 5), nrow = 1)) +
    labs(color = "")


p

ggsave(
    filename = "03_each_cells/ATII_Basal/Basal_DM.pdf",
    plot = p,
    width = 6,
    height = 6
)


### Monocle3
mono_obj <- convert_seurat_monocle(obj)


temp = obj@dr$umap@cell.embeddings
### find out which cells needs to kept
choosen = rep(TRUE, nrow(temp))
for(i in 1:ncol(temp)) {
    print(sum(choosen))
    summ = summary(temp[, i])
    print(summ)
    if (summ[1] / summ[2] > 20) {
        choosen = choosen & temp[, i] >= summ[2]
    }
    
    if (summ[6] / summ[5] > 20) {
        choosen = choosen & temp[, i] <= summ[5]
    }
}


# pData(cds)[,'Seurat_Cluster'] <- obj@ident[colnames(cds)]
# p1 = plot_cells(cds,alpha=0.5, color_cells_by='PatientID')
# p2 = plot_cells(cds,alpha=0.5, color_cells_by='batch')

# p <- cowplot::plot_grid(p1, p2)
# ggsave(filename = "monocle3_umap_paitent_batch.pdf", plot = p, width = 12, height = 6, dpi = 600, units = "in")


### using Benign or Stage I as root
stages = sort(unique(obj@meta.data$Stage))

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
mono_obj <- order_cells(mono_obj, root_pr_nodes = get_earliest_principal_node(mono_obj, "I"))
p1 = plot_cells(
    mono_obj,
    alpha=1,
    label_cell_groups=FALSE,
    label_leaves=TRUE,
    label_branch_points=TRUE,
    graph_label_size=3
) + labs(title = root[1]) + theme(legend.position = "none")
p2 = plot_cells(
    mono_obj,
    alpha=1,
    color_cells_by='Stage',
    label_cell_groups=FALSE,
    label_leaves=TRUE,
    label_branch_points=TRUE,
    graph_label_size=3
) + scale_colour_manual(values = stage_colors)

p <- cowplot::plot_grid(p1, p2)
p
# ggsave(filename = "monocle3_trajectory.pdf", plot = p, width = 12, height = 6, dpi = 600, units = "in")


p <- plot_cells(
    mono_obj,
    alpha=0.5,
    color_cells_by='Stage',
    label_cell_groups=FALSE,
    label_leaves=TRUE,
    label_branch_points=TRUE,
    graph_label_size=1.5
) + facet_grid(.~Stage) +
    scale_colour_manual(values = stage_colors) 

# ggsave(filename = "monocle3_trajectory_by_stage.pdf", plot = p, width = 6 * length(unique(cds$Stage)), height = 6, dpi = 600, units = "in")


####

obj <- readRDS("03_each_cells/ATII_Basal/ATII.rds")

dmp <- DiffusionMap(obj@dr$pca@cell.embeddings[, 1:10])
dpt <- DPT(dmp)

saveRDS(dpt, "03_each_cells/ATII_Basal/Basal_DPT.rds")

make_dc_plot(dpt, obj, output = "ATII")