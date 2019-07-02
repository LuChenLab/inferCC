#!/usr/bin/env Rscript
# @since 20190701
# Test metacell on bulk data

library(metacell)
library(openxlsx)
library(Seurat)
library(stringr)


setwd("/mnt/raid62/Lung_cancer_10x/")

nc_luad <- read.csv("/mnt/raid62/Lung_cancer_10x/TCGA/LUAD_raw_expresion.csv")
nc_luad[, 1] <- sapply(nc_luad[, 1], function(x) {
    str_split(x, "\\|")[[1]][1]
})

nc_luad = nc_luad[nc_luad$X %in% names(table(nc_luad$X)[table(nc_luad$X) == 1]), ]


rownames(nc_luad) <- nc_luad$X

nc_luad <- nc_luad[, colnames(nc_luad) != "X"]

output_dir = "/mnt/raid62/Lung_cancer_10x/MetaCell/bulk"
setwd(output_dir)



obj <- CreateSeuratObject(
    nc_luad
)

sce <- as.SingleCellExperiment(obj)

################################################
# MetaCell
#
################################################

base_dir = "TCGA_LUAD"
dir.create(base_dir, showWarnings = F, recursive = T)

scdb_init(base_dir, force_reinit=T)

### import expression data from matix
mat = scm_import_sce_to_mat(sce)

object = mat
save(
    object, 
    file=paste(
        output_dir,
        base_dir,
        paste0("mat.", base_dir, ".Rda"),
        sep = "/"
    )
)

# print(dim(mat@mat))

# set image output directory
figs_dir = paste(output_dir, base_dir, "figs", sep = "/")

if(!dir.exists(figs_dir)) dir.create(figs_dir)
scfigs_init(figs_dir)

# basic filtering
mcell_plot_umis_per_cell(base_dir)


nms = c(rownames(mat@mat), rownames(mat@ignore_gmat))
ig_genes = c(grep("^IGJ", nms, v=T), 
             grep("^IGH",nms,v=T),
             grep("^IGK", nms, v=T), 
             grep("^IGL", nms, v=T))

bad_genes = unique(c(grep("^MT-", nms, v=T), grep("^MTMR", nms, v=T), grep("^MTND", nms, v=T),"NEAT1","TMSB4X", "TMSB10", ig_genes))


mcell_mat_ignore_genes(new_mat_id=base_dir, mat_id=base_dir, bad_genes, reverse=F) 

mcell_mat_ignore_small_cells(base_dir, base_dir, 0)

# find feature genes
mcell_add_gene_stat(gstat_id=base_dir, mat_id=base_dir, force=T)


mcell_gset_filter_varmean(gset_id=paste(base_dir, "feats", sep = "_"), gstat_id=base_dir, T_vm=0.08, force_new=T)

mcell_gset_filter_cov(gset_id = paste(base_dir, "feats", sep = "_"), gstat_id=base_dir, T_tot=100, T_top3=2)


mcell_plot_gstats(gstat_id=base_dir, gset_id=paste(base_dir, "feats", sep = "_"))

# build balanced graph
mcell_add_cgraph_from_mat_bknn(mat_id=base_dir, 
                               gset_id = paste(base_dir, "feats", sep = "_"), 
                               graph_id=paste(base_dir, "graph", sep = "_"),
                               K=100,
                               dsamp=T
)

# resampling
mcell_coclust_from_graph_resamp(
    coc_id=paste(base_dir, "coc500", sep = "_"), 
    graph_id=paste(base_dir, "graph", sep = "_"),
    min_mc_size=20, 
    p_resamp=0.75, 
    n_resamp=500
)

# generate co-clustering graph
mcell_mc_from_coclust_balanced(
    coc_id=paste(base_dir, "coc500", sep = "_"), 
    mat_id=base_dir,
    mc_id=paste(base_dir, "mc", sep = "_"), 
    K=30, 
    min_mc_size=30, 
    alpha=2
)


# outliers heatmap
# mcell_plot_outlier_heatmap(mc_id=paste(base_dir, "mc", sep = "_"), mat_id = base_dir, T_lfc=3)


mcell_mc_split_filt(
    new_mc_id=paste(base_dir, "mc_f", sep = "_"), 
    mc_id=paste(base_dir, "mc", sep = "_"), 
    mat_id=base_dir,
    T_lfc=3, 
    plot_mats=F
)


mcell_gset_from_mc_markers(
    gset_id=paste(base_dir, "markers", sep = "_"), 
    mc_id=paste(base_dir, "mc_f", sep = "_")
)


### plot cells by different colors
marks_colors = read.csv("/mnt/raid62/Lung_cancer_10x/MetaCell/20190628_gene_markers.csv")[, 1:5]
marks_colors$group = paste(marks_colors$group, marks_colors$gene)
# 
# marks_colors = data.frame(
#     group=c("Naive", "Prolif", "Tfh", "Treg"),
#     gene=c("TCF7", "TOP2A", "CXCL13", "FOXP3"),
#     color=c("blue", "red", "green", "yellow"),
#     priority=c(1, 1, 1, 1),
#     T_fold=c(1, 1, 1, 1)
# )

# Monocytes	LYZ
# B cells	CD79A
# Mast	CLU
# CD8+	GZMK
# CD8+	CD8A
# CD8+	CD8B
# NK	GNLY
# NK	NKG7
# Treg	CTLA4
# Treg	TNFRSF4


# marks_colors = data.frame(
#     group=c("Monocytes", "B cells", "Mast", "Treg"),
#     gene=c("TCF7", "TOP2A", "CXCL13", "FOXP3"),
#     color=c("blue", "red", "green", "yellow"),
#     priority=c(1, 1, 1, 1),
#     T_fold=c(1, 1, 1, 1)
# )


# marks_colors = read.table(
#     system.file("extdata", "pbmc_mc_colorize.txt", package="metacell"), 
#     sep="\t", 
#     h=T, 
#     stringsAsFactors=F
# )

mc_colorize(paste(base_dir, "mc_f", sep = "_"), marker_colors=marks_colors)


mcell_mc_plot_marks(
    mc_id=paste(base_dir, "mc_f", sep = "_"),
    gset_id=paste(base_dir, "markers", sep = "_"),
    mat_id=base_dir,
    plot_cells=F
)


mcell_mc2d_force_knn(
    mc2d_id=paste(base_dir, "2dproj", sep = "_"),
    mc_id=paste(base_dir, "mc_f", sep = "_"), 
    graph_id=paste(base_dir, "graph", sep = "_")
)



mcell_mc2d_plot(
    mc2d_id = paste(base_dir, "2dproj", sep = "_"), 
    legend_pos="panel",
    plot_edges = TRUE
)

# mcell_mc2d_plot1(
#     mc2d_id = paste(base_dir, "2dproj", sep = "_"), 
#     plot_edges = TRUE
# )


mc_hc = mcell_mc_hclust_confu(
    mc_id=paste(base_dir, "mc_f", sep = "_"), 
    graph_id=paste(base_dir, "graph", sep = "_")
)


mc_sup = mcell_mc_hierarchy(
    mc_id=paste(base_dir, "mc_f", sep = "_"),
    mc_hc=mc_hc, 
    T_gap=0.04
)

mcell_mc_plot_hierarchy(
    mc_id=paste(base_dir, "mc_f", sep = "_"), 
    graph_id=paste(base_dir, "graph", sep = "_"), 
    mc_order=mc_hc$order, 
    sup_mc = mc_sup, 
    width=2800, 
    heigh=2000, 
    min_nmc=2
)

mcell_mc_export_tab(
    mc_id = paste(base_dir, "mc_f", sep = "_"), 
    gstat_id = base_dir, 
    mat_id = base_dir,
    metadata_fields=NULL
)
