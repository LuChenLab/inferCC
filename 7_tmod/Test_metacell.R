#!/usr/bin/env Rscript
# since 2019.06.24
# Have nothing to say about this one

library(metacell)

setwd("/mnt/raid62/Lung_cancer_10x/MetaCell/")

base_dir = "2018jz24"

dir.create(base_dir, showWarnings = F, recursive = T)

scdb_init(base_dir, force_reinit=T)

mcell_import_scmat_10x(
    base_dir, 
    "/mnt/raid62/Lung_cancer_10x/CellRanger/lung_cancer/2018jz24/outs/filtered_feature_bc_matrix_for_seurat"
)

mat = scdb_mat(base_dir)

print(dim(mat@mat))

# set image output directory
figs_dir = "/mnt/raid62/Lung_cancer_10x/MetaCell/2018jz24/figs"

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

mcell_mat_ignore_small_cells(base_dir, base_dir, 800)

# find feature genes
mcell_add_gene_stat(gstat_id=base_dir, mat_id=base_dir, force=T)


mcell_gset_filter_varmean(gset_id=paste(base_dir, "feats"), gstat_id=base_dir, T_vm=0.08, force_new=T)

mcell_gset_filter_cov(gset_id = paste(base_dir, "feats"), gstat_id=base_dir, T_tot=100, T_top3=2)


mcell_plot_gstats(gstat_id=base_dir, gset_id=paste(base_dir, "feats"))

# build balanced graph
mcell_add_cgraph_from_mat_bknn(mat_id=base_dir, 
    gset_id = paste(base_dir, "feats"), 
    graph_id=paste(base_dir, "graph"),
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
mcell_plot_outlier_heatmap(mc_id=paste(base_dir, "mc", sep = "_"), mat_id = base_dir, T_lfc=3)


mcell_mc_split_filt(
    new_mc_id=paste(base_dir, "mc_f", sep = "_"), 
    mc_id=paste(base_dir, "mc", sep = "_"), 
    mat_id=base_dir,
    T_lfc=3, 
    plot_mats=T
)


mcell_gset_from_mc_markers(
    gset_id=paste(base_dir, "markers", sep = "_"), 
    mc_id=paste(base_dir, "mc_f", sep = "_")
)


marks_colors = read.table(
    system.file("extdata", "pbmc_mc_colorize.txt", package="metacell"), 
    sep="\t", 
    h=T, 
    stringsAsFactors=F
)

mc_colorize(paste(base_dir, "mc_f", sep = "_"), marker_colors=marks_colors)


mcell_mc_plot_marks(
    mc_id=paste(base_dir, "mc_f", sep = "_"), 
    gset_id=paste(base_dir, "markers", sep = "_"),
     mat_id=base_dir
)


mcell_mc2d_force_knn(
    mc2d_id=paste(base_dir, "2dproj", sep = "_"),
    mc_id=paste(base_dir, "graph", sep = "_"), 
    graph_id=paste(base_dir, "graph", sep = "_")
)


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


library(Seurat)
library(stringr)

expression_matrix <- Read10X(data.dir = "/mnt/raid62/Lung_cancer_10x/CellRanger/lung_cancer/2018jz24/outs/filtered_feature_bc_matrix_for_seurat")
obj = CreateSeuratObject(raw.data = expression_matrix)

mito.genes <- grep(pattern = "^MT-", x = rownames(x = obj@raw.data), value = TRUE)
percent.mito <- Matrix::colSums(obj@raw.data[mito.genes, ])/Matrix::colSums(obj@raw.data)

obj <- AddMetaData(object = obj, metadata = percent.mito, col.name = "percent.mito")

obj <- NormalizeData(object = obj, normalization.method = "LogNormalize", scale.factor = 10000)
obj <- FindVariableGenes(object = obj, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
obj <- ScaleData(object = obj, vars.to.regress = c("nUMI", "percent.mito"))
obj <- RunPCA(object = obj, pc.genes = obj@var.genes, do.print = FALSE, pcs.compute = 100)


n.pcs = 15

obj <- RunTSNE(object = obj, dims.use = 1:n.pcs, do.fast = TRUE, reduction.use = "pca", reduction.name = "tsne")

obj <- FindClusters(
    object = obj, 
    reduction.type = "pca", 
    dims.use = 1:n.pcs, 
    resolution = 0.6, 
    print.output = 0, 
    save.SNN = FALSE, 
    force.recalc = TRUE
)


load(paste(base_dir, "mc.2018jz24_mc_f.Rda", sep = "/"))

p1 <- DimPlot(
    obj, 
    reduction.use = "tsne", 
    group.by = "res.0.6", 
    no.legend = TRUE, 
    do.label = TRUE, 
    label.size = 5, 
    do.return = TRUE, 
    pt.size = 0.01
)

p2 <- DimPlot(
    obj,
    cells.use = sapply(names(object@mc), function(x) {
        return(str_split(x, "-")[[1]][1])
    }), 
    reduction.use = "tsne", 
    group.by = "res.0.6", 
    no.legend = TRUE, 
    do.label = TRUE, 
    label.size = 5, 
    do.return = TRUE, 
    pt.size = 0.01
)

p <- cowplot::plot_grid(p1, p2, ncol = 2)

ggsave(
    "tsne_cluster.png", 
    plot=p, 
    width=12, 
    height=6, 
    dpi=600, 
    units="in"
)