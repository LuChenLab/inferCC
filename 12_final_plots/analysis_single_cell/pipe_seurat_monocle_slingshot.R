dpi = 600

# load packages
args = commandArgs(trailingOnly=TRUE)

# args = c(
#     '/mnt/raid62/Lung_cancer_10x/project_scripts/12_final_plots/analysis_single_cell/basic.R',
#     '/mnt/raid62/Lung_cancer_10x/final_plots/02_rds/01_split_cells/LUAD/Alveolar_II.rds',
#     '/mnt/raid62/Lung_cancer_10x/final_plots/03_each_cells/LUAD/Alveolar_II',
#     'DO'
# )

basic_script = args[1]
input_rds = args[2]
image_dir = args[3]
trajectory = args[4]


print(args)

source(basic_script)

set.seed(1)

dir.create(image_dir, showWarnings = FALSE, recursive=T)
img.units = 6

setwd(image_dir)


################

### 1. Load raw counts and meta info

if (file.exists("seurat.rds")) {
    selected_obj <- readRDS("seurat.rds")

    kee = import("kneed")
    n.pcs = kee$KneeLocator(
        1:ncol(selected_obj@dr$pca@cell.embeddings),
        apply(selected_obj@dr$pca@cell.embeddings, 2, sd), 
        curve='convex', 
        direction='decreasing'
    )
    n.pcs = n.pcs$knee
} else {
    obj <- readRDS(input_rds)
    print(obj)
    ### 2. check is the variables like Stage and Patient are highly variable
    if (trajectory == "DO") {

        meta = obj@meta.data
        meta = as.data.frame(meta %>% select(PatientID, Gender, Stage, Disease, Age))

        is_variable <- function(meta, targets = "Stage", p.value=0.05) {
            return (t.test(table(meta[, i]))$p.value > p.value)
        }

        m <- get_smote_match(meta)

        ### 3. create selected object and make barplots to check the SMOTE effect
        selected_obj <- CreateSeuratObject(
            obj@raw.data[, rownames(m)],
            meta.data = obj@meta.data[rownames(m),]
        )
        print(selected_obj)
        stats = as.data.frame(table(obj@meta.data$Stage))
        stats$label = "Before"
        stats$type = "Stage"

        temp = as.data.frame(table(selected_obj@meta.data$Stage))
        temp$label = "After"
        temp$type = "Stage"
        stats = rbind(stats, temp)

        temp = as.data.frame(table(obj@meta.data$PatientID))
        temp$label = "Before"
        temp$type = "Patient"
        stats = rbind(stats, temp)

        temp = as.data.frame(table(selected_obj@meta.data$PatientID))
        temp$label = "After"
        temp$type = "Patient"
        stats = rbind(stats, temp)


        p <- ggplot(stats, aes(x=Var1, y=log10(Freq), fill=label)) + 
            geom_bar(stat = "identity", position = "dodge") +
            facet_grid(.~type, scales = "free", space = "free") +
            theme(axis.text.x=element_text(angle = 45))

        ggsave("stats_SMOTE.pdf", plot = p, width = 12, height = 6, dpi=dpi, units = "in")
    } else {
        selected_obj = obj
    }


    ### 4. basic 

    #### 1. vlnplot
    mito.genes <- grep(pattern = "^MT-", x = rownames(x = selected_obj@raw.data), value = TRUE)
    percent.mito <- Matrix::colSums(selected_obj@raw.data[mito.genes, ])/Matrix::colSums(selected_obj@raw.data)

    # AddMetaData adds columns to object@meta.data, and is a great place to
    # stash QC stats
    selected_obj <- AddMetaData(object =selected_obj, metadata = percent.mito, col.name = "percent.mito")

    pdf("VlnPlot.pdf", width = 18, height = 6)
    VlnPlot(object = selected_obj, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
    dev.off()


    #### 2. GenePlot

    pdf("GenePlot.pdf", width = 12, height = 6)
    par(mfrow = c(1, 2))
    GenePlot(object = selected_obj, gene1 = "nUMI", gene2 = "percent.mito")
    GenePlot(object = selected_obj, gene1 = "nUMI", gene2 = "nGene")
    dev.off()


    #### 3. do PCA and harmony

    pcs = min(100, min(nrow(selected_obj@raw.data), ncol(selected_obj@raw.data)) - 1)

    selected_obj <- NormalizeData(object = selected_obj, normalization.method = "LogNormalize", scale.factor = 10000)
    selected_obj <- FindVariableGenes(object = selected_obj, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
    selected_obj <- ScaleData(object = selected_obj, vars.to.regress = c("nUMI", "percent.mito"))
    selected_obj <- RunPCA(object = selected_obj, pc.genes = selected_obj@var.genes, do.print = FALSE, pcs.compute = pcs)

    selected_obj <- RunHarmony(selected_obj, c("batch", "PatientID"))


    #### 4. make dotplot of SD of PCAs and Harmonys
    pdf("PCElbowPlot.pdf", width = 12, height = 6)
    PCElbowPlot(selected_obj, num.pc=pcs)  +
        scale_x_continuous(breaks=1:pcs) + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
    dev.off()


    temp <- as.data.frame(apply(selected_obj@dr$harmony@cell.embeddings, 2, sd))
    colnames(temp) <- c("sd")
    temp$Harmony = 1:nrow(temp)

    kee = import("kneed")
    n.pcs = kee$KneeLocator(
        1:ncol(selected_obj@dr$pca@cell.embeddings),
        apply(selected_obj@dr$pca@cell.embeddings, 2, sd), 
        curve='convex', 
        direction='decreasing'
    )
    n.pcs = n.pcs$knee

    p <- ggplot(temp, aes(x=Harmony, y=sd)) + 
        geom_point() + 
        ylab("Standard Dei of Harmony")  +
        scale_x_continuous(breaks=1:100) + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
        geom_vline(
            xintercept = n.pcs, 
            linetype="dashed", 
            color = "grey"
        )
    ggsave("HarmonyPlot.pdf", plot = p, width = 12, height = 6, dpi = dpi, units = "in")


    #### 5. tSNE and UMAP
    selected_obj <- RunTSNE(
        object = selected_obj, 
        dims.use = 1:n.pcs, 
        do.fast = TRUE, 
        reduction.use = "harmony", 
        reduction.name = "tsne",
        perplexity=min(30, floor(pcs / 3) + 1)
    )
    selected_obj <- RunUMAP(
        object = selected_obj, 
        dims.use = 1:n.pcs, 
        do.fast = TRUE, 
        reduction.use = "harmony",
        n_neighbors = min(30, pcs)
    )
    selected_obj <- FindClusters(object = selected_obj, reduction.type = "harmony", dims.use = 1:n.pcs, resolution = 0.8, print.output = 0, save.SNN = FALSE, force.recalc = TRUE)

    #### 6. change cluster id from 0 to 1, and save seurat object
    selected_obj@meta.data$res.0.8 <- as.numeric(selected_obj@meta.data$res.0.8) + 1
    saveRDS(selected_obj, "seurat.rds")

    #### 7. Dump seurat object data to several files for scanpy
    dump_seurat_for_scanpy(selected_obj, "scanpy")


    ### check cluster and batch effect by tsne and umap

    p <- myDimPlot(
        selected_obj, 
        reduction.use = "tsne", 
        group.by = "res.0.8", 
        no.legend = TRUE, 
        do.label = TRUE, 
        label.size = 5, 
        do.return = TRUE, 
        pt.size = 0.01
    )
    ggsave("tsne_cluster.pdf", plot=p, width=6, height=6, dpi=dpi, units="in")


    p <- myDimPlot(
        selected_obj, 
        reduction.use = "tsne", 
        group.by = "batch", 
        no.legend = FALSE, 
        do.label = FALSE, 
        do.return = TRUE, 
        pt.size = 0.01
    ) + theme(legend.position = c(0.1, 0.1))
    ggsave("tsne_batch.pdf", plot=p, width=6, height=6, dpi=dpi, units="in")


    p <- myDimPlot(
        selected_obj, 
        reduction.use = "umap", 
        group.by = "res.0.8", 
        no.legend = TRUE, 
        do.label = TRUE, 
        label.size = 5, 
        do.return = TRUE, 
        pt.size = 0.01
    )
    ggsave("umap_cluster.pdf", plot=p, width=6, height=6, dpi=dpi, units="in")

    p <- myDimPlot(
        selected_obj, 
        reduction.use = "umap", 
        group.by = "batch", 
        no.legend = FALSE, 
        do.label = FALSE, 
        do.return = TRUE, 
        pt.size = 0.01
    ) + theme(legend.position = c(0.1, 0.1))
    ggsave("umap_batch.pdf", plot=p, width=6, height=6, dpi=dpi, units="in")


    ### 

    # dir.create("tsne_by_disease", showWarnings=F)
    # batch = sort(unique(selected_obj@meta.data$Disease))
    # for (i in 1:length(batch)){
    #     print(i)
        
    #     stage = selected_obj@meta.data[selected_obj@meta.data$Disease == batch[i], ]
    #     p <- myDimPlot(
    #         selected_obj, 
    #         reduction.use = "tsne", 
    #         cells.highlight = rownames(stage), 
    #         do.return=TRUE,
    #         no.legend = FALSE,
    #         pt.size = 0.5,
    #         cols.highlight = rep(disease_colors[batch[i]], nrow(stage))
    #     )
        
    #     temp = selected_obj@dr$tsne@cell.embeddings
    #     p = p + geom_text(
    #         x = min(temp[, 1]) + (max(temp[, 1]) - min(temp[, 1])) * 0.01, 
    #         y = min(temp[, 2]) + (max(temp[, 2]) - min(temp[, 2])) * 0.01, 
    #         label = batch[i], 
    #         hjust = 0, 
    #         vjust = 0
    #     )
        
    #     ggsave(paste0("tsne_by_disease/", batch[i], ".pdf"), plot = p, width = 6, height = 6, units = "in", dpi = dpi)
    # }


    # dir.create("tsne_by_stage", showWarnings=F)
    # batch = sort(unique(selected_obj@meta.data$Stage))
    # for (i in 1:length(batch)){
    #     print(i)
        
    #     stage = selected_obj@meta.data[selected_obj@meta.data$Stage == batch[i], ]
    #     p <- myDimPlot(
    #         selected_obj, 
    #         reduction.use = "tsne", 
    #         cells.highlight = rownames(stage), 
    #         do.return=TRUE,
    #         no.legend = FALSE,
    #         pt.size = 0.5,
    #         cols.highlight = rep(stage_colors[batch[i]], nrow(stage))
    #     )
        
    #     temp = selected_obj@dr$tsne@cell.embeddings
    #     p = p + geom_text(
    #         x = min(temp[, 1]) + (max(temp[, 1]) - min(temp[, 1])) * 0.01, 
    #         y = min(temp[, 2]) + (max(temp[, 2]) - min(temp[, 2])) * 0.01, 
    #         label = batch[i], 
    #         hjust = 0, 
    #         vjust = 0
    #     )
        
    #     ggsave(paste0("tsne_by_stage/", batch[i], ".pdf"), plot = p, width = 6, height = 6, units = "in", dpi = dpi)
    # }


    # #### 3. tsne by patient

    # patients = sort(unique(selected_obj@meta.data$PatientID))

    # p <- make_tsne_plot(selected_obj, color="PatientID", colors = patient_colors, pt.size = 0.5, guide_ncol = 3) + 
    #     theme(
    #         legend.position = c(0.85, 0.85),
    #         legend.title = element_blank(),
    #         legend.background=element_blank(),
    #         legend.key = element_rect(colour = NA, fill = NA)
    #     )
    # ggsave("tsne_by_patient.pdf", plot = p, width = 6, height = 6, units = "in", dpi = dpi)


    # dir.create("tsne_by_patient", showWarnings=F)
    # for (i in 1:length(patients)){
    #     print(i)
        
    #     stage = selected_obj@meta.data[selected_obj@meta.data$PatientID == patients[i], ]
    #     p <- myDimPlot(
    #         selected_obj, 
    #         reduction.use = "tsne", 
    #         cells.highlight = rownames(stage), 
    #         do.return=TRUE,
    #         no.legend = FALSE,
    #         pt.size = 0.5,
    #         cols.highlight = rep(patient_colors[i], nrow(stage))
    #     )
        
    #     temp = selected_obj@dr$tsne@cell.embeddings
    #     p = p + geom_text(
    #         x = min(temp[, 1]) + (max(temp[, 1]) - min(temp[, 1])) * 0.01, 
    #         y = min(temp[, 2]) + (max(temp[, 2]) - min(temp[, 2])) * 0.01, 
    #         label = patients[i], 
    #         hjust = 0, 
    #         vjust = 0
    #     )
        
    #     ggsave(paste0("tsne_by_patient/",as.character(patients[i]), ".pdf"), plot = p, width = 6, height = 6, units = "in", dpi = dpi)
    # }
}

temp = str_detect(selected_obj@meta.data$Disease, "Normal")
selected_obj@meta.data$Stage[temp] = selected_obj@meta.data$Disease[temp]

#### 1. tsne by disease
p <- make_tsne_plot(selected_obj, pt.size = 0.5, colors = disease_colors, guide_ncol=4) + 
    theme(
        legend.position = "top",
        legend.title = element_blank(),
        legend.background=element_blank(),
        legend.justification = "center"
    ) + 
    guides(color = guide_legend(nrow = 1), override.aes = list(size = 3))


ggsave("tsne_by_disease.pdf", plot = p, width = 6, height = 6, dpi = dpi, units = "in")


#### 2. tsne by stage
p <- make_tsne_plot(selected_obj, color="Stage", colors = stage_colors, pt.size = 0.5, guide_ncol=4) + 
    theme(
        legend.position = "top",
        legend.title = element_blank(),
        legend.background=element_blank(),
        legend.justification = "center"
    ) +
    guides(color = guide_legend(nrow = 1), override.aes = list(size = 3))
ggsave("tsne_by_stage.pdf", plot = p, width = 6, height = 6, units = "in", dpi = dpi)


#### 4. make combine barplot
p <- make_combined_barplot(selected_obj, rel_height = c(0.3, 1), with_disease = length(unique(selected_obj@meta.data$Disease)) > 1)
ggsave("barplots.pdf", plot = p, width = 15, height = 10, dpi = dpi, units = "in", limitsize = FALSE)


### Find markers
if(file.exists("markers_by_stage.xlsx")) {
    markers_by_stage = read.xlsx("markers_by_stage.xlsx", sheet = 1, rowNames=T)
} else {
    markers_by_stage = find_markers(selected_obj, group.by = "Stage")

    wb = createWorkbook()
    addWorksheet(wb, "marker genes")
    writeData(wb, 1, markers_by_stage, rowNames = TRUE)
    saveWorkbook(wb, file = "markers_by_stage.xlsx", overwrite = T)
}

if(file.exists("markers_by_cluster.xlsx")) {
    markers_by_cluster = read.xlsx("markers_by_cluster.xlsx", sheet = 1, rowNames=T)
} else {
    markers_by_cluster = find_markers(selected_obj, group.by = "res.0.8")

    wb = createWorkbook()
    addWorksheet(wb, "marker genes")
    writeData(wb, 1, markers_by_cluster, rowNames = TRUE)
    saveWorkbook(wb, file = "markers_by_cluster.xlsx", overwrite = T)
}


#### 5. make dotplot and heatmap of top10 markers
# if (ncol(selected_obj@raw.data) > 30000){
#     cells = sample(1:ncol(obj@raw.data), 30000)
#     cells = colnames(selected_obj@raw.data)[cells]
# } else {
#     cells = colnames(selected_obj@raw.data)
# }

# top10 <- markers_by_cluster %>% dplyr::group_by(ident) %>% top_n(n = 10, wt = avg_logFC)

# pdf("heatmap_cluster.pdf", width = 8, height = 12)
# DoHeatmap(object = selected_obj, genes.use = top10$gene, slim.col.label = TRUE, remove.key = TRUE, cells.use = cells)
# dev.off()

# p <- DotPlot(selected_obj, genes.plot = unique(top10$gene), do.return = TRUE, x.lab.rot = TRUE) + coord_flip() + theme(axis.text.x = element_text(angle = 0))
# ggsave("dotplot_top10_cluster.pdf", p, dpi = dpi, width = 6, height = 12, units = "in")


# top10 <- markers_by_stage %>% dplyr::group_by(ident) %>% top_n(n = 10, wt = avg_logFC)

# pdf("heatmap_stage.pdf", width = 8, height = 12)
# DoHeatmap(object = selected_obj, genes.use = top10$gene, group.by = "Stage", slim.col.label = TRUE, remove.key = TRUE, cells.use = cells)
# dev.off()

# p <- DotPlot(selected_obj, genes.plot = unique(top10$gene), do.return = TRUE, x.lab.rot = TRUE, group.by = "Stage") + coord_flip() + theme(axis.text.x = element_text(angle = 0))
# ggsave("dotplot_top10_stage.pdf", p, dpi = dpi, width = 6, height = 12, units = "in")



### Do trajectory

### Monocle
# devtools::load_all("/mnt/raid61/Personal_data/huangfei/UCB/git_database/monocle3/")

# devtools::load_all("/mnt/raid62/Lung_cancer_10x/monocle3")

meta = selected_obj@meta.data

meta$Stage[!meta$Disease %in% c("LUSC", "LUAD")] = meta$Disease[!meta$Disease %in% c("LUSC", "LUAD")]


print(unique(meta$Stage))

num_cores = 10
low_thresh = 1 ### the threshold of low-expressed
reduced_element = 1 ### the con-founding factors need to be reduced (batch effect)
num_PCs =  n.pcs
relative_expr = TRUE ###TRUE when using counts; FALSE when treating PSI
norm_method = 'log' ### 'log' when using counts; 'none' when treating PSI
reduction_method = 'both' ### the method of dimension reduction: tSNE and UMAP, or both
RGE_method = 'DDRTree' ### method to learn cell trajectory (recommoned)
cluster_method = 'densityPeak'


temp = selected_obj@dr$umap@cell.embeddings
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

feature_info = data.frame(row.names=rownames(selected_obj@raw.data[,choosen]), gene_short_name=rownames(selected_obj@raw.data[,choosen]))

if (file.exists("monocle3.rds") && trajectory != "DO") {
    cds <- readRDS("monocle3.rds")
} else {

    ### create monocle object
    cds <- new_cell_data_set(
        as.matrix(selected_obj@raw.data[,choosen]),
        cell_metadata = selected_obj@meta.data[choosen,],
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
}


# pData(cds)[,'Seurat_Cluster'] <- selected_obj@ident[colnames(cds)]
# p1 = plot_cells(cds,alpha=0.5, color_cells_by='PatientID')
# p2 = plot_cells(cds,alpha=0.5, color_cells_by='batch')

# p <- cowplot::plot_grid(p1, p2)
# ggsave(filename = "monocle3_umap_paitent_batch.pdf", plot = p, width = 12, height = 6, dpi = 600, units = "in")


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


## Slingshot
if(file.exists("slingshot.rds") && trajectory != "DO") {
    sce = readRDS("slingshot.rds")
} else {
    sce <- slingshot(
        cds,
        clusterLabels = 'Stage', 
        reducedDim = 'UMAP', 
        start.clus=root, 
        end.clus = stages[[length(stages)]]
    )

    saveRDS(sce, file = "slingshot.rds")
}

# rownames(cols_for_stage) <- cols_for_stage$Cells

legends = unique(sort(selected_obj@meta.data$Stage))
cols_for_legend = sapply(legends, function(x){colors[x]})


cols_for_stage = as.data.frame(stage_colors)
cols_for_stage$Stage = rownames(cols_for_stage)
cols_for_stage = merge(selected_obj@meta.data, cols_for_stage, by = "Stage")

rownames(cols_for_stage) <- cols_for_stage$Cells

print(head(cols_for_stage))

tempUMAP = reducedDims(sce)$UMAP

print("slingshot")
pdf("slingpseudotime.pdf", width = 6, height = 6)

plot(
    tempUMAP, 
    col = as.character(cols_for_stage[rownames(tempUMAP), "stage_colors"]), 
    pch=16, 
    asp = 1, 
    cex = 0.8
)
lines(SlingshotDataSet(sce), lwd = 3)

# legend(
#     "topright", 
#     legends, 
#     col=as.character(cols_for_legend),
#     text.col="black", 
#     pch=16, 
#     bty="n", 
#     cex=1,
#     inset=c(-0.2,0)
# )

dev.off()


### plot first path
# dir.create("slingPseudotime", showWarnings=F)
# idx = 1
# for (i in grep("slingPseudotime", colnames(colData(sce)))) {
#     print(i)
#     pdf(paste0("slingPseudotime/", idx, ".pdf"), width = 6, height = 6)
    
#     colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
#     plot(reducedDims(sce)$UMAP, col = colors[cut(colData(sce)[, i], breaks=100)], pch=16, asp = 1, cex = 0.8)
#     lines(SlingshotDataSet(sce), lwd=2)
    
#     dev.off()
#     idx = idx + 1
# }



### plot by disease
pdf("slingpseudotime_disease.pdf", width = 6, height = 6)

cols_for_stage = as.data.frame(disease_colors)
colnames(cols_for_stage) <- "col"
cols_for_stage$Disease = rownames(cols_for_stage)
cols_for_stage = merge(selected_obj@meta.data, cols_for_stage, by = "Disease")

rownames(cols_for_stage) <- cols_for_stage$Cells

legends = sort(unique(selected_obj@meta.data$Disease))
cols_for_legend = sapply(legends, function(x){disease_colors[x]})

tempUMAP = reducedDims(sce)$UMAP
plot(tempUMAP, col = as.character(cols_for_stage[rownames(tempUMAP), "col"]), pch=16, asp = 1, cex = 0.8)
lines(SlingshotDataSet(sce), lwd = 3)


# legend("right", 
#     legends,
#     col=cols_for_stage, 
#     text.col="black", pch=16, bty="n", cex=1)

dev.off()


# ### KEGG、GO and DOSE
# egs <- get_entrzid(markers_by_stage)


# kegg_stage = NULL
# do_stage = NULL
# go_stage = NULL
# for(i in unique(egs$ident)) {
#     print(i)
#     eg <- egs[egs$ident == i,]

#     temp = do_kegg(eg, cluster = i)
#     kegg_stage = rbind(kegg_stage, temp)

#     temp = do_do(eg, cluster = i)
#     do_stage = rbind(do_stage, temp)

#     temp = do_go(eg, cluster = i)
#     go_stage = rbind(go_stage, temp)
# }


# egs <- get_entrzid(markers_by_cluster)


# kegg_clt = NULL
# do_clt = NULL
# go_clt = NULL
# for(i in unique(egs$ident)) {
#     print(i)
#     eg <- egs[egs$ident == i,]

#     temp = do_kegg(eg, cluster = i)
#     kegg_clt = rbind(kegg_clt, temp)

#     temp = do_do(eg, cluster = i)
#     do_clt = rbind(do_clt, temp)

#     temp = do_go(eg, cluster = i)
#     go_clt = rbind(go_clt, temp)
# }


# ### 4. save marker genes, kegg, go and dose results
# wb = createWorkbook()
# addWorksheet(wb, "marker genes")
# writeData(wb, 1, markers_by_stage, rowNames = TRUE)

# addWorksheet(wb, "KEGG")
# writeData(wb, 2, kegg_stage)

# addWorksheet(wb, "GO")
# writeData(wb, 3, go_stage)

# addWorksheet(wb, "DOSE")
# writeData(wb, 4, do_stage)

# saveWorkbook(wb, file = "markers_by_stage.xlsx", overwrite = T)


# wb = createWorkbook()
# addWorksheet(wb, "marker genes")
# writeData(wb, 1, markers_by_cluster, rowNames = TRUE)

# addWorksheet(wb, "KEGG")
# writeData(wb, 2, kegg_clt)

# addWorksheet(wb, "GO")
# writeData(wb, 3, go_clt)

# addWorksheet(wb, "DOSE")
# writeData(wb, 4, do_clt)

# saveWorkbook(wb, file = "markers_by_cluster.xlsx", overwrite = T)

