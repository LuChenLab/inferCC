dpi = 600

set.seed(1)

# load packages
source("/mnt/raid62/Lung_cancer_10x/scripts/each/basic.R")

use_python("/usr/bin/python")


# if the image direcotory not exists, create it
suppressPackageStartupMessages(load_pacakges())


use_python("/usr/bin/python")


# if the image direcotory not exists, create it
image_dir = here("02_figures_each_cell/Granulocyte_ADC")
input_rds = here("00_data_ingest/04_rds_generated/seurat_obj/Granulocyte.rds")


dir.create(image_dir, showWarnings = FALSE)
img.units = 6

setwd(image_dir)


################

### 1. Load raw counts and meta info

obj <- readRDS(input_rds)
obj <- chanage_meta_info(obj)


### 2. check is the variables like Stage and Patient are highly variable
meta = obj@meta.data[obj@meta.data$Batch != 3 & obj@meta.data$Disease == "ADC",]

m <- get_smote_match(meta)

### 3. create selected object and make barplots to check the SMOTE effect
selected_obj <- CreateSeuratObject(
    raw.data = obj@raw.data[, rownames(m)],
    meta.data = obj@meta.data[rownames(m),]
)


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

ggsave("stats_SMOTE.png", plot = p, width = 12, height = 6, dpi=dpi, units = "in")


### 4. basic 

#### 1. vlnplot
mito.genes <- grep(pattern = "^MT-", x = rownames(x = selected_obj@raw.data), value = TRUE)
percent.mito <- Matrix::colSums(selected_obj@raw.data[mito.genes, ])/Matrix::colSums(obj@raw.data)

# AddMetaData adds columns to object@meta.data, and is a great place to
# stash QC stats
selected_obj <- AddMetaData(object =selected_obj, metadata = percent.mito, col.name = "percent.mito")

png(filename = "VlnPlot.png", width = 18, height = 6, res = dpi, units = "in")
VlnPlot(object = selected_obj, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
dev.off()


#### 2. GenePlot

png(filename = "GenePlot.png", width = 12, height = 6, res = dpi, units = "in")
par(mfrow = c(1, 2))
GenePlot(object = selected_obj, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = selected_obj, gene1 = "nUMI", gene2 = "nGene")
dev.off()


#### 3. do PCA and harmony
selected_obj <- NormalizeData(object = selected_obj, normalization.method = "LogNormalize", scale.factor = 10000)
selected_obj <- FindVariableGenes(object = selected_obj, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
selected_obj <- ScaleData(object = selected_obj, vars.to.regress = c("nUMI", "percent.mito"))
selected_obj <- RunPCA(object = selected_obj, pc.genes = selected_obj@var.genes, do.print = FALSE, pcs.compute = 100)

selected_obj <- RunHarmony(selected_obj, c("batch", "PatientID"))


#### 4. make dotplot of SD of PCAs and Harmonys
png(filename = "PCElbowPlot.png", width = 12, height = 6, res = dpi, units = "in")
PCElbowPlot(selected_obj, num.pc=100)  +
    scale_x_continuous(breaks=1:100) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
dev.off()


temp <- as.data.frame(apply(selected_obj@dr$harmony@cell.embeddings, 2, sd))
colnames(temp) <- c("sd")
temp$Harmony = 1:nrow(temp)

p <- ggplot(temp, aes(x=Harmony, y=sd)) + 
    geom_point() + 
    ylab("Standard Dei of Harmony")  +
    scale_x_continuous(breaks=1:100) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave("HarmonyPlot.png", plot = p, width = 12, height = 6, dpi = dpi, units = "in")


#### 5. tSNE and UMAP
n.pcs = 10

selected_obj <- RunTSNE(object = selected_obj, dims.use = 1:n.pcs, do.fast = TRUE, reduction.use = "harmony", reduction.name = "tsne")
selected_obj <- RunUMAP(object = selected_obj, dims.use = 1:n.pcs, do.fast = TRUE, reduction.use = "harmony")
selected_obj <- FindClusters(object = selected_obj, reduction.type = "harmony", dims.use = 1:n.pcs, resolution = 0.6, print.output = 0, save.SNN = FALSE, force.recalc = TRUE)

#### 6. change cluster id from 0 to 1, and save seurat object
selected_obj@meta.data$res.0.6 <- as.numeric(selected_obj@meta.data$res.0.6) + 1
saveRDS(selected_obj, basename(input_rds))


#### 7. Dump seurat object data to several files for scanpy
dump_seurat_for_scanpy(selected_obj, "paga")


### check cluster and batch effect by tsne and umap

p <- myDimPlot(
    selected_obj, 
    reduction.use = "tsne", 
    group.by = "res.0.6", 
    no.legend = TRUE, 
    do.label = TRUE, 
    label.size = 5, 
    do.return = TRUE, 
    pt.size = 0.01
)
ggsave("tsne_cluster.png", plot=p, width=6, height=6, dpi=dpi, units="in")


p <- myDimPlot(
    selected_obj, 
    reduction.use = "tsne", 
    group.by = "batch", 
    no.legend = FALSE, 
    do.label = FALSE, 
    do.return = TRUE, 
    pt.size = 0.01
) + theme(legend.position = c(0.1, 0.1))
ggsave("tsne_batch.png", plot=p, width=6, height=6, dpi=dpi, units="in")


p <- myDimPlot(
    selected_obj, 
    reduction.use = "umap", 
    group.by = "res.0.6", 
    no.legend = TRUE, 
    do.label = TRUE, 
    label.size = 5, 
    do.return = TRUE, 
    pt.size = 0.01
)
ggsave("umap_cluster.png", plot=p, width=6, height=6, dpi=dpi, units="in")

p <- myDimPlot(
    selected_obj, 
    reduction.use = "umap", 
    group.by = "batch", 
    no.legend = FALSE, 
    do.label = FALSE, 
    do.return = TRUE, 
    pt.size = 0.01
) + theme(legend.position = c(0.1, 0.1))
ggsave("umap_batch.png", plot=p, width=6, height=6, dpi=dpi, units="in")


### 
#### 1. tsne by disease
p <- make_tsne_plot(selected_obj, pt.size = 0.5, colors=colors) + theme(
    legend.position = c(0.85, 0.85),
    legend.title = element_blank(),
    legend.background=element_blank(),
    legend.key = element_rect(colour = "transparent", fill = "transparent")
)

ggsave(filename = "tsne_by_disease.png", plot = p, width = 6, height = 6, dpi = dpi, units = "in")



#### 2. tsne by stage
p <- make_tsne_plot(selected_obj, color="Stage", pt.size = 0.5, colors=colors) + 
    theme(
        legend.position = c(0.85, 0.85),
        legend.title = element_blank(),
        legend.background=element_blank()
    )
ggsave("tsne_by_stage.png", plot = p, width = 6, height = 6, units = "in", dpi = dpi)



batch = sort(unique(selected_obj@meta.data$Stage))
for (i in 1:length(batch)){
    print(i)
    
    stage = selected_obj@meta.data[selected_obj@meta.data$Stage == batch[i], ]
    p <- myDimPlot(
        selected_obj, 
        reduction.use = "tsne", 
        cells.highlight = rownames(stage), 
        do.return=TRUE,
        no.legend = FALSE,
        pt.size = 0.5,
        cols.highlight = rep(colors[batch[i]], nrow(stage))
    )
    
    temp = selected_obj@dr$tsne@cell.embeddings
    p = p + geom_text(
        x = min(temp[, 1]) + (max(temp[, 1]) - min(temp[, 1])) * 0.01, 
        y = min(temp[, 2]) + (max(temp[, 2]) - min(temp[, 2])) * 0.01, 
        label = batch[i], 
        hjust = 0, 
        vjust = 0
    )
    
    ggsave(paste0("tsne_by_stage_", batch[i], ".png"), plot = p, width = 6, height = 6, units = "in", dpi = dpi)
}


#### 3. tsne by patient

patients = sort(unique(selected_obj@meta.data$PatientID))

p <- make_tsne_plot(selected_obj, color="PatientID", pt.size = 0.5, guide_ncol = 3, colors=colors) + 
    theme(
        legend.position = c(0.85, 0.85),
        legend.title = element_blank(),
        legend.background=element_blank(),
        legend.key = element_rect(colour = NA, fill = NA)
    ) +
    guides(
        color = guide_legend(
            ncol=4, 
            override.aes = list(size = 4, alpha = 1),
            keywidth=0.05,
            keyheight=0.05,
            default.unit="inch"
        )
    )
ggsave("tsne_by_patient.png", plot = p, width = 6, height = 6, units = "in", dpi = dpi)


for (i in 1:length(patients)){
    print(i)
    
    stage = selected_obj@meta.data[selected_obj@meta.data$PatientID == patients[i], ]
    p <- myDimPlot(
        selected_obj, 
        reduction.use = "tsne", 
        cells.highlight = rownames(stage), 
        do.return=TRUE,
        no.legend = FALSE,
        pt.size = 0.5,
        cols.highlight = rep(patient_colors[i], nrow(stage))
    )
    
    temp = selected_obj@dr$tsne@cell.embeddings
    p = p + geom_text(
        x = min(temp[, 1]) + (max(temp[, 1]) - min(temp[, 1])) * 0.01, 
        y = min(temp[, 2]) + (max(temp[, 2]) - min(temp[, 2])) * 0.01, 
        label = patients[i], 
        hjust = 0, 
        vjust = 0
    )
    
    ggsave(paste0("tsne_by_patient_",as.character(patients[i]), ".png"), plot = p, width = 6, height = 6, units = "in", dpi = dpi)
}


#### 4. make combine barplot
p <- make_combined_barplot(selected_obj, rel_height = c(0.3, 1), with_disease = F)
ggsave("barplots.png", plot = p, width = 15, height = 10, dpi = dpi, units = "in", limitsize = FALSE)


### Find markers
markers_by_stage = find_markers(selected_obj, group.by = "Stage", n.cores=1)
markers_by_cluster = find_markers(selected_obj, group.by = "res.0.6", n.cores=1)


#### 5. make dotplot and heatmap of top10 markers
if (ncol(selected_obj@raw.data) > 30000){
    cells = sample(1:ncol(obj@raw.data), 30000)
    cells = colnames(selected_obj@raw.data)[cells]
} else {
    cells = colnames(selected_obj@raw.data)
}


top10 <- markers_by_cluster %>% group_by(ident) %>% top_n(n = 10, wt = avg_logFC)

png("heatmap_cluster.png", width = 8, height = 12, res = dpi, units = "in")
DoHeatmap(object = selected_obj, genes.use = top10$gene, slim.col.label = TRUE, remove.key = TRUE, cells.use = cells, group.by="res.0.6")
dev.off()

p <- DotPlot(selected_obj, genes.plot = unique(top10$gene), do.return = TRUE, x.lab.rot = TRUE) + coord_flip() + theme(axis.text.x = element_text(angle = 0))
ggsave("dotplot_top10_cluster.png", p, dpi = dpi, width = 6, height = 12, units = "in")


top10 <- markers_by_stage %>% group_by(ident) %>% top_n(n = 10, wt = avg_logFC)

png("heatmap_stage.png", width = 8, height = 12, res = dpi, units = "in")
DoHeatmap(object = selected_obj, genes.use = top10$gene, group.by = "Stage", slim.col.label = TRUE, remove.key = TRUE, cells.use = cells)
dev.off()

p <- DotPlot(selected_obj, genes.plot = unique(top10$gene), do.return = TRUE, x.lab.rot = TRUE, group.by = "Stage") + coord_flip() + theme(axis.text.x = element_text(angle = 0))
ggsave("dotplot_top10_stage.png", p, dpi = dpi, width = 6, height = 12, units = "in")



### Monocle
devtools::load_all("/mnt/raid61/Personal_data/huangfei/UCB/git_database/monocle3/")

num_cores = 10
low_thresh = 1 ### the threshold of low-expressed
reduced_element = 1 ### the con-founding factors need to be reduced (batch effect)
num_PCs =  cal_pcs(obj)   ### the numbers of used PCs
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
    if (summ[1] / summ[2] > 30) {
        choosen = choosen & temp[, i] >= summ[2]
    }
    
    if (summ[6] / summ[5] > 30) {
        choosen = choosen & temp[, i] <= summ[5]
    }
}


feature_info = data.frame(row.names=rownames(selected_obj@raw.data[,choosen]), gene_short_name=rownames(selected_obj@raw.data[,choosen]))

### create monocle object
print("create monocle object")
pd <-  new('AnnotatedDataFrame', data = selected_obj@meta.data[choosen,])
#### pre-process the matrix according to the format

fd = new('AnnotatedDataFrame', data = feature_info)

cds <- newCellDataSet(as.matrix(selected_obj@raw.data[,choosen]),
                      phenoData = pd,
                      featureData = fd,
                      lowerDetectionLimit = low_thresh,
                      expressionFamily = negbinomial.size())

print("normalize")
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

cds <- preprocessCDS(cds, num_dim = n.pcs,
                     relative_expr = relative_expr,
                     norm_method = 'log', verbse=T,
                     residualModelFormulaStr = paste0('~',reduced_element))


### replace pca with harmony
harmony_mat <- selected_obj@dr$harmony@cell.embeddings[choosen,1:n.pcs]
cds@normalized_data_projection <- harmony_mat[,]

### dimesion reduction
cds <- reduceDimension(cds, reduction_method = 'tSNE', norm_method = norm_method)


### replace umap with Seurat umap results
reducedDimA(cds)  <- t(selected_obj@dr$umap@cell.embeddings[choosen, ]) -> reducedDimS(cds)
cds <- partitionCells(cds)
cds <- learnGraph(cds,  RGE_method = RGE_method)


cds <- clusterCells(cds,verbose = F,
                    method = cluster_method,
                    res = 1,
                    louvain_iter = 1,
                    cores = num_cores)


pData(cds)[,'Seurat_Cluster'] <- selected_obj@ident[colnames(cds)]
p1 = plot_cell_clusters(cds,alpha=0.5,color_by='PatientID')
p2 = plot_cell_clusters(cds,alpha=0.5,color_by='batch')

p <- cowplot::plot_grid(p1, p2)
ggsave(filename = "monocle3_umap_paitent_batch.png", plot = p, width = 12, height = 6, dpi = 600, units = "in")


### using Benign or Stage I as root
root = "I"

node_ids = get_correct_root_state(cds,cell_phenotype='Stage', root)
cds <- orderCells(cds, root_pr_nodes = node_ids)
p1 = plot_cell_trajectory(cds,alpha=0.5) + labs(title = root[1])
p2 = plot_cell_trajectory(cds,alpha=0.5,color_by='Stage')

p <- cowplot::plot_grid(p1, p2)

ggsave(filename = "monocle3_trajectory.png", plot = p, width = 12, height = 6, dpi = 600, units = "in")


p <- plot_cell_trajectory(cds,alpha=0.5,color_by='Stage') + facet_grid(.~Stage)
ggsave(filename = "monocle3_trajectory_by_stage.png", plot = p, width = 6 * length(unique(cds$Stage)), height = 6, dpi = 600, units = "in")


### find pseudotime temporal genes
pr_graph_test <- principalGraphTest(cds, k=3, cores=10)

pr_graph_signif = setDT(rownames_to_column(pr_graph_test, var = 'gene_symbol'))

pr_graph_test$gene_symbol <- rownames(pr_graph_test)

### choose the significant junctions
pr_graph_signif = pr_graph_signif[morans_I>0.5&qval<0.05][order(qval,-morans_test_statistic)]
cluster_marker_res <- find_cluster_markers(cds, pr_graph_test, group_by = 'Stage', morans_I_threshold = 0.5)
cluster_marker_res

write.xlsx(cluster_marker_res, "monocle3_marker_genes.xlsx")

### need the expressed percentage > 0.2
cluster_marker_res <- setDT(cluster_marker_res)[percentage > 0.25]

top_marker <- cluster_marker_res[order(-specificity,-percentage,-morans_I)][,head(.SD,10),by=Group][order(Group)]
top_marker <- unique(as.character(rev(top_marker$gene_short_name)))

png("monocle3_heatmap_stage_top_marker.png", width = 8, height = 12, res = dpi, units = "in")
DoHeatmap(object = selected_obj, genes.use = top_marker, group.by = "Stage", slim.col.label = TRUE, remove.key = TRUE, cells.use = colnames(selected_obj@raw.data[, choosen]))
dev.off()

if (length(top_marker) > 1) {
    p <- plot_markers_by_group(cds, rev(top_marker), group_by = 'Stage', ordering_type='none')
    ggsave(filename = "monocle3_dotplot_marker_genes.png", plot = p, width = 8, height = 12, dpi = 600, units = "in")
    
    
    p <- plot_cell_clusters(cds, markers = top_marker, ncol=5)
    
    width = 15
    if (length(top_marker) < 5) {
        width = 3 * length(top_marker)
    }
    
    ggsave(filename = "monocle3_umap_marker_genes.png", plot = p, width = width, height = 3 * ceiling(length(top_marker) / 5), dpi = 600, units = "in")
}

saveRDS(cds, "monocle3.rds")


### Slingshot
library(slingshot)
if(sum(choosen) < ncol(selected_obj@raw.data)) {
    temp = exportCDS(cds, export_to = "Seurat" , export_all = TRUE)
    sling = as.SingleCellExperiment(temp)
} else {
    sling = as.SingleCellExperiment(selected_obj)
}

stages = sort(unique(selected_obj@meta.data$Stage[choosen]))

print(stages)

reducedDims(sling)$UMAP <- selected_obj@dr$umap@cell.embeddings[choosen,]
sce <- slingshot(sling, clusterLabels = 'Stage', reducedDim = 'UMAP', start.clus="I", end.clus = "IV")


legends = c("I", "II", "III",'IV')
cols_for_legend = sapply(legends, function(x){colors[x]})


cols_for_stage = as.data.frame(colors)
cols_for_stage$Stage = rownames(cols_for_stage)
cols_for_stage = merge(selected_obj@meta.data, cols_for_stage, by = "Stage")

rownames(cols_for_stage) <- cols_for_stage$Cells


png("slingpseudotime.png", width = 6, height = 6, res = dpi, units = "in")

tempUMAP = reducedDims(sce)$UMAP
plot(tempUMAP, col = as.character(cols_for_stage[rownames(tempUMAP), "colors"]), pch=16, asp = 1, cex = 0.8)
lines(SlingshotDataSet(sce), lwd = 3)

legend("bottomright", legends, col=cols_for_legend, 
       text.col="black", pch=16, bty="n", cex=1)

dev.off()



### plot single path
idx = 1
for (i in grep("slingPseudotime", colnames(colData(sce)))) {
    print(i)
    png(paste0("slingpseudotime_", idx, ".png"), width = 6, height = 6, res = dpi, units = "in")
    
    colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
    plot(reducedDims(sce)$UMAP, col = colors[cut(colData(sce)[, i], breaks=100)], pch=16, asp = 1, cex = 0.8)
    lines(SlingshotDataSet(sce), lwd=2)
    
    dev.off()
    idx = idx + 1
}



### make heatmap by top genes
#### plot heatmap

require(gam)
slingshot_topgenes = c()

idx = 1
for (i in grep("slingPseudotime", colnames(colData(sce)))) {
    print(i)
    t <- colData(sce)[, i]
    
    # for time, only look at the 100 most variable genes
    Y <- assays(sce)$logcounts
    var1K <- names(sort(apply(Y,1,var),decreasing = TRUE))[1:100]
    Y <- Y[var1K,]
    
    # fit a GAM with a loess term for pseudotime
    gam.pval <- apply(Y,1,function(z){
        d <- data.frame(z=z, t=t)
        tmp <- gam(z ~ lo(t), data=d)
        p <- summary(tmp)[4][[1]][1,5]
        p
    })
    
    topgenes <- names(sort(gam.pval, decreasing = FALSE))[1:100]
    p <- DoHeatmap(selected_obj, genes.use = topgenes, group.by = "Stage", slim.col.label = TRUE, remove.key = TRUE, cells.use = colnames(selected_obj@raw.data[, choosen]), do.plot = F)
    ggsave(paste0("slingshot_pseudotime_", idx,"_heatmap.png"), width = 8, height = 12, units = "in", dpi = dpi)
    
    idx = idx + 1
    slingshot_topgenes = c(slingshot_topgenes, topgenes)
}


saveRDS(sce, file = "slingshot.rds")



### KEGG、GO and DOSE

egs <- get_entrzid(markers_by_stage)

registerDoMC(1)
kegg_stage <- foreach(i=unique(egs$ident), .combine=rbind) %dopar% {
    print(i)
    eg <- egs[egs$ident == i,]
    do_kegg(eg, cluster = i)
}

registerDoMC(1)
do_stage <- foreach(i=unique(egs$ident), .combine=rbind) %dopar% {
    print(i)
    eg <- egs[egs$ident == i,]
    do_do(eg, cluster = i)
}

registerDoMC(1)
go_stage <- foreach(i=unique(egs$ident), .combine=rbind) %dopar% {
    print(i)
    eg <- egs[egs$ident == i,]
    do_go(eg, cluster = i)
}


# eg <- bitr(scanpy[scanpy$cluster == cluster[1], "names"], fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

egs <- get_entrzid(markers_by_cluster)

registerDoMC(1)
kegg_cluster <- foreach(i=unique(egs$ident), .combine=rbind) %dopar% {
    print(i)
    eg <- egs[egs$ident == i,]
    do_kegg(eg, cluster = i)
}

registerDoMC(1)
do_cluster <- foreach(i=unique(egs$ident), .combine=rbind) %dopar% {
    eg <- egs[egs$ident == i,]
    do_do(eg, cluster = i)
}

registerDoMC(1)
go_cluster <- foreach(i=unique(egs$ident), .combine=rbind) %dopar% {
    print(i)
    eg <- egs[egs$ident == i,]
    do_go(eg, cluster = i)
}


#### 4. save marker genes, kegg, go and dose results
wb = createWorkbook()
addWorksheet(wb, "marker genes")
writeData(wb, 1, markers_by_stage, rowNames = TRUE)

addWorksheet(wb, "KEGG")
writeData(wb, 2, kegg_stage)

addWorksheet(wb, "GO")
writeData(wb, 3, go_stage)

addWorksheet(wb, "DOSE")
writeData(wb, 4, do_stage)

saveWorkbook(wb, file = "annotation_results_by_stage.xlsx", overwrite = T)


wb = createWorkbook()
addWorksheet(wb, "marker genes")
writeData(wb, 1, markers_by_cluster, rowNames = TRUE)

addWorksheet(wb, "KEGG")
writeData(wb, 2, kegg_cluster)

addWorksheet(wb, "GO")
writeData(wb, 3, go_cluster)

addWorksheet(wb, "DOSE")
writeData(wb, 4, do_cluster)

saveWorkbook(wb, file = "annotation_results_by_cluster.xlsx", overwrite = T)



####### 
##### 

library(TCGAbiolinks)
library(dplyr)
library(DT)
library(SummarizedExperiment)


dir.create("TCGA", showWarnings = F)

# markers_by_stage = read.xlsx("annotation_results_by_stage.xlsx", rowNames = T)
# markers_by_cluster = read.xlsx("annotation_results_by_cluster.xlsx", rowNames = T)

# markers, dataframe
get_top <- function(data, n = 10) {
    res = c()
    for(i in unique(data$ident)) {
        temp = data[data$ident == i, ]
        temp = temp[order(temp$avg_logFC), "gene"]
        res = c(res, temp[1:n])
    }
    
    return(na.omit(res))
}

top10 <- as.data.frame(markers_by_cluster %>% group_by(ident) %>% top_n(n = 10, wt = avg_logFC))$gene
top10 <- c(top10, as.data.frame(markers_by_stage %>% group_by(ident) %>% top_n(n = 10, wt = avg_logFC))$gene )

top10 <- unique(c(top10, top_marker, slingshot_topgenes))

setwd(here("TCGA"))

# query <- GDCquery(
#     project = "TCGA-LUAD",
#     data.category = "Gene expression",
#     data.type = "Gene expression quantification",
#     platform = "Illumina HiSeq", 
#     file.type  = "results",
#     experimental.strategy = "RNA-Seq",
#     legacy = TRUE
# )
# 
# 
# data <- GDCprepare(query)

data = readRDS("LUAD.rds")

# function to make suvival plots
# :param TCGAbiolinks prepared RangedSummarizedExperiments
generate_tcga_suvplots <- function(data, genes, n.cores=10) {
    genes_known = as.character(sapply(rownames(data), function(x){strsplit(x, "\\|", perl=T)[[1]][1]}))
    
    data_raw <- assay(data, "raw_count")
    genes = intersect(genes_known, genes)
    expression_genes = data_raw[which(genes_known %in% genes),]
    expression_genes = t(expression_genes)
    colnames(expression_genes) <- as.character(sapply(colnames(expression_genes), function(x){strsplit(x, "\\|", perl=T)[[1]][1]}))
    
    expression_genes = merge(as.data.frame(colData(data)), expression_genes, by = 0)
    
    registerDoMC(n.cores)
    
    foreach(i=genes, .combine=rbind) %dopar% {
        expression_genes[, i] <- ifelse(expression_genes[, i] > median(expression_genes[, i]),'high','low')
        
        TCGAanalyze_survival(
            expression_genes, 
            i,
            main = i, 
            height = 10, 
            width=10,
            filename = paste(image_dir, "TCGA", paste0(i, ".pdf"), sep = "/")
        )
    }
}


generate_tcga_suvplots(data, top10, n.cores=5)

