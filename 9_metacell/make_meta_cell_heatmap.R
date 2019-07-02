library(Seurat)
library(reshape2)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(openxlsx)


make_heatmap_seurat <- function(obj, data, group.by = "res.0.6") {
    meta = obj@meta.data[order(obj@meta.data[, group.by]), ]

    temp = obj@scale.data[
        markers$gene[markers$gene %in% rownames(data)], 
        rownames(meta)
    ]
    
    if (group.by == "Stage") {
        ra = rowAnnotation(
            stage = factor(markers$ident[markers$gene %in% rownames(data)])
        )
        ha = HeatmapAnnotation(
            stage = factor(meta[, group.by])
        )
    } else {
        ra = rowAnnotation(
            cluster = factor(markers$ident[markers$gene %in% rownames(data)])
        )
        ha = HeatmapAnnotation(
            cluster = factor(meta[, group.by])
        )
    }

    
    Heatmap(temp, right_annotation = ra, top_annotation = ha, cluster_rows = F, cluster_columns=F, name = "mat", show_column_names = F)

}


make_heatmap <- function(data, markers) {
    temp = data[markers$gene[markers$gene %in% rownames(data)], ]
    ra = rowAnnotation(cluster = as.factor(markers$ident[markers$gene %in% rownames(data)]))
    
    Heatmap(temp, right_annotation = ra, cluster_rows = F, name = "mat")
}


args = commandArgs(trailingOnly = T)
rds = args[1]
scmat = args[2]
mc = args[3]
gstat = args[4]
output_dir = args[5]

dir.create(output_dir, showWarnings = F, recursive = T)


obj <- readRDS(rds)


load(scmat)
scmat = object


load(mc)
mc = object

load(gstat)
gstat = object

out_df = rbind(tapply(colSums(scmat@mat[, names(mc@mc)]), mc@mc, mean), table(mc@mc))
out_df = cbind(rep("", nrow(out_df)), out_df)


fp_max = apply(mc@mc_fp, 1, max)
fp_tot = gstat[intersect(rownames(mc@mc_fp), rownames(gstat)), "tot"]


f = fp_max > 0 & fp_tot > 0

lfp = round(log2(mc@mc_fp[f,]), 2)
colnames(lfp) <- paste("M", colnames(lfp))



markers <- read.xlsx(paste(dirname(rds), "annotation_results_by_cluster.xlsx", sep = "/"), rowNames = T)

markers <- as.data.frame(markers %>% filter(p_val_adj < 0.05) %>% group_by(ident) %>% top_n(10, wt = avg_logFC))


png(filename = paste(output_dir, "metacell_top10_cluster.png", sep = "/"), width = 8, height = nrow(markers) / 6, res = 600, units = "in")
make_heatmap(lfp, markers)
dev.off()

png(filename = paste(output_dir, "obj_top10_cluster.png", sep = "/"), width = 8, height = nrow(markers) / 6, res = 600, units = "in")
make_heatmap_seurat(obj, lfp)
dev.off()

markers <- read.xlsx(paste(dirname(rds), "annotation_results_by_stage.xlsx", sep = "/"), rowNames = T)

markers <- as.data.frame(markers %>% filter(p_val_adj < 0.05) %>% group_by(ident) %>% top_n(10, wt = avg_logFC))


png(filename = paste(output_dir, "metacell_top10_stage.png", sep = "/"), width = 8, height = nrow(markers) / 6, res = 600, units = "in")
make_heatmap(lfp, markers)
dev.off()


png(filename = paste(output_dir, "obj_top10_stage.png", sep = "/"), width = 8, height = nrow(markers) / 6, res = 600, units = "in")
make_heatmap_seurat(obj, lfp, group.by = "Stage")
dev.off()


clt = as.data.frame(mc@mc)
colnames(clt) = "M"

clt$Stage = obj@meta.data[rownames(clt), "Stage"]
clt$clt = obj@meta.data[rownames(clt), "res.0.6"]


temp = clt %>% 
    select(M, clt) %>% 
    unique() %>%
    group_by(M, clt) %>% 
    add_tally() %>% 
    group_by(M) %>% 
    mutate(perc = n / sum(n)) %>%
    select(M, clt, perc) %>% 
    unique()

temp = dcast(temp, M~clt, value.var = "perc")

temp$M <- paste0("M", temp$M)

rownames(temp) <- temp$M

temp = temp[, colnames(temp) != "M"]


png(filename = paste(output_dir, "metacell_cluster_overlap.png", sep = "/"), width = 8, height = 8, res = 600, units = "in")


if (length(unique(melt(temp))) == 2) {
    temp[is.na(temp)] = 0
}

Heatmap(
    temp, 
    cluster_rows = F, 
    cluster_columns = F, 
    name = "perc",
    cell_fun = function(j, i, x, y, width, height, fill) {
        if(!is.na(temp[i, j ]))
            grid.text(sprintf("%.2f", temp[i, j]), x, y, gp = gpar(fontsize = 10))
    }
)
dev.off()



temp = clt %>% 
    select(M, Stage) %>% 
    unique() %>%
    group_by(M, Stage) %>% 
    add_tally() %>% 
    group_by(M) %>% 
    mutate(perc = n / sum(n)) %>%
    select(M, Stage, perc) %>% 
    unique()

temp = dcast(temp, M~Stage, value.var = "perc")

temp$M <- paste0("M", temp$M)

rownames(temp) <- temp$M

temp = temp[, colnames(temp) != "M"]


png(filename = paste(output_dir, "metacell_stage_overlap.png", sep = "/"), width = 8, height = 8, res = 600, units = "in")

if (length(unique(melt(temp))) == 2) {
    temp[is.na(temp)] = 0
}

Heatmap(
    temp, 
    cluster_rows = F, 
    cluster_columns = F, 
    name = "perc",
    cell_fun = function(j, i, x, y, width, height, fill) {
        if(!is.na(temp[i, j ]))
            grid.text(sprintf("%.2f", temp[i, j]), x, y, gp = gpar(fontsize = 10))
    }
)
dev.off()


