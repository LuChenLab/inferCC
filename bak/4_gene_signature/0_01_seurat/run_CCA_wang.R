### Created by Zhang at 2019.01.23
### This script is used to run CCA of wang
library(openxlsx)
library(Seurat)
library(logging)
library(gridExtra)

setwd(output)


meta <- read.xlsx("LungCancer10X_table1.xlsx", sheet = 3)

# meta <- meta[order(meta$ID),]

## load from Cell ranger results
#### 后期还需要取出下batch effect
# ScaleData(object = control.combined, vars.to.regress = c("batchid", "nUMI", "percent.mito"))

load_from_cell_ranger <- function(input, group, min.cells=3, min.genes=500) {
    obj <- Read10X(data.dir = input)
    obj <- CreateSeuratObject(
        raw.data=obj, 
        min.cells = min.cells,
        min.genes=min.genes, 
        project="LungCancer10X"
    )
    obj@meta.data$group <- group
    obj@meta.data$sample_id <- basename(input)
    obj <- NormalizeData(obj)
    obj <- ScaleData(obj)
    obj <- FindVariableGenes(obj, do.plot = F)
    return(obj)
}



load_by_meta_info <- function(meta, input) {
    idx = 1
    objects <- list()
    meta <- na.omit(meta)
    
    for (i in 1:nrow(meta)) {
        
        input_dir = paste(input, meta[i, "sample_id"], sep = "/")
        loginfo(input_dir)
        
        objects[[idx]] <- load_from_cell_ranger(
            input = input_dir, 
            group = meta[i, "tissue"]
        )
        
        idx = idx + 1
        
    }
    
    return (objects)
}


find_genes_use <- function(objects, top=2000) {
    genes.use = c()
    
    for (i in objects) {
        genes.use = c(genes.use, head(rownames(i@hvg.info), 2000))
    }
    
    genes.use = names(which(table(genes.use) > 1))

    for (i in objects) {
        genes.use = genes.use[genes.use %in% rownames(i@scale.data)]
    }
    
    return(genes.use)
}


get_cell_ids <- function(objects) {
    cell_ids = c()
    for(i in objects) {
        cell_ids = unique(c(cell_ids, i@meta.data$sample_id))
    }
    
    return(cell_ids)
}

run_CCA <- function(objects, cell_ids=NA, ccs = 15) {
    genes.use <- find_genes_use(objects)
    
    if(length(cell_ids) == 1 && is.na(cell_ids)) {
        cell_ids = get_cell_ids(objects)
    }
    
    obj <- RunMultiCCA(
        objects, 
        genes.use = genes.use, 
        num.ccs = 15,
        add.cell.ids = cell_ids
        )
    
    return(obj)
}

## run CCA between groups
objects_per_group <- load_by_meta_info(meta, input)
objects_per_group <- run_CCA(objects_per_group)
saveRDS(objects_per_group, file = "objects_cca_between_groups.rds")


## run CCA inside each group, eg: all normal samples
objects <- list()

idx = 1
for (i in unique(meta$tissue)) {
    temp_meta = meta[meta$tissue == i,]
    
    temp_data <- load_by_meta_info(temp_meta, input)
    
    if (length(temp_data) > 2) {
        objects[[idx]] <- run_CCA(temp_data)
    } else if (length(temp_data) == 2) {
        objects[[idx]] <- RunCCA(
            temp_data[[1]], 
            temp_data[[2]],
            add.cell.id1 = temp_data[[1]]@meta.data$sample_id,
            add.cell.id2 = temp_data[[2]]@meta.data$sample_id
            )
    } else {
        objects[[idx]] <- temp_data[[1]]
    }
    
    idx = idx + 1
}

saveRDS(object = objects, file = "objects.rds")

## the run CCA needs 
recalculate <- function(objects) {
    for (i in 1:length(objects)) {
        obj <- NormalizeData(objects[[i]])
        obj <- ScaleData(obj)
        objects[[i]] <- FindVariableGenes(obj, do.plot = F)
    }
    return(objects)
}

# check whether which is which
obj_cca <- run_CCA(
    obj,
    cell_ids = c("Adenocarcinoma", "Normal", "Lymph", "Squamous_cell_carcinoma", "Bronichi")    
)

saveRDS(object = obj_cca, file = "objects_cca.rds")

obj <- objects[[1]]
for ( i in 2:length(objects)) {
    obj <- MergeSeurat(
        obj, objects[[i]],
        add.cell.id1 = obj@meta.data$sample_id,
        add.cell.id2 = objects[[i]]@meta.data$sample_id
    )
}

saveRDS(object = objects, file = "objects_merged.rds")


### Vlnplot for QC
# visualize results of CCA plot CC1 versus CC2 and look at a violin plot
p1 <- DimPlot(
    object = obj_cca,
    reduction.use = "cca",
    group.by = "batch",
    pt.size = 0.5,
    do.return = TRUE
)

p2 <- VlnPlot(object = objects, features.plot = "CC1", group.by = "group",
              do.return = TRUE)
p <- plot_grid(p1, p2)

ggsave(filename="CCA.png", plot=p)


obj <- AlignSubspace(obj, reduction.type = "cca", grouping.var = "group", dims.align = 1:20)
obj <- RunTSNE(obj, reduction.use = "cca.aligned", dims.use = 1:20, do.fast = T)
obj <- FindClusters(obj, reduction.type = "cca.aligned", resolution = 0.6, dims.use = 1:20)

saveRDS(obj, "clustered_between_group.rds")

## load markers genes
# markers <- read.xlsx("/mnt/data5/zhangyiming/seurat/results/LungCancer10X_table1.xlsx", sheet=4)
# markers <- markers[markers$gene %in% gene.use, ]
# 
pdf("tSEN_between_group.pdf", height = 20, width = 15)
p1 <- TSNEPlot(obj, do.return = T, pt.size = 0.5, group.by = "group")
p2 <- TSNEPlot(obj, do.label = T, do.return = T, pt.size = 0.5)
grid.arrange(p1, p2, ncol=1)
dev.off()

library(dplyr)
library(Seurat)

# min.pct = 0.5, logfc.threshold = 1s
All.markers = FindAllMarkers(object = obj, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 1)
All.top10.markers = All.markers  %>% group_by(cluster) %>% top_n(10, avg_logFC)
length(All.top10.markers$gene)

write.table(All.markers, "all_markers_between_groups.txt", quote = F, row.names = F, sep = "\t")
write.table(All.top10.markers, "markers_between_groups.txt", quote = F, row.names = F, sep = "\t")

# gene list from 王成弟
# "SCGB1A1", "KRT5", "P63", "IGTA6", "NGFR", "KRT14", "DNAI1", "MUC5AC", "ASGR1", "ATP6V0D2", "ATP6V1C2", "CFTR", "OXI1", "GM933", "MOXD1", "P2RY14", "STAP1"   

pdf("features_each_group.pdf", height = 10, width = 10)
FeaturePlot(
    object = obj,
    features.plot = c("SCGB1A1", "KRT5", "NGFR", "KRT14", "DNAI1", "MUC5AC", "ASGR1", "ATP6V0D2", "ATP6V1C2", "CFTR", "MOXD1", "P2RY14", "STAP1"),
    min.cutoff = "q9",
    cols.use = c("lightgrey", "blue"),
    pt.size = 0.5
    )
dev.off()


pdf("features_between_group.pdf", height = 10, width = 10)
FeaturePlot(
    object = obj,
    features.plot = c("SCGB1A1", "KRT5",  "NGFR", "KRT14", "DNAI1", "MUC5AC", "ASGR1", "ATP6V0D2", "ATP6V1C2", "CFTR", "MOXD1", "P2RY14", "STAP1"   ),
    min.cutoff = "q9",
    cols.use = c("lightgrey", "blue"),
    pt.size = 0.5
)
dev.off()

# c("GYPA", "ALAS2", "HBB", "MS4A1","CD79A","CD3D","CD44","CD4","CD8A", "FOXP3", "NKG7", "LYZ", "PTPRC", "PF4", "CD34", "ENG", "NT5E", "ITGB1","MME","GP9")
genes <- c("HBB", "MS4A1","CD79A","CD3D","CD44","CD4","CD8A", "FOXP3", "NKG7", "LYZ", "PTPRC", "PF4", "CD34", "ENG", "NT5E", "ITGB1","MME","GP9")

pdf("features_each_group_blood.pdf", height = 10, width = 10)
FeaturePlot(
    object = obj,
    features.plot = genes,
    min.cutoff = "q9",
    cols.use = c("lightgrey", "blue"),
    pt.size = 0.5
)
dev.off()


genes <- c(
    "MS4A1",
    "CD3D",
    "CD44",
    "CD4",
    "CD8A", 
    "FOXP3", 
    "NKG7", 
    "LYZ",
    "PTPRC",
    "CD34", 
    "ENG",
    "NT5E",
    "ITGB1",
    "SCGB1A1",
    "KRT5",
    "ASGR1",
    "MUC5AC",
    "STAP1",
    "DNAI1",
    "MUC5AC"
)

pdf("features_each_group_known.pdf", height = 15, width = 15)
FeaturePlot(
    object = obj,
    features.plot = genes,
    min.cutoff = "q9",
    cols.use = c("lightgrey", "blue"),
    pt.size = 0.5
)
dev.off()



pdf("features_each_group_tumor.pdf", height = 2, width = 5)
FeaturePlot(
    object = obj,
    features.plot = c("TP53", "KRAS"),
    min.cutoff = "q9",
    cols.use = c("lightgrey", "blue"),
    pt.size = 0.5
)
dev.off()



pdf("features_each_group_fiblast.pdf", height = 10, width = 10)
FeaturePlot(
    object = obj,
    features.plot = c("CXCL12", "MUSTN1", "COL5A1", "LYSMD2", "CYGB", "FKBP9", "TIMP3", "PDPN", "CILP", "PID1"),
    min.cutoff = "q9",
    cols.use = c("lightgrey", "blue"),
    pt.size = 0.5
)
dev.off()


current.clusters <- c(
    7, 
    13, 
    3, 20, 
    5, 8, 6, 10,
    4, 9, 2, 12, 22,
    0, 1, 17, 14,
    15, 11
    )
cell.ids <- c(
    "Basal",
    "Club",
    "B", "B",
    "NK", "NK", "NK", "NK",
    "Mono", "Mono", "Mono", "Mono", "Mono",
    "T", "T", "T", "T",
    "Fiblast", "Fiblast"
)
obj1 <- obj
obj1@ident <- plyr::mapvalues(x = obj1@ident, from = current.clusters, to = cell.ids)
pdf("tSNE_anno.pdf",width = 10,height = 7)
TSNEPlot(object = obj1, do.label = TRUE, pt.size = 1.5, label.size=5)
dev.off()


pdf("heatmap_top10.pdf", height = 10, width = 20)
DoHeatmap(object = obj, genes.use = All.top10.markers$gene, slim.col.label = TRUE, remove.key = TRUE)
dev.off()

# 
# 
# #### make p value density plot
# library(ggplot2)
# data <- read.xlsx("/Users/zhangyiming/Library/Mobile Documents/com~apple~CloudDocs/gene list/mouse_cell_altas.xlsx")
# 
# data$q <- p.adjust(data$p_val, method = "fdr")
# data$logq <- -log10(data$q)
# 
# data[is.infinite(data$logq), "logq"] = -1
# max(data$logq)
# 
# data[data$logq == -1, "logq"] <- 307
# 
# summary(data$logq>=2)
# dim(data)
# head(sort(table(data[data$logq>=2,"gene"]),decreasing = T))
# length(unique(data[,"gene"]))
# 
# 
# log10(0.01)
# 
# 
# p <- ggplot(data=data, aes(x=logq)) + 
#     geom_density() +
#     labs(title="P value distribution of MCA data", x="-log10(q value)") +
#     theme_bw() + # 先用全白的主题，在手动添加自己想要的内容
#     theme(
#         text=element_text(size=35),    # 全图的字号  
#         legend.title=element_blank(),    # 不要图注的题目
#         axis.line = element_line(colour = "black"),    # 黑色边界线
#         plot.title = element_text(lineheight=.8, face="bold"),   # 图题的大小
#         axis.text.x = element_text(angle = 45, hjust = 1),    # 图的x坐标的方向
#         legend.key.size = unit(3,"line")    # 调整boxplot的legend中box的大小
#     )
# p
# ggsave(filename = "/Users/zhangyiming/Library/Mobile Documents/com~apple~CloudDocs/lung_cancer_10X/results/MCA_p_value_distribution.pdf", dpi=300)
# 
# ??pheatmap


# Labyrinthine trophoblast (Placenta)
# Fibroblast (Adult-Muscle)
# Mouse embryonic fibroblast (Trophoblast-Stem-Cell)
# microglia (Retina)
# Amacrine 
library(MASS)
input_dir = "seurat/results/cluster/"
for (i in list.dirs(input_dir)) {
    if ( i != input_dir) {
        print(i)
        
        in_file = paste(i, "2_clustered.rds", sep = "/")
        out_file = paste(i, "expression.mtx", sep = "/")
        
        data <- readRDS(in_file)
        
        write.table(as.matrix(data@raw.data), file = out_file, sep = "\t", quote = F, row.names = T, col.names = T)
    }
    
}

