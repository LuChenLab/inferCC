### Created by Zhang at 2019.01.22
### This scripts is used to  run multiple CCA by split it into multiple groups
###

library(openxlsx)
library("Seurat")
library("logging")

basicConfig()


setwd(output)

## laod meta info
meta <- meta[order(meta$ID),]

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
            group = meta[i, "tissue_type"]
            )
        
        idx = idx + 1
        
    }
    
    return (objects)
}



### Collect genes.use
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

run_CCA <- function(objects, ccs=15) {
    genes.use <- find_genes_use(objects)
    
    cell_ids = c()
    for(i in objects) {
        cell_ids = c(cell_ids, i@meta.data$sample_id)
    }
    
    obj <- RunMultiCCA(
        objects, 
        genes.use = genes.use, 
        num.ccs = ccs,
        add.cell.ids = unique(cell_ids)
        )
    
    return(obj)
}


objects <- load_by_meta_info(meta, input)

saveRDS(objects, "objects_between_groups.rds")

objects <- run_CCA(objects)

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

saveRDS(obj.list, file = paste(output, "objects_list.rds", sep = "/"))


## run tSNE each groups
for ( i in data) {
    print(unique(i@meta.data$sample_id))
}

cell_ids <- c(
    "Tuberculosis",
    "Normal",
    "Squamous_cell_carcinoma",
    "Pseudotumor",
    "Tumor",
    "Adenocarcinoma"
)

obj_cca <- run_CCA(data, cell_ids)

saveRDS(obj_cca, "cca_each_group.rds")

obj <- AlignSubspace(obj, reduction.type = "cca", grouping.var = "group", dims.align = 1:15)
obj <- RunTSNE(obj, reduction.use = "cca.aligned", dims.use = 1:15, do.fast = T)
obj <- FindClusters(obj, reduction.type = "cca.aligned", resolution = 0.6, dims.use = 1:15)

saveRDS(obj, "clustered_each_group.rds")


pdf("tSEN_each_group.pdf", height = 20, width = 15)
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
    11, 
    10, 
    13, 
    0, 2,
    12, 4,
    1, 7, 8, 9,
    15,
    2, 14, 3
)
cell.ids <- c(
    "Basal",
    "Club",
    "B", 
    "NK", "NK",
    "T", "T", 
    "Mono", "Mono", "Mono", "Mono",
    "Fiblast",
    "Tc", "Tc", "Tc"
)
obj1 <- obj
obj1@ident <- plyr::mapvalues(x = obj1@ident, from = current.clusters, to = cell.ids)

png("tSNE_anno.png",width = 10,height = 15, res = 300)
p1 <- TSNEPlot(obj, do.return = T, pt.size = 0.5, group.by = "group")
p2 <- TSNEPlot(object = obj1, do.return=T, do.label = TRUE, label.size=5, pt.size = 0.5)
grid.arrange(p1, p2, ncol=1)
dev.off()


pdf("heatmap_top10.pdf", height = 10, width = 20)
DoHeatmap(object = obj, genes.use = All.top10.markers$gene, slim.col.label = TRUE, remove.key = TRUE)
dev.off()

# genes.use = c()
# for (i in 1:length(ob.list)) {
#     genes.use = c(genes.use, head(rownames(ob.list[[i]]@hvg.info), 2000))
# }
# genes.use = names(which(table(genes.use) > 1))
# for (i in 1:length(ob.list)) {
#     genes.use = genes.use[genes.use %in% rownames(ob.list[[i]]@scale.data)]
# }
# length(genes.use) #1958 genes