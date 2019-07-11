
immu = c(
    "B cells",
    "CD4+",
    "Dendritic",
    "Exhaust T",
    "Erythroid precursor",
    "Granulocyte",
    "Mast",
    "Monocytes",
    "NK",
    "T_cells",
    "Treg"
)

library(Seurat)

obj <- readRDS("seurat_obj_ori.rds")


cells = rownames(obj@meta.data)[as.character(obj@meta.data$cell_name) %in% immu & obj@meta.data$Batch != 3]

expr <- readRDS("all_cell_expr.rds")
new_obj <- CreateSeuratObject(
    expr[, cells],
    meta.data = obj@meta.data[cells, ]
)


saveRDS(new_obj, "seurat_obj_ori_immu.rds")