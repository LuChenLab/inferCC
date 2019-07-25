
immu = c(
    "B cells",
    "CD4",
    "CD8",
    "Dendritic",
    "Exhaust T",
    "Endothelial",
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


sink("outfile.txt")
cat("hello")
cat("\n")
cat("world")
sink()


outfile = "02_rds/expr.csv"
cols = rownames(obj@meta.data)
rows = rownames(obj@raw.data)

cat(",", file = outfile, sep = ",")
cat(cols, file = outfile, sep = ",", append=T)
cat("\n", file = outfile, sep = "", append=T)

for(i in 1:length(rows)) {
    cat(paste0(rows[i], ","), file = outfile, sep = ",", append = T)
    cat(as.vector(obj@raw.data[i, cols]), file = outfile, sep = ",", append=TRUE)
    cat("\n", file = outfile, sep = "", append=T) 
}


cells = rownames(obj@meta.data)[as.character(obj@meta.data$cell_name) %in% immu & obj@meta.data$Batch != 3]

# expr <- readRDS("all_cell_expr.rds")
new_obj <- CreateSeuratObject(
    expr[,cells],
    meta.data = obj@meta.data[cells, ]
)


saveRDS(new_obj, "seurat_obj_immu.rds")

expr = readRDS("all_cell_expr.rds")

cells = rownames(obj@meta.data)[obj@meta.data$Batch != 3]

# expr <- readRDS("all_cell_expr.rds")
new_obj <- CreateSeuratObject(
    expr[, cells],
    meta.data = obj@meta.data[cells, ]
)

saveRDS(new_obj, "seurat_obj_hx.rds")


## Immu LUAD
cells = rownames(obj@meta.data)[as.character(obj@meta.data$cell_name) %in% immu & obj@meta.data$Batch != 3 & obj@meta.data$Disease == "LUAD"]

# expr <- readRDS("all_cell_expr.rds")
new_obj <- CreateSeuratObject(
    expr[, cells],
    meta.data = obj@meta.data[cells, ]
)
saveRDS(new_obj, "seurat_obj_immu_luad.rds")


## Immu LUSC
cells = rownames(obj@meta.data)[as.character(obj@meta.data$cell_name) %in% immu & obj@meta.data$Batch != 3 & obj@meta.data$Disease == "LUSC"]

# expr <- readRDS("all_cell_expr.rds")
new_obj <- CreateSeuratObject(
    expr[, cells],
    meta.data = obj@meta.data[cells, ]
)
saveRDS(new_obj, "seurat_obj_immu_lusc.rds")



for(i in 1:3) {
    set.seed(i)

    cells = sample(1:nrow(obj@meta.data), 100000)

    new_obj <- CreateSeuratObject(
        obj@raw.data[, cells],
        meta.data = obj@meta.data[cells, ]
    )

    saveRDS(new_obj, paste0("seurat_obj_ori_R", i, ".rds")) 
}


for(i in unique(obj@meta.data$Disease)) {
    cells = rownames(obj@meta.data)[obj@meta.data$Disease == i]
    new_obj <- CreateSeuratObject(
        obj@raw.data[, cells],
        meta.data=obj@meta.data[cells, ]
    )

    saveRDS(new_obj, paste0("seurat_obj_immu_", i, ".rds"))
}