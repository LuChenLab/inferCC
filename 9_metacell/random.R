obj <- readRDS("20190610_seurat_obj_selected_patients_harmony_formatted_meta.rds")

cells = rownames(obj@meta.data)[obj@meta.data$Batch != 3]


set.seed(1)
cells1 = cells[sample(1:length(cells), 100000)]

set.seed(2)
cells2 = cells[sample(1:length(cells), 100000)]


set.seed(3)
cells3 = cells[sample(1:length(cells), 100000)]


set.seed(4)
cells4 = cells[sample(1:length(cells), 100000)]


set.seed(5)
cells5 = cells[sample(1:length(cells), 100000)]


new_obj <- CreateSeuratObject(
    obj@raw.data[, cells1],
    meta.data = obj@meta.data[cells1, ]
)
saveRDS(
    new_obj, 
    "/mnt/raid62/Lung_cancer_10x/00_data_ingest/04_rds_generated/metacell_random_selected_1.rds"
)


new_obj <- CreateSeuratObject(
    obj@raw.data[, cells2],
    meta.data = obj@meta.data[cells2, ]
)
saveRDS(
    new_obj, 
    "/mnt/raid62/Lung_cancer_10x/00_data_ingest/04_rds_generated/metacell_random_selected_2.rds"
)


new_obj <- CreateSeuratObject(
    obj@raw.data[, cells3],
    meta.data = obj@meta.data[cells3, ]
)
saveRDS(
    new_obj, 
    "/mnt/raid62/Lung_cancer_10x/00_data_ingest/04_rds_generated/metacell_random_selected_3.rds"
)


new_obj <- CreateSeuratObject(
    obj@raw.data[, cells4],
    meta.data = obj@meta.data[cells4, ]
)
saveRDS(
    new_obj, 
    "/mnt/raid62/Lung_cancer_10x/00_data_ingest/04_rds_generated/metacell_random_selected_4.rds"
)



new_obj <- CreateSeuratObject(
    obj@raw.data[, cells5],
    meta.data = obj@meta.data[cells5, ]
)
saveRDS(
    new_obj, 
    "/mnt/raid62/Lung_cancer_10x/00_data_ingest/04_rds_generated/metacell_random_selected_5.rds"
)

