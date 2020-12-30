
library(ComplexHeatmap)
library(doMC)
library(stringr)

files = list.files("LungCancer10x/03_each_cells/", pattern = "seurat.rds", full.names = T, recursive = T)


stage_colors = c(
    "I"="#65A9A3", 
    "II"="#4A933E", 
    "III"="#EC7A21", 
    "IV"="#D73F47", 
    "LUAD_Normal"="#FECC1B", 
    "LUSC_Normal"="#778793"
)

disease_colors = c(
    "LUAD" = "#0084D1",
    "LUAD_Normal"="#FECC1B", 
    "Normal"="#73BDFF",
    "LUSC_Normal"="#778793", 
    "LUSC"="#A0C807", 
    "LELC"="#6A0019",
    "CPI"="#C5000B"
)


registerDoMC(10)
foreach(i = files, .errorhandling = "pass") %dopar% {
    
    cell <- basename(dirname(i))
    type <- basename(dirname(dirname(i)))
    
    obj <- readRDS(i)
    temp <- melt(as.matrix(obj@scale.data[res$gene, ]))
    
    output_dir <- paste0("LungCancer10x/corr_gene/", type)
    dir.create(output_dir, showWarnings = F, recursive = T)

    temp$PatientID <- obj@meta.data[temp$Var2, "PatientID"]

    temp <- temp %>%
        group_by(Var1, PatientID) %>%
        mutate(value = mean(value)) %>%
        dplyr::select(Var1, value, PatientID) %>%
        unique()

    temp <- dcast(temp, Var1 ~ PatientID, fun.aggregate = mean, value.var = "value")

    rownames(temp) <- temp$Var1
    temp <- temp[, colnames(temp) != "Var1"]


    temp_cor <- cor(temp)

    temp_meta <- unique(obj@meta.data[, c("PatientID", "Stage")])
    rownames(temp_meta) <- temp_meta$PatientID

    ba <- HeatmapAnnotation(
        Stage=temp_meta[colnames(temp_cor), "Stage"],
        col = list(
            Stage=stage_colors
        )
    )

    ra <- rowAnnotation(
        Stage=temp_meta[rownames(temp_cor), "Stage"],
        col = list(
            Stage=stage_colors
        ),
        show_legend = F,
        show_annotation_name = F
    )

    h <- Heatmap(
        temp_cor,
        name = "Corr",
        right_annotation = ra,
        bottom_annotation = ba
    )

    pdf(paste0(output_dir, "/", cell, ".pdf"), width = 8, height = 6)
    draw(h, merge_legend = T)
    dev.off()
}



foreach(i = files[str_detect(files, "total")], .errorhandling = "pass") %dopar% {
    
    cell <- basename(dirname(i))
    type <- basename(dirname(dirname(i)))
    
    
    obj <- readRDS(i)
    temp <- melt(as.matrix(obj@scale.data[res$gene, ]))
    
    output_dir <- paste0("LungCancer10x/corr_gene/", type)
    dir.create(output_dir, showWarnings = F, recursive = T)
    
    temp$PatientID <- obj@meta.data[temp$Var2, "PatientID"]
    temp$SampleID <- obj@meta.data[temp$Var2, "SampleID"]
    temp$Disease <- obj@meta.data[temp$Var2, "Disease"]
    
    temp <- temp %>%
        group_by(Var1, SampleID) %>%
        mutate(value = mean(value)) %>%
        dplyr::select(Var1, value, SampleID, PatientID, Disease) %>%
        unique()
    
    temp <- dcast(temp, Var1 ~ SampleID, fun.aggregate = mean, value.var = "value")
    
    rownames(temp) <- temp$Var1
    temp <- temp[, colnames(temp) != "Var1"]
    
    
    temp_cor <- cor(temp)
    
    temp_meta <- unique(obj@meta.data[, c("PatientID", "Stage", "SampleID", "Disease")])
    rownames(temp_meta) <- temp_meta$SampleID
    
    ba <- HeatmapAnnotation(
        Stage=temp_meta[colnames(temp_cor), "Stage"],
        Disease=temp_meta[colnames(temp_cor), "Disease"], 
        col = list(
            Stage=stage_colors,
            Disease=disease_colors
        )
    )
    
    ra <- rowAnnotation(
        Stage=temp_meta[rownames(temp_cor), "Stage"],
        Disease=temp_meta[rownames(temp_cor), "Disease"], 
        col = list(
            Stage=stage_colors,
            Disease=disease_colors
        ),
        show_legend = F,
        show_annotation_name = F
    )
    
    h <- Heatmap(
        temp_cor,
        name = "Corr",
        right_annotation = ra,
        bottom_annotation = ba
    )
    
    pdf(paste0(output_dir, "/", cell, ".pdf"), width = 8, height = 6)
    draw(h, merge_legend = T)
    dev.off()
}