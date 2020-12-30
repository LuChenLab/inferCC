options(stringsAsFactors = F)
library(Seurat)
library(scatterpie)
library(dplyr)
library(ggthemes)
library(ggrastr)
library(ComplexHeatmap)
library(circlize)
library(wesanderson)
library(openxlsx)
library(stringr)
library(Mfuzz)
library(reticulate)
library(ggradar)
library(reshape2)
library(ggrepel)
library(ggpubr)
library(ggthemes)


use_python("program/pyenv/shims/python")


gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
}


stage_colors = c("I"="#65A9A3", "II"="#4A933E", "III"="#EC7A21", "IV"="#D73F47", "LUAD_Normal"="#FECC1B", "LUSC_Normal"="#778793")
disease_colors = c("LUAD" = "#0084D1", "LUAD_Normal"="#FECC1B", "Normal"="#73BDFF", "LUSC_Normal"="#778793", "LUSC"="#A0C807", "LELC"="#6A0019", "CPI"="#C5000B")
tissue_colors = c("Bronichi"="#FB290F", "Lymph node"="#488F17", "Tumor"="#FC8210")


cluster_cols = c(
    as.character(wes_palette("Zissou1")),
    as.character(wes_palette("Royal2")),
    as.character(wes_palette("Royal1")),
    as.character(wes_palette("Darjeeling2")),
    as.character(wes_palette("Darjeeling1"))
)

patient_colors = cluster_cols[1:length(unique(obj@meta.data$PatientID))]
names(patient_colors) <- unique(obj@meta.data$PatientID)

colors = c(stage_colors, disease_colors, tissue_colors, patient_colors)
labels_size = 12
title_size = 15


files <- list.files("03_each_cells/", pattern = "results.csv", full.names = T, recursive = T)
files <- files[str_detect(files, "mfuzz")]
files <- files[!str_detect(files, "(Alveolar|Basal)")]

registerDoMC(10)
foreach(f = files) %dopar% {
    input_dir = dirname(f)
    
    if (!file.exists(f)) {
        return ()
    }
    
    res <- read.csv(f, row.names = 1, stringsAsFactors = F)
    # res$Clt <- sapply(res$Clt, function(x) {
    #     clt = c(
    #         "1"="I",
    #         "2"="II",
    #         "3"="III"
    #     )
    #     clt[x]
    # })
    
    res$Clt = paste0("M.", res$Clt)
    
    
    meta = obj@meta.data[order(obj@meta.data$Stage, obj@meta.data$PatientID), ]
    
    col_fun = colorRamp2(c(-5, 0, 5), c("blue", "white", "red"))
    
    
    stage_colors <- c(
        "I"="#65A9A3", 
        "II"="#4A933E",
        "III"="#EC7A21",
        "IV"="#D73F47"
    )
    
    
    
    # top annotation
    ta = HeatmapAnnotation(
        Stage = meta$Stage,
        PatientID = meta$PatientID,
        col = list(
            Stage =  stage_colors,
            PatientID = patient_colors
        ),
        show_legend = T,
        show_annotation_name = T,
        annotation_legend_param = list(
            Stage=list(
                direction = "horizontal", 
                nrow = 1,
                labels_gp = gpar(fontsize = labels_size),
                title_gp = gpar(fontsize = title_size)
            ),
            PatientID=list(
                direction = "horizontal",
                nrow = 4,
                labels_gp = gpar(fontsize = labels_size),
                title_gp = gpar(fontsize = title_size)
            )
        )
    )
    
    
    module_colors = wes_palette("GrandBudapest2", length(unique(res$Clt)), type = "continuous")
    names(module_colors) <- unique(res$Clt)
    
    # left annotation
    la = rowAnnotation(
        Module = res$Clt,
        col = list(
            Module = module_colors
        ),
        show_legend = T,
        show_annotation_name = T,
        annotation_name_rot = 0,
        annotation_legend_param=list(
            Module=list(
                direction = "horizontal", 
                nrow = 1,
                labels_gp = gpar(fontsize = labels_size),
                title_gp = gpar(fontsize = title_size)
            )
        )
    )
    
    h <- Heatmap(
        name = "Expr",
        obj@scale.data[res$gene, rownames(meta)],
        cluster_rows = F,
        cluster_columns = F,
        show_column_names = F,
        show_row_names = F,
        col = col_fun,
        top_annotation = ta,
        left_annotation = la,
        border = T,
        column_title_gp = gpar(fontsize=0),
        row_title_gp = gpar(fontsize = 0),
        column_split = meta$Stage,
        row_split = res$Clt,
        heatmap_legend_param = list(
            direction = "horizontal",
            labels_gp = gpar(fontsize = labels_size),
            title_gp = gpar(fontsize = title_size)
        )
    )
    
    pdf(paste0(input_dir, "/mfuzz.pdf"), height = 6, width = 6)
    draw(
        h, 
        column_title = "stage module",
        column_title_gp = gpar(fontsize = 20),
        merge_legend = TRUE, 
        heatmap_legend_side = "bottom", 
        annotation_legend_side = "bottom"
    )
    
    dev.off()
}



format_radar_expr_mat <- function(
    obj, 
    module,
    legend.position="bottom",
    plot.title="",
    legend.title = "",
    axis.label.size = 8,
    batch=F,
    Stage = NULL
) {
    meta = as.data.frame(obj@meta.data)
    
    if (!is.null(Stage)) {
        meta = meta[meta$Stage == Stage, ]
    }
    
    if (batch) {
        meta = meta[meta$Batch == 3, ]
        meta$Group = meta$PatientID
    } else {
        meta = meta[meta$Batch != 3, ]
        meta$Group = meta$PatientID # paste0(meta$PatientID, "(", meta$Stage, ")")
    }
    
    
    cells = meta[, "Group", drop = F]
    cells$Cells = rownames(cells)
    
    mat = obj@scale.data[as.character(module$gene), ]
    
    mat = melt(as.matrix(mat))
    mat = merge(cells, mat, by.x = "Cells", by.y = "Var2")
    
    module$Clt = as.character(module$Clt)
    # module$Clt = sapply(module$Clt, function(x) { paste0("M.", x) })
    
    module = merge(module, mat, by.x = "gene", by.y = "Var1", all.x = T)
    
    
    module = module %>% group_by(Clt, Group, gene) %>%
        mutate(m = mean(value)) %>%
        dplyr::select(Clt, Group, m) %>%
        unique() %>%
        group_by(Clt, Group) %>%
        mutate(value = mean(m)) %>%
        dplyr::select(Clt, Group, value) %>%
        unique()
    
    # return(module)
    module = as.data.frame(module)
    
    # return(module)
    #print(module)
    # 
    # module
    p <- ggradar(
        dcast(module, Group~Clt),
        centre.y = floor(min(module$value)),
        values.radar = round(seq(from=floor(min(module$value)), to=max(module$value), length.out = 3), 2),
        grid.max = ceiling(max(module$value)),
        grid.label.size = 8,
        legend.position = legend.position,
        legend.title = legend.title,
        plot.title = plot.title,
        axis.label.size = axis.label.size,
        group.line.width = 0.5,
        group.point.size = 4,
        axis.label.offset = 1.15,
        base.size = 10
    )
    p = p + theme(
        legend.direction="horizontal",
        legend.spacing.x = unit(0.2, 'cm'),
        legend.position = c(0.5, 0.01),
        legend.background = element_blank(),
        legend.text = element_text(size = 20),
        plot.title = element_text(size = 25, hjust = 0.5)
    ) +
        guides(
            color = guide_legend(override.aes = list(ncol = 3))
        )
    return(p)
}




for (i in unique(obj@meta.data$Stage)) {
    p = format_radar_expr_mat(
        obj, 
        res,
        legend.position="bottom",
        # plot.title=basename(dirname(dirname(i))),
        plot.title = paste0("Basal ", "(", i, ")"),
        axis.label.size = 10,
        Stage=i
    )
    
    ggsave(
        filename = here(paste0("mfuzz/mfuzz_radar_", i, ".pdf")),
        plot = p,
        width = 6,
        height = 6
    )
}



make_heatmap <- function(obj, res) {
    res$Clt = paste0("M.", res$Clt)
    
    meta = obj@meta.data[order(obj@meta.data$Stage, obj@meta.data$PatientID), ]
    col_fun = colorRamp2(c(-5, 0, 5), c("blue", "white", "red"))
    
    stage_colors <- c(
        "I"="#65A9A3", 
        "II"="#4A933E",
        "III"="#EC7A21",
        "IV"="#D73F47"
    )
    
    # top annotation
    ta = HeatmapAnnotation(
        Stage = meta$Stage,
        PatientID = meta$PatientID,
        col = list(
            Stage =  stage_colors,
            PatientID = patient_colors
        ),
        show_legend = T,
        show_annotation_name = T,
        annotation_legend_param = list(
            Stage=list(
                direction = "horizontal", 
                nrow = 1,
                labels_gp = gpar(fontsize = labels_size),
                title_gp = gpar(fontsize = title_size)
            ),
            PatientID=list(
                direction = "horizontal",
                nrow = 4,
                labels_gp = gpar(fontsize = labels_size),
                title_gp = gpar(fontsize = title_size)
            )
        )
    )
    
    
    module_colors = wes_palette("GrandBudapest2", length(unique(res$Clt)), type = "continuous")
    names(module_colors) <- unique(res$Clt)
    
    # left annotation
    la = rowAnnotation(
        Module = res$Clt,
        col = list(
            Module = module_colors
        ),
        show_legend = T,
        show_annotation_name = T,
        annotation_name_rot = 0,
        annotation_legend_param=list(
            Module=list(
                direction = "horizontal", 
                nrow = 1,
                labels_gp = gpar(fontsize = labels_size),
                title_gp = gpar(fontsize = title_size)
            )
        )
    )
    
    h <- Heatmap(
        name = "Expr",
        obj@scale.data[res$gene, rownames(meta)],
        cluster_rows = F,
        cluster_columns = F,
        show_column_names = F,
        show_row_names = F,
        col = col_fun,
        top_annotation = ta,
        left_annotation = la,
        border = T,
        column_title_gp = gpar(fontsize=0),
        row_title_gp = gpar(fontsize = 0),
        column_split = meta$Stage,
        row_split = res$Clt,
        heatmap_legend_param = list(
            direction = "horizontal",
            labels_gp = gpar(fontsize = labels_size),
            title_gp = gpar(fontsize = title_size)
        )
    )
    
    draw(
        h, 
        column_title = "stage module",
        column_title_gp = gpar(fontsize = 20),
        merge_legend = TRUE, 
        heatmap_legend_side = "bottom", 
        annotation_legend_side = "bottom"
    )
}


expr <- readRDS("02_rds/seurat_obj_nm.rds")


meta <- read.xlsx("02_rds/NM_meta.xlsx")

rownames(meta) <- meta$SampleID
meta$PatientID <- paste0("DL", meta$PatientID)

expr@meta.data$PatientID <- as.character(expr@meta.data$PatientID)
expr@meta.data$Disease <- as.character(expr@meta.data$Disease)
expr@meta.data$PatientID[expr@meta.data$Batch == 3] <- meta[as.character(expr@meta.data$SampleID[expr@meta.data$Batch == 3]), "PatientID"]
expr@meta.data$Disease[expr@meta.data$Batch == 3] <- meta[as.character(expr@meta.data$SampleID[expr@meta.data$Batch == 3]), "Disease"]
expr@meta.data$Stage[expr@meta.data$Batch == 3] <- meta[as.character(expr@meta.data$SampleID[expr@meta.data$Batch == 3]), "Stage"]



for(i in c("LUAD", "LUSC")) {
    for (j in unique(expr@meta.data$cell_name)) {
        tryCatch({
            temp_meta <- expr@meta.data[expr@meta.data$cell_name == j & expr@meta.data$Disease == i, ]
            
            temp_obj <- CreateSeuratObject(
                expr@raw.data[, rownames(temp_meta)],
                meta.data = temp_meta
            )
            
            temp_obj@scale.data = expr@scale.data[, rownames(temp_meta)]
            
            # temp_obj@meta.data$PatientID <- paste(temp_obj@meta.data$PatientID, temp_obj@meta.data$Stage, sep = "-")
            
            j <- str_replace(j, "\\s+", "_")
            print(paste(i, j))
            
            csv_file <- paste0("03_each_cells/", i, "/", j, "/mfuzz/results.csv")
            
            if (file.exists(csv_file)) {
                res <- read.csv(csv_file, row.names = 1, stringsAsFactors = F)
                # res$Clt <- paste0("M.", res$Clt)
                
                pdf(paste0("mfuzz/DL/", i, "_", j, ".pdf"), width = 6, height = 6)
                make_heatmap(temp_obj, res)
                dev.off()
            
                p <- format_radar_expr_mat(
                    temp_obj,
                    res,
                    legend.position="bottom",
                    # plot.title=basename(dirname(dirname(i))),
                    plot.title = paste(j, "(DL2018)"),
                    axis.label.size = 10,
                    batch = T
                )

                ggsave(
                    filename = paste0("mfuzz/DL/", i, "_", j, "_radar.pdf"),
                    plot = p,
                    width = 6,
                    height = 6
                )
            }
        }, error=function(e) {
            next()
        }, warning=function(e) {
            
        }, finally = function() {
            next
        })
    }
    
}


temp_meta <- expr@meta.data[expr@meta.data$cell_name == "Alveolar II" & expr@meta.data$Disease == "LUAD", ]

temp_obj <- CreateSeuratObject(
    expr@raw.data[, rownames(temp_meta)],
    meta.data = temp_meta
)

temp_obj@scale.data = expr@scale.data[, rownames(temp_meta)]

temp_obj@meta.data$PatientID <- paste(temp_obj@meta.data$PatientID, temp_obj@meta.data$Stage, sep = "-")

p <- format_radar_expr_mat(
    temp_obj, 
    res,
    legend.position="bottom",
    # plot.title=basename(dirname(dirname(i))),
    plot.title = "ATII (DL2018)",
    axis.label.size = 10,
    batch = T
)

p

ggsave(
    filename = paste0("mfuzz/mfuzz_radar_DL.pdf"),
    plot = p,
    width = 6,
    height = 6
)