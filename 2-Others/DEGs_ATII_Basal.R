options(stringsAsFactors = F)
library(openxlsx)
library(Seurat)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggcorrplot)
library(ComplexHeatmap)
library(circlize)
library(stringr)
library(wesanderson)
library(dplyr)
library(reshape2)
library(doMC)


meta <- read.csv("02_rds/meta_after_singleR.csv", row.names = 1, stringsAsFactors = F)
nm_meta <- read.xlsx("02_rds/NM_meta.xlsx")
rownames(nm_meta) <- nm_meta$SampleID

meta[meta$Batch == 3, "Stage"] <- nm_meta[meta$SampleID[meta$Batch == 3], "Stage"]
meta[meta$Batch == 3, "Disease"] <- nm_meta[meta$SampleID[meta$Batch == 3], "Disease"]


cells <- rownames(meta)[meta$cell_short %in% c("ATII", "Basal")]

expr <- readRDS("02_rds/all_cell_expr.rds")


obj <- CreateSeuratObject(
    expr[, cells],
    meta.data = meta[cells, ]
)


obj <- NormalizeData(object = obj, normalization.method = "LogNormalize", scale.factor = 10000)
obj <- FindVariableGenes(object = obj, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
obj <- ScaleData(object = obj, vars.to.regress = c("nUMI", "percent.mito"))

saveRDS(obj, "03_each_cells/select_features/atii_basal.rds")

obj <- readRDS("03_each_cells/select_features/atii_basal.rds")


short_names <- c(
    "APII"="ATII",
    "Alveolar II"="ATII",
    "Basal"="Basal",
    "B cells"="B",
    "CD4+"="CD4",
    "CD4"="CD4",
    "CD8+"="CD8",
    "CD8"="CD8",
    "Ciliated"="Cilia",
    "Club"="Club",
    "DC"="DC",
    "Dendritic"="DC",
    "Endothelial"="EC",
    "Epithelial"="Epi",
    "Exhaust T"="Exh T",
    "Fibroblasts"="Fib",
    "Granulocyte"="Gran",
    "Mast"="Mast",
    "Mascrophages"="Mφ",
    "Macrophages"="Mφ",
    "Monocytes"="Mo",
    "Neuroendocrine"="NE",
    "NK"="NK",
    "Tregs"="Tregs"
)

short_names <- as.data.frame(short_names)


meta <- obj@meta.data
meta$Disease = as.character(meta$Disease)
meta$Disease <- sapply(meta$Disease, function(x) {
    disease = c(
        "LUAD"="AD",
        "LUSC"="SC",
        "LUAD_Normal"="NL(AD)",
        "LUSC_Normal"="NL(SC)",
        "LULC"="LC",
        "LULC_Normal"="NL(LC)",
        "Pleiomorphic"="UPS",
        "Pleiomorphic Normal"="NL(UPS)"
    )
    
    as.character(disease[x])
})

meta$cell_name <- short_names[meta$cell_name, 1]

expr <- melt(as.matrix(obj@scale.data))
expr$ident <- paste(
    meta[expr$Var2, "cell_name"], 
    meta[expr$Var2, "Disease"],
    meta[expr$Var2, "Stage"],
    sep = "-"
)


expr1 <- dcast(expr, Var1~ident, fun.aggregate = mean, value.var = "value")
rownames(expr1) <- expr1$Var1
expr1 <- expr1[, colnames(expr) != "Var1"]

colors = list(
    Stage = c(
        "I"="#65A9A3", 
        "II"="#4A933E", 
        "III"="#EC7A21", 
        "IV"="#D73F47", 
        "LUAD_Normal"="#FECC1B", 
        "LUSC_Normal"="#778793"
    ),
    Disease = c(
        "AD" = "#0084D1", 
        "NL(AD)"="#FECC1B", 
        "NL(AD"="#778793", 
        "Normal"="#73BDFF", 
        "NL(SC)"="#778793", 
        "SC"="#A0C807",
        "LC"="#91822B",
        "NL(LC)"="#EDCAB0",
        "UPS"="#EBAAA4",
        "NL(UPS)"="#B43018"
    ),
    Batch = c("1"="#E9B0B7", "2"="#90D0E0", "3"="#769982"),
    Cell = c(
        "ATII"="#FF0000",
        "Basal"="#F2AD00",
        "CD4"="#F98400",
        "CD8"="#5BBCD6",
        "Cilia"="#85D4E3",
        "Club"="#F4B5BD",
        "DC"="#9C964A",
        "EC"="#CDC08C",
        "Epi"="#FAD77B",
        "Fib"="#ECCBAE",
        "Mast"="#D69C4E",
        "Mo"="#ABDDDE",
        "NK"="#000000",
        "NE"="#F3DF6C",
        "Tregs"="#CEAB07",
        "Mφ"="#046C9A"
    )
)

cluster_cols = c(
    as.character(wes_palette("Zissou1")),
    as.character(wes_palette("Royal2")),
    as.character(wes_palette("Royal1")),
    as.character(wes_palette("Darjeeling2")),
    as.character(wes_palette("Darjeeling1")),
    wes_palette("BottleRocket2"),
    wes_palette("Rushmore1")
)


patient_colors = cluster_cols[1:length(unique(obj@meta.data$PatientID))]
names(patient_colors) <- unique(obj@meta.data$PatientID)


expr_cor <- cor(expr1)
ba <- HeatmapAnnotation(
    Disease = sapply(colnames(expr_cor), function(x) { str_split(x, "-")[[1]][2] }),
    Stage = sapply(colnames(expr_cor), function(x) { str_split(x, "-")[[1]][3] }),
    Cell=sapply(colnames(expr_cor), function(x) { str_split(x, "-")[[1]][1] }),
    col = colors
)


ra <- rowAnnotation(
    Disease = sapply(rownames(expr_cor), function(x) { str_split(x, "-")[[1]][2] }),
    Stage = sapply(rownames(expr_cor), function(x) { str_split(x, "-")[[1]][3] }),
    Cell=sapply(rownames(expr_cor), function(x) { str_split(x, "-")[[1]][1] }),
    col = colors,
    show_annotation_name = F,
    show_legend = F
)

scale_dot_size = 0.5

pdf("03_each_cells/select_features//ATII_Basal_gene_corr.pdf", width = 10, height = 6)
Heatmap(
    expr_cor,
    name = "Corr",
    show_column_names = F,
    bottom_annotation = ba,
    right_annotation = ra,
    rect_gp = gpar(type = "none"),
    cell_fun = function(j, i, x, y, width, height, fill) {
        s = min(unit.c(convertWidth(width, "cm"), convertHeight(height, "cm")))
        grid.circle(
            x = x,
            y = y, 
            r = s * (abs(expr_cor[i, j]) + 0.1) * scale_dot_size,
            gp = gpar(col = "white", fill = fill)
        )
    }
)
dev.off()


expr_cor <- cor(expr, method = "spearman")
ba <- HeatmapAnnotation(
    Disease = sapply(colnames(expr_cor), function(x) { str_split(x, "-")[[1]][2] }),
    Stage = sapply(colnames(expr_cor), function(x) { str_split(x, "-")[[1]][3] }),
    Cell=sapply(colnames(expr_cor), function(x) { str_split(x, "-")[[1]][1] }),
    col = colors
)


ra <- rowAnnotation(
    Disease = sapply(rownames(expr_cor), function(x) { str_split(x, "-")[[1]][2] }),
    Stage = sapply(rownames(expr_cor), function(x) { str_split(x, "-")[[1]][3] }),
    Cell=sapply(rownames(expr_cor), function(x) { str_split(x, "-")[[1]][1] }),
    col = colors,
    show_annotation_name = F,
    show_legend = F
)

pdf("03_each_cells/select_features//ATII_Basal_gene_corr_spearman.pdf", width = 10, height = 6)
Heatmap(
    expr_cor,
    name = "Corr",
    show_column_names = F,
    bottom_annotation = ba,
    right_annotation = ra,
    rect_gp = gpar(type = "none"),
    cell_fun = function(j, i, x, y, width, height, fill) {
        s = min(unit.c(convertWidth(width, "cm"), convertHeight(height, "cm")))
        grid.circle(
            x = x,
            y = y, 
            r = s * (abs(expr_cor[i, j]) + 0.1) * scale_dot_size,
            gp = gpar(col = "white", fill = fill)
        )
    }
)
dev.off()


### Markers between ATII and Basal
markers <- FindMarkers(
    obj, 
    ident.1 = rownames(obj@meta.data)[obj@meta.data$cell_name == "Alveolar II"], 
    ident.2 = rownames(obj@meta.data)[obj@meta.data$cell_name == "Basal"],
    logfc.threshold = 0
)

markers$gene <- rownames(markers)

write.csv(markers, "03_each_cells/select_features/atii_basal_markers.csv")

markers_de <- markers[abs(markers$avg_logFC) > 0.25 & markers$p_val_adj < 0.05, ]


eg <- bitr(markers_de$gene, fromType="SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)


go <- enrichGO(
    eg$ENTREZID,
    OrgDb = org.Hs.eg.db,
    ont = "BP"
)

gosim <- simplify(go)

write.csv(gosim, "03_each_cells/select_features/atii_basal_markers_go.csv")

target_genes <- c(
    "PDCD1", "CD28", "CTLA4",
    "ICOS", "BTLA", "LAG3",
    "TNFRSF9", "TNFRSF4",
    "CD27", "CD40LG", "HAVCR2",
    "GZMA", "GZMB", "GZMM", 
    "GZMH", "GZMK", "SIGLEC15"
)


pdf("03_each_cells/select_features/basal_target_genes_vln_cluster.pdf", width = 16, height = 20)
VlnPlot(basal, features.plot = target_genes, group.by = "res.0.8")
dev.off()

pdf("03_each_cells/select_features/basal_target_genes_vln_stage.pdf", width = 16, height = 20)
VlnPlot(basal, features.plot = target_genes, group.by = "Stage")
dev.off()


find_markers_of_disease <- function(obj, disease=TRUE, group.by = "Stage", n.core = 4) {
    rule = str_detect(obj@meta.data$Disease, "Normal")
    
    if (disease) {
        rule = !rule
    }
    meta = obj@meta.data[rule, ]
    
    registerDoMC(n.core)
    res = foreach(i = unique(meta$Stage), .combine = "rbind", .errorhandling = "pass") %dopar% {
        cells1 = rownames(meta)[meta$Stage == i]
        cells2 = rownames(meta)[meta$Stage != i]
        
        temp <- FindMarkers(obj, ident.1 = cells1, ident.2 = cells2, logfc.threshold = 0)
        temp$ident = i
        temp$gene = rownames(temp)
        temp
    } 
    
    res
}

markers_basal <- find_markers_of_disease(basal)
markers_atii <- find_markers_of_disease(atii)

write.csv(markers_basal, "03_each_cells/select_features/markers_basal_stage.csv")
write.csv(markers_atii, "03_each_cells/select_features/markers_atii_stage.csv")

markers_basal <- read.csv("03_each_cells/select_features/markers_basal_stage.csv", row.names = 1)
markers_atii <- read.csv("03_each_cells/select_features/markers_atii_stage.csv", row.names = 1)

####
make_stage_corr_heatmap <- function(obj, markers, logfc=0.25, pval=0.05, topn=50, absolute = FALSE, features = NULL) {
    
    if (is.null(features)) {
        markers <- markers[markers$avg_logFC > logfc & markers$p_val_adj < pval, ]
        
        if(absolute) {
            markers <- markers[order(abs(markers$avg_logFC), decreasing = T), ]
        } else {
            markers <- markers[order(markers$avg_logFC, decreasing = T), ]
        }
        
        markers <- markers %>%
            group_by(ident) %>%
            slice(1:topn) %>%
            as.data.frame()
        
        features - markers$gene
    } 
  
    data = melt(as.matrix(obj@scale.data[features, ]))
    
    meta = obj@meta.data
    meta$Disease = as.character(meta$Disease)
    meta$Disease <- sapply(meta$Disease, function(x) {
        disease = c(
            "LUAD"="AD",
            "LUSC"="SC",
            "LUAD_Normal"="NL(AD)",
            "LUSC_Normal"="NL(SC)",
            "LULC"="LC",
            "LULC_Normal"="NL(LC)",
            "Pleiomorphic"="UPS",
            "Pleiomorphic Normal"="NL(UPS)"
        )
        
        as.character(disease[x])
    })
    
    
    data$ident <- paste(meta[data$Var2, "PatientID"], meta[data$Var2, "Disease"], meta[data$Var2, "Stage"], sep = "-")
    
    data <- data %>%
        group_by(Var1, ident) %>%
        mutate(value = mean(value)) %>%
        dplyr::select(Var1, ident, value) %>%
        unique() %>%
        as.data.frame()
    
    data <- dcast(data, Var1~ident, fun.aggregate = mean, value.var = "value")
    rownames(data) <- data$Var1
    data <- data[, colnames(data) != "Var1"]
    
    data = cor(data, method = "spearman")
    
    colors = list(
        Stage = c(
            "I"="#65A9A3", 
            "II"="#4A933E", 
            "III"="#EC7A21", 
            "IV"="#D73F47", 
            "LUAD_Normal"="#FECC1B", 
            "LUSC_Normal"="#778793"
        ),
        Disease = c(
            "AD" = "#0084D1", 
            "NL(AD)"="#FECC1B", 
            "NL(AD"="#778793", 
            "Normal"="#73BDFF", 
            "NL(SC)"="#778793", 
            "SC"="#A0C807",
            "LC"="#91822B",
            "NL(LC)"="#EDCAB0",
            "UPS"="#EBAAA4",
            "NL(UPS)"="#B43018"
        ),
        Batch = c("1"="#E9B0B7", "2"="#90D0E0", "3"="#769982"),
        Cell = c(
            "ATII"="#FF0000",
            "Basal"="#F2AD00",
            "CD4"="#F98400",
            "CD8"="#5BBCD6",
            "Cilia"="#85D4E3",
            "Club"="#F4B5BD",
            "DC"="#9C964A",
            "EC"="#CDC08C",
            "Epi"="#FAD77B",
            "Fib"="#ECCBAE",
            "Mast"="#D69C4E",
            "Mo"="#ABDDDE",
            "NK"="#000000",
            "NE"="#F3DF6C",
            "Tregs"="#CEAB07",
            "Mφ"="#046C9A"
        )
    )
    
    print(sapply(rownames(data), function(x) { str_split(x, "-")[[1]][2] }))
    ba <- HeatmapAnnotation(
        Disease = sapply(rownames(data), function(x) { str_split(x, "-")[[1]][2] }),
        Stage = sapply(rownames(data), function(x) { str_split(x, "-")[[1]][3] }),
        col = colors
    )
    
    
    ra <- rowAnnotation(
        Disease = sapply(rownames(data), function(x) { str_split(x, "-")[[1]][2] }),
        Stage = sapply(rownames(data), function(x) { str_split(x, "-")[[1]][3] }),
        col = colors,
        show_annotation_name = F,
        show_legend = F
    )
    
    # rownames(data) <- sapply(rownames(data), function(x) {
    #     str_split(x, " ")[[1]][1]
    # })
    
    Heatmap(
        data,
        name = "Corr",
        col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
        show_column_names = F,
        bottom_annotation = ba,
        right_annotation = ra,
        clustering_method_rows = "ward.D",
        clustering_method_columns = "ward.D"
    )
    
}

pdf("03_each_cells/select_features/atii_top20_stage_markers.pdf", width = 8, height = 6)
make_stage_corr_heatmap(atii, markers_atii, logfc=0.25, topn=20)
dev.off()

pdf("03_each_cells/select_features/basal_top20_stage_markers.pdf", width = 8, height = 6)
make_stage_corr_heatmap(basal, markers_basal, logfc=0.25, topn=20)
dev.off()


### markers between stage i and its normal tissue
find_markers_disease_normal_by_stage <- function(obj, n.cores = 4) {
    meta = obj@meta.data
    
    registerDoMC(n.cores)
    res = foreach(i = unique(meta$Stage), .combine="rbind", .errorhandling = "pass") %dopar% {
        temp_meta = meta[meta$Stage == i, ]
        
        cells1 = rownames(temp_meta)[!str_detect(temp_meta$Disease, "Normal")]
        cells2 = rownames(temp_meta)[str_detect(temp_meta$Disease, "Normal")]
        
        if (length(cells1) >= 5 && length(cells2) >= 5) {
            temp <- FindMarkers(obj, ident.1 = cells1, ident.2 = cells2, logfc.threshold = 0)
            temp$ident = i
            temp$gene = rownames(temp)
        } else {
            temp = NULL
        }
        temp
    }
    res
}


markers_atii_disea_normal <- find_markers_disease_normal_by_stage(atii)
write.csv(markers_atii_disea_normal, "03_each_cells/select_features/markers_atii_disease_normal_by_stage.csv")


pdf("03_each_cells/select_features/atii_top20_disease_normal_by_stage_markers.pdf", width = 8, height = 6)
make_stage_corr_heatmap(atii, markers_atii_disea_normal, logfc=0.25, topn=40, absolute = T)
dev.off()


### top disease and normal markers by LUAD and LUSC

atii_lusc <- read.xlsx("03_each_cells/total/Alveolar_II/normal_cancer/LUSC.xlsx", rowNames = T)
atii_luad <- read.xlsx("03_each_cells/total/Alveolar_II/normal_cancer/LUAD.xlsx", rowNames = T)

atii_lusc$ident = "LUSC"
atii_luad$ident = "LUAD"

atii_lusc$gene <- rownames(atii_lusc)
atii_luad$gene <- rownames(atii_luad)

atii_disease_markers <- rbind(atii_lusc, atii_luad)


pdf("03_each_cells/select_features/atii_all_genes.pdf", width = 10, height = 8)
make_stage_corr_heatmap(atii, atii_disease_markers, absolute = T, features = rownames(atii@raw.data))
dev.off()


pdf("03_each_cells/select_features/atii_disease_normal_markers.pdf", width = 8, height = 6)
make_stage_corr_heatmap(atii, atii_disease_markers, absolute = T, features = atii_disease_markers$gene[atii_disease_markers$p_val_adj < 0.05 & abs(atii_disease_markers$avg_logFC) > 0.25])
dev.off()


### find genes that show stage and normal disease spec

markers_atii_stage <- markers_atii[markers_atii$avg_logFC > 0.25 & markers_atii$p_val_adj < 0.05, ]
markers_atii_normal <- markers_atii_disea_normal[markers_atii_disea_normal$avg_logFC < -0.25 & markers_atii_disea_normal$p_val_adj < 0.05, ]

unique_markers <- unique(c(markers_atii_stage$gene, markers_atii_normal$gene))

write.csv(as.matrix(atii@raw.data[unique_markers, ]), "03_each_cells/total/Alveolar_II/scanpy/raw_for_features.csv")
write.csv(atii@scale.data[unique_markers, ], "03_each_cells/total/Alveolar_II/scanpy/scale_for_features.csv")


importance <- read.csv("03_each_cells/select_features/atii_importance.csv", row.names = 1)

pdf("03_each_cells/select_features/atii_wo_low_importance_markers.pdf", width = 10, height = 8)
make_stage_corr_heatmap(atii, features = importance$feature[importance$importance > 0 & importance$cumulative_importance > 0.99], markers = NULL)
dev.off()


pdf("03_each_cells/select_features/atii_all_importance_markers.pdf", width = 10, height = 8)
make_stage_corr_heatmap(atii, features = importance$feature[importance$importance > 0], markers = NULL)
dev.off()


### ATII LUAD stage module and not expressed in LUSC or normal

atii_module <- read.csv("03_each_cells/LUAD/Alveolar_II/mfuzz/results.csv", row.names = 1, stringsAsFactors = F)

atii_module$Clt = sapply(atii_module$Clt, function(x) {
    stages = c("1"="III", "2"="IV", "3"="I")
    
    paste0("M.", stages[as.character(x)])
})



make_stage_module_heatmap <- function(obj, module, col_fun = NULL, labels_size=15, title_size = 20, cluster=F, do.return=F, show_row_names = F) {
    if (is.null(col_fun)) {
        col_fun = colorRamp2(c(-5, 0, 5), c("blue", "white", "red"))
    } 
    
    meta <- obj@meta.data
    
    meta$Disease <- as.character(meta$Disease)
    meta$Disease[meta$Disease == "LUAD"] <- "AD"
    meta$Disease[meta$Disease == "LUSC"] <- "SC"
    meta$Disease[meta$Disease == "LUAD_Normal"] <- "NL(AD)"
    meta$Disease[meta$Disease == "LUSC_Normal"] <- "NL(SC)"
    meta$Disease[meta$Disease == "Pleiomorphic Normal"] <- "NL(UPS)"
    meta$Disease[meta$Disease == "Pleiomorphic"] <- "UPS"
    meta$Disease[meta$Disease == "LULC"] <- "LC"
    meta$Disease[meta$Disease == "LULC_Normal"] <- "NL(LC)"
    
    meta <- meta[order(meta$Disease, meta$Stage, meta$PatientID), ]

    patient_colors = cluster_cols[1:length(unique(meta$PatientID))]
    names(patient_colors) <- unique(meta$PatientID)
    
    # top annotation
    ta = HeatmapAnnotation(
        Stage = meta$Stage,
        Disease = meta$Disease,
        PatientID = meta$PatientID,
        col = list(
            Stage =  colors[["Stage"]],
            PatientID = patient_colors,
            Disease = colors[["Disease"]]
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
            ),
            Disease=list(
                direction = "horizontal",
                nrow = 1,
                labels_gp = gpar(fontsize = labels_size),
                title_gp = gpar(fontsize = title_size)
            )
        )
    )
    
    if (!cluster) {
        module_colors = as.character(wes_palette("GrandBudapest2", length(unique(module$Clt)), type = "continuous"))
        names(module_colors) <- unique(module$Clt)
        # left annotation
        la = rowAnnotation(
            Module = module$Clt,
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
            obj@scale.data[module$gene, rownames(meta)],
            cluster_rows = F,
            cluster_columns = F,
            show_column_names = F,
            show_row_names = show_row_names,
            col = col_fun,
            top_annotation = ta,
            left_annotation = la,
            border = T,
            column_title_gp = gpar(fontsize=0),
            row_title_gp = gpar(fontsize = 0),
            column_split = meta$Stage,
            row_split = module$Clt,
            heatmap_legend_param = list(
                direction = "horizontal",
                labels_gp = gpar(fontsize = labels_size),
                title_gp = gpar(fontsize = title_size)
            )
        )
    } else {
        h <- Heatmap(
            name = "Expr",
            obj@scale.data[module, rownames(meta)],
            cluster_rows = T,
            cluster_columns = F,
            show_column_names = F,
            show_row_names = show_row_names,
            col = col_fun,
            top_annotation = ta,
            column_split = meta$Disease,
            border = T,
            column_title_gp = gpar(fontsize=0),
            row_title_gp = gpar(fontsize = 0),
            heatmap_legend_param = list(
                direction = "horizontal",
                labels_gp = gpar(fontsize = labels_size),
                title_gp = gpar(fontsize = title_size)
            )
        )
    }

    if(do.return) {
        return(h)
    } else {
        draw(
            h, 
            column_title_gp = gpar(fontsize = 20),
            merge_legend = TRUE, 
            heatmap_legend_side = "bottom", 
            annotation_legend_side = "bottom"
        )
    }
}

pdf("03_each_cells/select_features/ATII_stage_module_scale_heatmap.pdf", width = 10, height = 8)
make_stage_module_heatmap(atii, atii_module)
dev.off()


pdf("03_each_cells/select_features/ATII_all_importance_scale_heatmap.pdf", width = 10, height = 8)
h <- make_stage_module_heatmap(atii, importance$feature[importance$importance > 0], cluster=T)
dev.off()

pdf("03_each_cells/select_features/ATII_wo_low_importance_scale_heatmap.pdf", width = 10, height = 30)
make_stage_module_heatmap(atii, importance$feature[importance$importance > 0 & importance$cumulative_importance > 0.99], cluster=T, show_row_names = T)
dev.off()

saveRDS(h, "03_each_cells/select_features/ATII_importance_gene_heatmap.rds")


### custome ATII importance gene expression 
h <- readRDS("03_each_cells/select_features/ATII_importance_gene_heatmap.rds")

meta <- atii@meta.data

meta$Disease <- as.character(meta$Disease)
meta$Disease[meta$Disease == "LUAD"] <- "AD"
meta$Disease[meta$Disease == "LUSC"] <- "SC"
meta$Disease[meta$Disease == "LUAD_Normal"] <- "NL(AD)"
meta$Disease[meta$Disease == "LUSC_Normal"] <- "NL(SC)"

meta <- meta[order(meta$Disease, meta$Stage, meta$PatientID), ]

labels_size = 15
title_size = 20

temp_gene = h@ht_list$Expr@row_names_param$labels[h@ht_list$Expr@row_order]

# top annotation
ta = HeatmapAnnotation(
    Stage = meta$Stage,
    Disease = meta$Disease,
    PatientID = meta$PatientID,
    col = list(
        Stage =  colors[["Stage"]],
        PatientID = patient_colors,
        Disease = colors[["Disease"]]
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
        ),
        Disease=list(
            direction = "horizontal",
            nrow = 1,
            labels_gp = gpar(fontsize = labels_size),
            title_gp = gpar(fontsize = title_size)
        )
    )
)

temp_gene_mark <- round(length(temp_gene) / 2)
temp_gene_mark <- temp_gene[(temp_gene_mark - 120): (temp_gene_mark + 30)]

temp_gene_mark <- c("ACSL4", "ALPL", "AQP4", "BMP2", "SFTPC", "C4BPA")

ra = rowAnnotation(
    foo = anno_mark(
        at = which(temp_gene %in% temp_gene_mark), 
        labels = temp_gene[temp_gene %in% temp_gene_mark],
        which = "row",
        side = "left"
    )
)

library(edgeR)
library(hett)

pdf("03_each_cells/select_features/ATII_all_importance_scale_heatmap.pdf", width = 10, height = 8)
temp = atii@scale.data[temp_gene, ]
draw(
    Heatmap(
        name = "Expr",
        as.matrix(temp),
        cluster_rows = F,
        cluster_columns = F,
        show_column_names = F,
        show_row_names = F,
        col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
        top_annotation = ta,
        left_annotation = ra,
        column_split = meta$Disease,
        border = T,
        column_title_gp = gpar(fontsize=0),
        row_title_gp = gpar(fontsize = 0),
        heatmap_legend_param = list(
            direction = "horizontal",
            labels_gp = gpar(fontsize = 15),
            title_gp = gpar(fontsize = 20)
        )
    ),
    column_title_gp = gpar(fontsize = 20),
    merge_legend = TRUE, 
    heatmap_legend_side = "bottom", 
    annotation_legend_side = "bottom"
)
dev.off()


expr <- readRDS("02_rds/seurat_obj_nm.rds")
meta <- read.xlsx("02_rds/NM_meta.xlsx")

rownames(meta) <- meta$SampleID
meta$PatientID <- paste0("DL", meta$PatientID)

expr@meta.data$PatientID <- as.character(expr@meta.data$PatientID)
expr@meta.data$Disease <- as.character(expr@meta.data$Disease)
expr@meta.data$PatientID[expr@meta.data$Batch == 3] <- meta[as.character(expr@meta.data$SampleID[expr@meta.data$Batch == 3]), "PatientID"]
expr@meta.data$Disease[expr@meta.data$Batch == 3] <- meta[as.character(expr@meta.data$SampleID[expr@meta.data$Batch == 3]), "Disease"]
expr@meta.data$Stage[expr@meta.data$Batch == 3] <- meta[as.character(expr@meta.data$SampleID[expr@meta.data$Batch == 3]), "Stage"]


temp_meta <- expr@meta.data[expr@meta.data$cell_name == "Alveolar II", ]

temp_obj <- CreateSeuratObject(
    expr@raw.data[, rownames(temp_meta)],
    meta.data = temp_meta
)

temp_obj@scale.data = expr@scale.data[, rownames(temp_meta)]



pdf("03_each_cells/select_features/ATII_all_importance_scale_heatmap_DL.pdf", width = 10, height = 8)
temp_gene_mark <- round(length(temp_gene) / 2)
make_stage_module_heatmap(obj, temp_gene[(temp_gene_mark - 120): (temp_gene_mark + 30)], cluster=T)
dev.off()

### make basal
importance <- read.csv("03_each_cells/select_features/basal_importance.csv", row.names = 1, stringsAsFactors = F)

make_stage_module_heatmap(basal, importance$feature[importance$importance > 0], cluster=T)



#####
obj <- readRDS("02_rds/seurat_obj_multi_reduction.rds")

meta <- read.csv("02_rds/meta_after_singleR.csv", row.names = 1, stringsAsFactors = F)
nm_meta <- read.xlsx("02_rds/NM_meta.xlsx")
rownames(nm_meta) <- nm_meta$SampleID

meta[meta$Batch == 3, "Stage"] <- nm_meta[meta$SampleID[meta$Batch == 3], "Stage"]
meta[meta$Batch == 3, "Disease"] <- nm_meta[meta$SampleID[meta$Batch == 3], "Disease"]

# meta <- unique(meta[, c("SampleID", "Stage", "Disease", "PatientID")])
# rownames(meta) <- meta$SampleID

meta$Disease = as.character(meta$Disease)
meta$Disease <- sapply(meta$Disease, function(x) {
    disease = c(
        "LUAD"="AD",
        "LUSC"="SC",
        "LUAD_Normal"="NL(AD)",
        "LUSC_Normal"="NL(SC)",
        "LULC"="LC",
        "LULC_Normal"="NL(LC)",
        "Pleiomorphic"="PL",
        "Pleiomorphic Normal"="NL(PL)"
    )
    
    as.character(disease[x])
})


set.seed(42)

meta <- meta[sample(1:nrow(meta), size = 10000), ]


data <- obj@scale.data[, rownames(meta)]

data <- melt(as.matrix(data))
data$ident <- paste(meta[data$Var2, "cell_name"], meta[data$Var2, "Disease"], meta[data$Var2, "Stage"])



data <- data %>%
    group_by(Var1, ident) %>%
    mutate(value = mean(value)) %>%
    dplyr::select(Var1, ident, value) %>%
    unique() %>%
    as.data.frame()

data <- dcast(data, Var1~ident, fun.aggregate = mean, value.var = "value")
rownames(data) <- data$Var1
data <- data[, colnames(data) != "Var1"]

data_corr <- cor(as.data.frame(data), method = "spearman")


short_names <- c(
    "Alveolar II"="ATII",
    "Basal"="Basal",
    "B cells"="B",
    "CD4+"="CD4",
    "CD4"="CD4",
    "CD8+"="CD8",
    "CD8"="CD8",
    "Ciliated"="Cilia",
    "Club"="Club",
    "DC"="DC",
    "Dendritic"="DC",
    "Endothelial"="EC",
    "Epithelial"="Epi",
    "Exhaust T"="Exh T",
    "Fibroblasts"="Fib",
    "Granulocyte"="Gran",
    "Mast"="Mast",
    "Mascrophages"="Mφ",
    "Macrophages"="Mφ",
    "Monocytes"="Mo",
    "Neuroendocrine"="NE",
    "NK"="NK",
    "Tregs"="Tregs"
)

short_names <- as.data.frame(short_names)

cell_types = short_names[sapply(colnames(data), function(x) { str_split(x, "\\.")[[1]][1] }), 1]
cell_color = c(
    wes_palette("Darjeeling1"), 
    wes_palette("Moonrise3"),
    wes_palette("Darjeeling2"),
    wes_palette("Moonrise1")
)[1:length(unique(cell_types))]


names(cell_color) = unique(cell_types)


colors = list(
    Stage = c(
        "I"="#65A9A3", 
        "II"="#4A933E", 
        "III"="#EC7A21", 
        "IV"="#D73F47", 
        "LUAD_Normal"="#FECC1B", 
        "LUSC_Normal"="#778793"
    ),
    Disease = c(
        "AD" = "#0084D1", 
        "NL(AD)"="#FECC1B", 
        "NL(AD"="#778793", 
        "Normal"="#73BDFF", 
        "NL(SC)"="#778793", 
        "SC"="#A0C807",
        "LC"="#91822B",
        "NL(LC)"="#EDCAB0",
        "PL"="#EBAAA4",
        "NL(PL)"="#B43018"
    ),
    Batch = c("1"="#E9B0B7", "2"="#90D0E0", "3"="#769982"),
    Cell = c(
        "ATII"="#FF0000",
        "Basal"="#F2AD00",
        "CD4"="#F98400",
        "CD8"="#5BBCD6",
        "Cilia"="#85D4E3",
        "Club"="#F4B5BD",
        "DC"="#9C964A",
        "EC"="#CDC08C",
        "Epi"="#FAD77B",
        "Fib"="#ECCBAE",
        "Mast"="#D69C4E",
        "Mo"="#ABDDDE",
        "NK"="#000000",
        "NE"="#F3DF6C",
        "Tregs"="#CEAB07",
        "Mφ"="#046C9A"
    )
)



ba <- HeatmapAnnotation(
    Disease = meta[sapply(colnames(data_corr), function(x) { str_split(x, " ")[[1]][2] }), "Disease"],
    Stage = meta[sapply(colnames(data_corr), function(x) { str_split(x, " ")[[1]][3] }), "Stage"],
    Cell=short_names[sapply(colnames(data_corr), function(x) { str_split(x, " ")[[1]][1] }), 1],
    col = colors,
    show_annotation_name = F,
    show_legend = F
)


ra <- rowAnnotation(
    Disease = meta[sapply(rownames(data_corr), function(x) { str_split(x, " ")[[1]][2] }), "Disease"],
    Stage = meta[sapply(rownames(data_corr), function(x) { str_split(x, " ")[[1]][3] }), "Stage"],
    Cell = short_names[sapply(rownames(data_corr), function(x) { str_split(x, " ")[[1]][1] }), 1],
    col = colors
)



Heatmap(
    data_corr,
    name = "Corr",
    show_column_names = F,
    bottom_annotation = ba,
    right_annotation = ra
)


library(stringr)
library(sctransform)

data = read.csv("Alveolar_II/scanpy/raw.csv", row.name = 1)
colnames(data) <- sapply(colnames(data), function(x) {
    x = str_replace(x, "^X", "")
    str_replace(x, "\\.", "-")
})

meta = read.csv("Alveolar_II/scanpy/meta.csv", row.name = 1)


obj <- CreateSeuratObject(
    data,
    meta.data = meta
)

obj <- PercentageFeatureSet(obj, pattern = "^MT-", col.name = "percent.mt")
obj <- SCTransform(obj, vars.to.regress = "Batch", verbose = T, variable.features.n=nrow(data))

saveRDS(obj@assays$SCT@scale.data, "Alveolar_II/scanpy/scale_data_sct.rds")

expr <- readRDS("03_each_cells/total/Alveolar_II/scanpy/scale_data_sct.rds")

bak_scale = atii@scale.data
atii@scale.data = as.matrix(expr)

importance <- read.csv("03_each_cells/total/Alveolar_II/scanpy/scale_data_sct_importance.csv", row.names = 1, stringsAsFactors = F)

pdf("03_each_cells/select_features/atii_wo_low_importance_markers.pdf", width = 10, height = 8)
make_stage_corr_heatmap(atii, features = importance$feature[importance$importance > 0 & importance$cumulative_importance > 0.99], row.names(atii@scale.data), markers = NULL)
dev.off()


make_stage_module_heatmap(
    atii, 
    intersect(importance$feature[importance$importance > 0], row.names(atii@scale.data)), 
    cluster=T,
    col_fun = colorRamp2(c(-3, 0, 3), c("blue", "white", "red"))
)



Heatmap(expr, cluster_rows = F, cluster_columns = F, show_row_names = F, show_column_names = F, col = colorRamp2(c(-2, 0, 6), c("blue", "white", "red")))


atii <- readRDS("03_each_cells/total/Alveolar_II/seurat.rds")


temp_gene = h@ht_list$Expr@row_names_param$labels[h@ht_list$Expr@row_order]
temp_gene_mark <- round(length(temp_gene) / 2)
temp_gene_mark <- temp_gene[(temp_gene_mark - 120): (temp_gene_mark + 30)]
temp_gene_mark <- c("ACSL4", "ALPL", "AQP4", "BMP2", "SFTPC", "C4BPA")


make_stage_module_heatmap(atii, temp_gene_mark, cluster = T)


atii@meta.data$ident <- paste(atii@meta.data[, "Disease"], atii@meta.data[, "Stage"], sep = "-")


expr <- readRDS("03_each_cells/total/Alveolar_II/scanpy/scale_data_sct.rds")

bak_scale = atii@scale.data
atii@scale.data = as.matrix(expr)


atii@meta.data$PatientID <- as.character(atii@meta.data$PatientID)

markers = NULL
for(i in unique(atii@meta.data$ident)) {
    print(i)
    ident.1 = rownames(atii@meta.data)[atii@meta.data$ident == i & atii@meta.data$PatientID != "PA11"]
    ident.2 = rownames(atii@meta.data)[atii@meta.data$ident != i & atii@meta.data$PatientID != "PA11"]
    
    print(length(ident.1))
    print(length(ident.2))
    
    if (length(ident.1) >=3 && length(ident.2) >= 3) {
        temp = FindMarkers(atii, ident.1 = ident.1, ident.2 = ident.2)
        temp$ident = i
        temp$gene = rownames(temp)
        markers = rbind(markers, temp)
    }
}



importance <- read.csv("03_each_cells/select_features/atii_importance.csv", row.names = 1)


make_stage_module_heatmap(atii, unique(markers[markers$p_val_adj < 0.05 & markers$avg_logFC > .5, "gene"]), cluster = T)


make_stage_module_heatmap(atii, unique(importance$feature[importance$importance > 0 & importance$cumulative_importance > 0.99]), cluster = T)


temp <- melt(as.matrix(atii@scale.data[unique(markers[markers$p_val_adj < 0.05 & markers$avg_logFC > 1, "gene"]), ]))


temp$ident = paste(atii@meta.data[as.character(temp$Var2), "Disease"], atii@meta.data[as.character(temp$Var2), "Stage"], sep = "-")

temp <- temp %>%
    group_by(ident, Var1) %>%
    mutate(value = mean(value)) %>%
    dplyr::select(Var1, ident, value) %>%
    unique() %>%
    as.data.frame()



ggplot(temp, aes(x=Var1, y=value, color = ident, group = 1)) +
    geom_line() +
    facet_grid(.~ident, scales = "free_x")


saveRDS(atii@raw.data, "Lung_cancer_10x/ATII/raw_data.rds")
saveRDS(atii@meta.data, "Lung_cancer_10x/ATII/meta_data.rds")
saveRDS(atii@scale.data, "Lung_cancer_10x/ATII/scale_data.rds")

saveRDS(markers, "Lung_cancer_10x/ATII/markers.rds")