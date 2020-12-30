library(dtangle)
library(ggplot2)
library(openxlsx)
library(here)
library(reshape2)
library(ggplot2)
library(stringr)
library(ComplexHeatmap)
library(circlize)


options(stringsAsFactors = F)


short_names <- c(
    "APII"="ATII",
    "Alveolar II"="ATII",
    "Basal"="Basal",
    "B.cells"="B",
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


## prepare single cell reference
# scExpr <- readRDS(here("02_rds/all_cell_expr.rds"))
# meta <- read.csv(here("02_rds/meta_after_singleR.csv"), row.names = 1)
# meta <- meta[meta$Batch != 3 & meta$cell_name != "", ]
# scExpr <- scExpr[, rownames(meta)]
# 
# ### Get mean by cell types
# scExpr <- melt(as.matrix(scExpr))
# 
# 
# write.csv(scExpr, "02_rds/raw_data.csv")


references <- read.table("raw_data_mean_by_cell.tsv", sep = "\t", header = T)
rownames(references) <- make.unique(references$gene)
colnames(references) <- short_names[colnames(references), 1]
references <- references[, !is.na(colnames(references)) & colnames(references) != "Mφ"]


### dtangle on bulk
expr <- read.xlsx("09_bulk/RNA_seq/dtangle/RNA_seq/RSEM.xlsx", rowNames = T)

gene_pos <- read.table("09_bulk/RNA_seq/dtangle/RNA_seq/gene_pos.txt", header = F, row.names = 5)
gene_pos$V6 <- make.unique(gene_pos$V6)
rownames(expr) <- gene_pos[rownames(expr), "V6"]


common_genes <- intersect(rownames(expr), rownames(references))


# read markers
markers = read.xlsx("20190701_gene_markers.xlsx", sheet = 3)
markers <- na.omit(markers[, 1:2])
markers <- unique(markers)

markers <- markers[markers$Cells %in% c(
    "Alveolar II", "B cells", "Basal", "CD4+",                    
    "CD8+", "Ciliated", "Club",
    "Dendritic", "Endothelial", "Epithelial", 
    "Fibroblasts", "Granulocyte",
    "Mast", "Monocytes", "Neuroendocrine", "NK",
    "T cells", "Treg"   
), ]


res = dtangle(
    as.matrix(t(log2(expr[common_genes, ] + 1))), 
    references = as.matrix(t(log2(references[common_genes, ] + 1))),
    markers = which(common_genes %in% markers$Markers)
)


head(res$estimates)

saveRDS(res, "09_bulk/RNA_seq/dtangle/estimate.rds")



### bulk on stage
bulk_meta <- read.xlsx("09_bulk/RNATAC50commonSample-1119.xlsx")
bulk_meta <- bulk_meta[!is.na(bulk_meta$SampleID), ]

rownames(bulk_meta) <- bulk_meta$SampleID
bulk_meta$Stage <- str_replace(bulk_meta$Stage, "[^IV]", "")



### 
res <- readRDS("09_bulk/RNA_seq/dtangle/es")

immu = c(
    "B",
    "CD4",
    "CD8",
    "DC",
    "EC",
    "Gran",
    "Mast",
    "Mo",
    "NK",
    "Tregs",
    "Mφ"
)

res <- readRDS("09_bulk/RNA_seq/dtangle/estimate.rds")

temp <- melt(as.matrix(res$estimates))
temp$Disease <- str_replace_all(temp$Var1, "\\d+", "")
temp$Type <- ifelse(temp$Var2 %in% immu, "Immune", "Non-immune") 


p <- ggplot(temp, aes(x=Var2, y=value, color = Disease)) +
    geom_boxplot() +
    facet_grid(.~Type, scale = "free", space = "free") +
    theme_bw() +
    theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 20),
        axis.text.y = element_text(size = 15),
        strip.text = element_text(size = 25),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.position = "top"
    ) +
    labs(x = "", y = "")

ggsave(
    filename = "09_bulk/RNA_seq/dtangle/dtangle_boxplot.pdf",
    plot = p,
    width = 12,
    height = 8
)



###

make_heatmap_bulk <- function(
    data,
    cluster_method = "ward.D2",
    border = TRUE,
    scale_dot_size = 0.9,
    stage=NULL
) {
    library(ComplexHeatmap)
    library(stringr)
    library(circlize)
    
    colors = list(
        Stage = c("I"="#65A9A3", "II"="#4A933E", "III"="#EC7A21", "IV"="#D73F47", "LUAD_Normal"="#FECC1B", "LUSC_Normal"="#778793", "BSPN"="#73BDFF"),
        Disease = c("LUAD" = "#0084D1", "NL(AD)"="#FECC1B", "NL(AD"="#778793", "BSPN"="#73BDFF", "NL(SC)"="#778793", "LUSC"="#A0C807"),
        Batch = c("1"="#E9B0B7", "2"="#90D0E0")
    )
    
    if (is.null(stage)) {
        ra <- rowAnnotation(
            Disease = str_replace_all(rownames(data), "\\d+", ""),
            Stage = str_replace_all(bulk_meta[rownames(data), "Stage"], "[^IV]", ""),
            col = colors,
            show_annotation_name = F,
            show_legend = T
        )
        
        ba <- HeatmapAnnotation(
            Disease = str_replace_all(colnames(data), "\\d+", ""),
            Stage = str_replace_all(bulk_meta[colnames(data), "Stage"], "[^IV]", ""),
            col = colors,
            show_annotation_name = T,
            show_legend = F
        )
    } else {
        ra <- rowAnnotation(
            Stage = sapply(rownames(data), function(x) { stage[x] }),
            col = colors,
            show_annotation_name = F,
            show_legend = T
        )
        
        ba <- HeatmapAnnotation(
            Stage = sapply(colnames(data), function(x) { stage[x] }),
            col = colors,
            show_annotation_name = T,
            show_legend = F
        )
    }
    
    
    
    col_fun = colorRamp2(c(-1, 0, 1), c("green", "white", "red"))

    Heatmap(
        data,
        name = "Corr",
        col = col_fun,
        border = border,
        right_annotation = ra,
        show_row_names = T,
        show_column_names = F,
        bottom_annotation = ba,
        clustering_method_rows = cluster_method,
        clustering_method_columns = cluster_method,
        # rect_gp = gpar(col = "white", lwd = 2)
        rect_gp = gpar(type = "none"),
        cell_fun = function(j, i, x, y, width, height, fill) {
            s = min(unit.c(convertWidth(width, "cm"), convertHeight(height, "cm")))
            grid.circle(
                x = x,
                y = y, 
                r = s * (data[i, j] + 0.1) * scale_dot_size,
                gp = gpar(col = "white", fill = fill)
            )
        }
    )
}

est <- res$estimates

pdf("09_bulk/RNA_seq/dtangle/dtangle_heatmap_all_sample.pdf", width = 10, height = 8)
make_heatmap_bulk(cor(t(est)), scale_dot_size = 0.4)
dev.off()


###
immu <- c(
    "B",
    "CD4",
    "CD8",
    "DC",
    "EC",
    "Mast",
    "Mo",
    "NK",
    "Tregs"
)

pdf("09_bulk/RNA_seq/dtangle/dtangle_heatmap_immune.pdf", width = 10, height = 8)
make_heatmap_bulk(cor(t(est[, colnames(est) %in% immu])), scale_dot_size = 0.4)
dev.off()


pdf("09_bulk/RNA_seq/dtangle/dtangle_heatmap_non_immune.pdf", width = 10, height = 8)
make_heatmap_bulk(cor(t(est[, !colnames(est) %in% immu])), scale_dot_size = 0.4)
dev.off()



### Stage
for (i in unique(bulk_meta$Stage)) {
    
    if (!is.na(i)) {
        print(i)
        
        pdf(paste0("09_bulk/RNA_seq/dtangle/dtangle_heatmap_", i, ".pdf"), width = 10, height = 8)
        draw(
            make_heatmap_bulk(
                cor(t(est[intersect(rownames(est), bulk_meta$SampleID[bulk_meta$Stage == i]), ])), 
                scale_dot_size = 0.4
            )
        )
        dev.off()
        
        
        pdf(paste0("09_bulk/RNA_seq/dtangle/dtangle_heatmap_", i, "_non_immune.pdf"), width = 10, height = 8)
        draw(
            make_heatmap_bulk(
                cor(t(est[intersect(rownames(est), bulk_meta$SampleID[bulk_meta$Stage == i]), !colnames(est) %in% immu])), 
                scale_dot_size = 0.4
            )
        )
        dev.off()
        
        
        pdf(paste0("09_bulk/RNA_seq/dtangle/dtangle_heatmap_", i, "_immune.pdf"), width = 10, height = 8)
        draw(
            make_heatmap_bulk(
                cor(t(est[intersect(rownames(est), bulk_meta$SampleID[bulk_meta$Stage == i]), colnames(est) %in% immu])), 
                scale_dot_size = 0.4
            )
        )
        dev.off()
    }
}


### Disease
bulk_meta$Disease = str_replace_all(rownames(bulk_meta), "\\d+", "")
stages <- bulk_meta$Stage
stages[is.na(stages)] <- "BSPN"
names(stages) <- rownames(bulk_meta)

for (i in unique(bulk_meta$Disease)) {
    
    if (!is.na(i)) {
        print(i)
        
        pdf(paste0("09_bulk/RNA_seq/dtangle/dtangle_heatmap_", i, ".pdf"), width = 10, height = 8)
        draw(
            make_heatmap_bulk(
                cor(t(est[intersect(rownames(est), bulk_meta$SampleID[bulk_meta$Disease == i]), ])), 
                scale_dot_size = 0.4,
                stage = stages
            )
        )
        dev.off()
        
        
        pdf(paste0("09_bulk/RNA_seq/dtangle/dtangle_heatmap_", i, "_non_immune.pdf"), width = 10, height = 8)
        draw(
            make_heatmap_bulk(
                cor(t(est[intersect(rownames(est), bulk_meta$SampleID[bulk_meta$Disease == i]), !colnames(est) %in% immu])), 
                scale_dot_size = 0.4,
                stage = stages
            )
        )
        dev.off()
        
        
        pdf(paste0("09_bulk/RNA_seq/dtangle/dtangle_heatmap_", i, "_immune.pdf"), width = 10, height = 8)
        draw(
            make_heatmap_bulk(
                cor(t(est[intersect(rownames(est), bulk_meta$SampleID[bulk_meta$Disease == i]), colnames(est) %in% immu])), 
                scale_dot_size = 0.4,
                stage = stages
            )
        )
        dev.off()
    }
}


##### 

### Dtangle on TCGA

# read_tcga <- function(path) {
#     data = read.csv(path, stringsAsFactors = F)
#     rownames(data) <- make.unique(sapply(data$X, function(x) { str_split(x, "\\|")[[1]][1] }))
#     data[, colnames(data) != "X"]
# }
# 
# luad <- read_tcga("TCGA/LUAD_raw_expresion.csv")
# lusc <- read_tcga("TCGA/LUSC_raw_expresion.csv")
# 
# tcga <- cbind(luad, lusc)
# 
# luad <- readRDS("TCGA/LUAD_clinical.rds")
# lusc <- readRDS("TCGA/LUSC_clinical.rds")
# 
# clinical <- rbind(luad, lusc)
# rownames(clinical) <- clinical$submitter_id
# 
# colnames(tcga) <- sapply(colnames(tcga), function(x) {
#     submit <- str_split(x, "\\.")[[1]]
#     submit <- paste(submit[1], submit[2], submit[3], sep = "-")
#     paste(x, clinical[submit, "disease"], str_to_upper(str_replace_all(clinical[submit, "tumor_stage"], "[^iv]", "")), sep = "-")
# })
# 
# write.csv(tcga, "TCGA/Lung.csv")
tcga <- read.csv("TCGA/Lung.csv", row.names = 1, stringsAsFactors = F)


references <- read.table("raw_data_mean_by_cell.tsv", sep = "\t", header = T)
rownames(references) <- make.unique(references$gene)
references <- references[, colnames(references) != "gene"]
colnames(references) <- short_names[colnames(references), 1]
references <- references[, colnames(references) %in% short_names[, 1] & colnames(references) != "Mφ"]


common_genes <- intersect(rownames(tcga), rownames(references))


# read markers
markers = read.xlsx("20190701_gene_markers.xlsx", sheet = 3)
markers <- na.omit(markers[, 1:2])
markers <- unique(markers)

markers <- markers[markers$Cells %in% c(
    "Alveolar II", "B cells", "Basal", "CD4+",                    
    "CD8+", "Ciliated", "Club",
    "Dendritic", "Endothelial", "Epithelial", 
    "Fibroblasts", "Granulocyte",
    "Mast", "Monocytes", "Neuroendocrine", "NK",
    "T cells", "Treg"   
), ]


res = dtangle(
    as.matrix(t(log2(tcga[common_genes, ] + 1))), 
    references = as.matrix(t(log2(references[common_genes, ] + 1))),
    markers = which(common_genes %in% markers$Markers)
)


head(res$estimates)

saveRDS(res, "09_bulk/RNA_seq/dtangle/tcga_estimate.rds")


immu = c(
    "B",
    "CD4",
    "CD8",
    "DC",
    "EC",
    "Gran",
    "Mast",
    "Mo",
    "NK",
    "Tregs",
    "Mφ"
)

res1 <- readRDS("09_bulk/RNA_seq/dtangle/tcga_estimate.rds")

temp <- melt(as.matrix(res$estimates))

temp$Disease <- sapply(as.character(temp$Var1), function(x) {
    str_split(x, "\\.")[[1]][8]
})

temp$Type <- ifelse(temp$Var2 %in% immu, "Immune", "Non-immune") 


p <- ggplot(temp, aes(x=Var2, y=value, color = Disease)) +
    geom_boxplot() +
    facet_grid(.~Type, scale = "free", space = "free") +
    theme_bw() +
    theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 20),
        axis.text.y = element_text(size = 15),
        strip.text = element_text(size = 25),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.position = "top"
    ) +
    labs(x = "", y = "") +
    scale_y_log10()

p

ggsave(
    filename = "09_bulk/RNA_seq/dtangle/dtangle_tcga_boxplot.pdf",
    plot = p,
    width = 12,
    height = 8
)




make_heatmap_tcga <- function(
    data,
    cluster_method = "ward.D2",
    border = TRUE,
    scale_dot_size = 0.9,
    stage=NULL
) {
    library(ComplexHeatmap)
    library(stringr)
    library(circlize)
    
    colors = list(
        Stage = c("I"="#65A9A3", "II"="#4A933E", "III"="#EC7A21", "IV"="#D73F47", "LUAD_Normal"="#FECC1B", "LUSC_Normal"="#778793", "BSPN"="#73BDFF", ""),
        Disease = c("LUAD" = "#0084D1", "NL(AD)"="#FECC1B", "NL(AD"="#778793", "BSPN"="#73BDFF", "NL(SC)"="#778793", "LUSC"="#A0C807"),
        Batch = c("1"="#E9B0B7", "2"="#90D0E0")
    )
    
    row_sel = sapply(rownames(data), function(x) { str_split(x, "\\.")[[1]][9] })
    row_sel = row_sel != ""
    
    col_sel = sapply(colnames(data), function(x) { str_split(x, "\\.")[[1]][9] })
    col_sel = col_sel != ""
  
    print(dim(data))
    data = data[row_sel, col_sel]
    print(dim(data))
    
    ra <- rowAnnotation(
        Disease = sapply(rownames(data), function(x) { str_split(x, "\\.")[[1]][8] }),
        Stage = sapply(rownames(data), function(x) { str_split(x, "\\.")[[1]][9] }),
        col = colors,
        show_annotation_name = F,
        show_legend = T
    )
    
    ba <- HeatmapAnnotation(
        Disease = sapply(colnames(data), function(x) { str_split(x, "\\.")[[1]][8] }),
        Stage = sapply(colnames(data), function(x) { str_split(x, "\\.")[[1]][9] }),
        col = colors,
        show_annotation_name = T,
        show_legend = F
    )
    
    col_fun = colorRamp2(c(-1, 0, 1), c("green", "white", "red"))

    Heatmap(
        data,
        name = "Corr",
        col = col_fun,
        border = border,
        right_annotation = ra,
        show_row_names = F,
        show_column_names = F,
        bottom_annotation = ba,
        cluster_rows = T,
        cluster_columns = T
        # clustering_method_rows = cluster_method,
        # clustering_method_columns = cluster_method,
        # rect_gp = gpar(type = "none"),
        # cell_fun = function(j, i, x, y, width, height, fill) {
        #     s = min(unit.c(convertWidth(width, "cm"), convertHeight(height, "cm")))
        #     grid.circle(
        #         x = x,
        #         y = y,
        #         r = s * (data[i, j] + 0.1) * scale_dot_size,
        #         gp = gpar(col = "white", fill = fill)
        #     )
        # }
    )
}

est <- res$estimates

Heatmap(
    cor(t(est)),
    show_row_names = F,
    show_column_names = F
)

pdf("09_bulk/RNA_seq/dtangle/dtangle_tcga_heatmap_all_sample.pdf", width = 10, height = 8)
make_heatmap_tcga(cor(t(est)), scale_dot_size = 0.4)
dev.off()


###
immu <- c(
    "B",
    "CD4",
    "CD8",
    "DC",
    "EC",
    "Mast",
    "Mo",
    "NK",
    "Tregs"
)

pdf("09_bulk/RNA_seq/dtangle/dtangle_heatmap_immune.pdf", width = 10, height = 8)
make_heatmap_bulk(cor(t(est[, colnames(est) %in% immu])), scale_dot_size = 0.4)
dev.off()


pdf("09_bulk/RNA_seq/dtangle/dtangle_heatmap_non_immune.pdf", width = 10, height = 8)
make_heatmap_bulk(cor(t(est[, !colnames(est) %in% immu])), scale_dot_size = 0.4)
dev.off()



### bulk on stage
bulk_meta <- read.xlsx("09_bulk/RNA_seq/dtangle/RNATAC50commonSample-1119.xlsx")
bulk_meta <- bulk_meta[!is.na(bulk_meta$SampleID), ]
bulk_meta <- bulk_meta[2:nrow(bulk_meta), ]
rownames(bulk_meta) <- bulk_meta$SampleID
bulk_meta$Stage <- str_replace(bulk_meta$Stage, "[^IV]", "")


### Stage
for (i in unique(bulk_meta$Stage)) {
    
    if (!is.na(i)) {
        print(i)
        
        pdf(paste0("09_bulk/RNA_seq/dtangle/dtangle_heatmap_", i, ".pdf"), width = 10, height = 8)
        draw(
            make_heatmap_bulk(
                cor(t(est[intersect(rownames(est), bulk_meta$SampleID[bulk_meta$Stage == i]), ])), 
                scale_dot_size = 0.4
            )
        )
        dev.off()
        
        
        pdf(paste0("09_bulk/RNA_seq/dtangle/dtangle_heatmap_", i, "_non_immune.pdf"), width = 10, height = 8)
        draw(
            make_heatmap_bulk(
                cor(t(est[intersect(rownames(est), bulk_meta$SampleID[bulk_meta$Stage == i]), !colnames(est) %in% immu])), 
                scale_dot_size = 0.4
            )
        )
        dev.off()
        
        
        pdf(paste0("09_bulk/RNA_seq/dtangle/dtangle_heatmap_", i, "_immune.pdf"), width = 10, height = 8)
        draw(
            make_heatmap_bulk(
                cor(t(est[intersect(rownames(est), bulk_meta$SampleID[bulk_meta$Stage == i]), colnames(est) %in% immu])), 
                scale_dot_size = 0.4
            )
        )
        dev.off()
    }
}


### Disease
bulk_meta$Disease = str_replace_all(rownames(bulk_meta), "\\d+", "")
stages <- bulk_meta$Stage
stages[is.na(stages)] <- "BSPN"
names(stages) <- rownames(bulk_meta)

for (i in unique(bulk_meta$Disease)) {
    
    if (!is.na(i)) {
        print(i)
        
        pdf(paste0("09_bulk/RNA_seq/dtangle/dtangle_heatmap_", i, ".pdf"), width = 10, height = 8)
        draw(
            make_heatmap_bulk(
                cor(t(est[intersect(rownames(est), bulk_meta$SampleID[bulk_meta$Disease == i]), ])), 
                scale_dot_size = 0.4,
                stage = stages
            )
        )
        dev.off()
        
        
        pdf(paste0("09_bulk/RNA_seq/dtangle/dtangle_heatmap_", i, "_non_immune.pdf"), width = 10, height = 8)
        draw(
            make_heatmap_bulk(
                cor(t(est[intersect(rownames(est), bulk_meta$SampleID[bulk_meta$Disease == i]), !colnames(est) %in% immu])), 
                scale_dot_size = 0.4,
                stage = stages
            )
        )
        dev.off()
        
        
        pdf(paste0("09_bulk/RNA_seq/dtangle/dtangle_heatmap_", i, "_immune.pdf"), width = 10, height = 8)
        draw(
            make_heatmap_bulk(
                cor(t(est[intersect(rownames(est), bulk_meta$SampleID[bulk_meta$Disease == i]), colnames(est) %in% immu])), 
                scale_dot_size = 0.4,
                stage = stages
            )
        )
        dev.off()
    }
}

