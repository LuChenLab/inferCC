options(stringsAsFactors = F)
library(dplyr)
library(stringr)
library(reshape2)
library(rtracklayer)
library(doMC)
library(openxlsx)
library(wesanderson)
library(ComplexHeatmap)
library(ggrastr)
library(Seurat)
library(UBL)

meta <- read.csv("02_rds/meta_after_singleR.csv", row.names = 1, stringsAsFactors = F)

# load inferCNV results
files <- list.files("11_CNV/", pattern = "run.final.infercnv_obj", recursive = T, full.names = T)


### infercnb_obj with gene chrom
gene_pos <- read.table("09_bulk/RNA_seq/gene_pos.txt")
rownames(gene_pos) <- make.unique(as.character(gene_pos$V6))


### cluste cnv

cell_color = c(
    wes_palette("Darjeeling1"), 
    wes_palette("Moonrise3"),
    wes_palette("Darjeeling2"),
    wes_palette("Moonrise1")
)[1:length(unique(meta$cell_name1))]

names(cell_color) = unique(meta$cell_name1)

disease_colors = c(
    "LUAD" = "#0084D1", 
    "Normal (LUAD)"="#FECC1B", 
    "LUAD_Normal"="#FECC1B", 
    "Normal"="#73BDFF", 
    "Normal (LUSC)"="#778793", 
    "LUSC_Normal"="#778793", 
    "LUSC"="#A0C807", 
    "BSPN"="#6A0019", 
    "BSPN_Normal"="#C5000B"
)


chrom_color <- c(
    wes_palette("Darjeeling1"), 
    wes_palette("Darjeeling2"), 
    wes_palette("Moonrise1"), 
    wes_palette("Moonrise2"), 
    wes_palette("Moonrise3"),
    wes_palette("Royal1"),
    wes_palette("Royal2"),
    wes_palette("Chevalier1"),
    wes_palette("FantasticFox1"),
    wes_palette("GrandBudapest1"),
    wes_palette("GrandBudapest2")
)
chrom_color <- chrom_color[1:length(unique(gene_pos$V1))]
names(chrom_color) <- unique(gene_pos$V1)


stage_color =  c(
    "I"="#65A9A3", 
    "II"="#4A933E", 
    "III"="#EC7A21", 
    "IV"="#D73F47", 
    "LUAD_Normal"="#FECC1B", 
    "LUSC_Normal"="#778793",
    "Normal(LUAD)"="#FECC1B", 
    "Normal(LUSC)"="#778793"
    
)


#####
temp_gene <- gene_pos[gene_pos$V6 %in% infercnv_genes, ]
temp_gene$V1 <- paste0("chr", temp_gene$V1)


write.table(
    temp_gene[, c("V1", "V2", "V3", "V6", "V5", "V4")],
    "11_CNV/infercnv_required_genes_hg38.bed", 
    sep = "\t", 
    quote = F,
    row.names = F,
    col.names = F
)

files = list.files("09_bulk/DNA/computedMatrix/", pattern = ".csv", full.names = T)



registerDoMC(10)
cnv = foreach(i = files, .combine="rbind") %dopar% {
    # print(i)
    temp = read.csv(i, header = F)
    
    sample_info <- str_split(basename(i), "[._]")[[1]]
    
    temp$SampleID <- sample_info[1]
    temp$SampleType <- str_to_title(sample_info[2])
    
    colnames(temp)[1:5] <- c("chrom", "start", "end", "strand", "gene_name")
    
    temp
}


saveRDS(cnv, "09_bulk/DNA/computedMatrix/wgs.rds")

wgs <- readRDS("09_bulk/DNA/computedMatrix/wgs.rds")


bulk_meta <- read.xlsx("09_bulk/RNATAC50commonSample-1119.xlsx")
bulk_meta <- bulk_meta[!is.na(bulk_meta$SampleID), ]
rownames(bulk_meta) <- bulk_meta$SampleID
bulk_meta$Stage <- str_replace(bulk_meta$Stage, "[^IV]", "")
bulk_meta$Stage[is.na(bulk_meta$Stage)] <- "BSPN"



make_point_plot_of_wgs <- function(gene) {
    temp_wgs = wgs[wgs$gene_name == gene, ]
    rownames(temp_wgs) <- paste(temp_wgs$SampleID, temp_wgs$SampleType, sep = "-")
    temp_data = melt(as.matrix(temp_wgs[, str_detect(colnames(temp_wgs), "^V")]))
    temp_data$Var2 <- as.numeric(as.character(str_replace_all(temp_data$Var2, "^V", "")))
    
    temp_data$SampleType = sapply(temp_data$Var1, function(x) {str_split(x, "-")[[1]][2]})
    temp_data$SampleID = sapply(temp_data$Var1, function(x) { str_replace(str_split(x, "-")[[1]][1], "\\d+", "") })
    
    temp_data$SampleType = factor(temp_data$SampleType, levels = c("Tumor", "Normal"))
    
    p <- ggplot(
        temp_data[temp_data$SampleID %in% c("LUAD", "LUSC"), ], 
        aes(x=Var2, y = value, color = SampleType)
    ) +
        geom_point(size = 0.4) +
        facet_grid(.~SampleID+SampleType) +
        geom_smooth(method = loess) +
        labs(x = "", y = "", color = "") +
        theme_bw() +
        theme(
            axis.text.x = element_blank(),
            legend.position = "bottom",
            legend.text = element_text(size = 15),
            strip.text = element_text(size = 15),
            axis.text.y = element_text(size = 10)
        )
    
    return (p)
}


target_genes = c(
    "CYP2B7P", "CKS1B", "IGHG1",
    "MUC1", "PIGR", "SFTPA1", 
    "SFTPA2", "PIGR", "VIM", 
    "SPINK1", "PTMA", "IGFBP7"
)

registerDoMC(10)
res = foreach(i = target_genes, .errorhandling = "pass") %dopar% {
    ggsave(
        filename = paste0("09_bulk/WGS/", i, ".pdf"),
        plot = make_point_plot_of_wgs(i),
        width = 6,
        height = 4
    )
}


ggsave(
    filename = paste0("09_bulk/WGS/UBB.pdf"),
    plot = make_point_plot_of_wgs("UBB"),
    width = 6,
    height = 4
)

ggsave(
    filename = paste0("09_bulk/WGS/CCND1.pdf"),
    plot = make_point_plot_of_wgs("CCND1"),
    width = 6,
    height = 4
)


make_point_plot_of_wgs("H2AFZ")

# 
# ra <- rowAnnotation(
#     Sample=str_replace(temp_wgs$SampleID, "\\d+", ""),
#     Type = temp_wgs$SampleType,
#     Stage = bulk_meta[as.character(temp_wgs$SampleID), "Stage"],
#     col = list(
#         Sample = c(
#             "LUAD" = "#0084D1", 
#             "Normal"="#73BDFF", 
#             "LUSC"="#A0C807", 
#             "BSPN"="#6A0019", 
#             "BSPN_Normal"="#C5000B"
#         ),
#         Type = c(
#             "Tumor"="#E9B0B7",
#             "Normal"="#90D0E0"
#         ),
#         Stage = c(
#             "I"="#65A9A3", 
#             "II"="#4A933E", 
#             "III"="#EC7A21", 
#             "IV"="#D73F47", 
#             "LUAD_Normal"="#FECC1B", 
#             "LUSC_Normal"="#778793",
#             "BSPN"="#6A0019"
#             
#         )
#     ),
#     annotation_legend_param = list(
#         Sample = list(
#             labels_gp = gpar(fontsize = 12),
#             title_gp = gpar(fontsize = 15)
#         ),
#         Type = list(
#             labels_gp = gpar(fontsize = 12),
#             title_gp = gpar(fontsize = 15)
#         )
#     )
# )
# 
# draw(Heatmap(
#     as.matrix(scale(temp_wgs[, str_detect(colnames(temp_wgs), "^V")])),
#     name = "",
#     cluster_rows = T,
#     cluster_columns = F,
#     show_row_names = F,
#     show_column_names = F,
#     right_annotation = ra,
#     column_title = "CCDC33",
#     column_title_gp = gpar(fontsize = 30)
#     # heatmap_legend_param = list(
#     #     labels_gp = gpar(fontsize = 18),
#     #     title_gp = gpar(fontsize = 20)
#     # )
# ))
# 
# 
# h <- readRDS("03_each_cells/select_features/ATII_importance_gene_heatmap.rds")
# temp_gene = h@ht_list$Expr@row_names_param$labels[h@ht_list$Expr@row_order]
# temp_gene_mark <- round(length(temp_gene) / 2)
# temp_gene_mark <- temp_gene[(temp_gene_mark - 120): (temp_gene_mark + 30)]
# 
# 
# for (i in intersect(wgs$gene_name, new_markers$gene[new_markers$ident == "LUAD-I"])) {
#     print(i)
#     temp_wgs = wgs[wgs$gene_name == i, ]
#     
#     ra <- rowAnnotation(
#         Sample=str_replace(temp_wgs$SampleID, "\\d+", ""),
#         Type = temp_wgs$SampleType,
#         Stage = bulk_meta[as.character(temp_wgs$SampleID), "Stage"],
#         col = list(
#             Sample = c(
#                 "LUAD" = "#0084D1", 
#                 "Normal"="#73BDFF", 
#                 "LUSC"="#A0C807", 
#                 "BSPN"="#6A0019", 
#                 "BSPN_Normal"="#C5000B"
#             ),
#             Type = c(
#                 "Tumor"="#E9B0B7",
#                 "Normal"="#90D0E0"
#             ),
#             Stage = c(
#                 "I"="#65A9A3", 
#                 "II"="#4A933E", 
#                 "III"="#EC7A21", 
#                 "IV"="#D73F47", 
#                 "LUAD_Normal"="#FECC1B", 
#                 "LUSC_Normal"="#778793",
#                 "BSPN"="#6A0019"
#                 
#             )
#         ),
#         annotation_legend_param = list(
#             Sample = list(
#                 labels_gp = gpar(fontsize = 12),
#                 title_gp = gpar(fontsize = 15)
#             ),
#             Type = list(
#                 labels_gp = gpar(fontsize = 12),
#                 title_gp = gpar(fontsize = 15)
#             )
#         )
#     )
#     
#     pdf(paste0("09_bulk/DNA/WGS/", i, ".pdf"), width = 10, height = 6)
#     draw(Heatmap(
#         as.matrix(scale(temp_wgs[, str_detect(colnames(temp_wgs), "^V")])),
#         name = "",
#         cluster_rows = T,
#         cluster_columns = F,
#         show_row_names = F,
#         show_column_names = F,
#         right_annotation = ra,
#         column_title = i,
#         column_title_gp = gpar(fontsize = 30)
#         # heatmap_legend_param = list(
#         #     labels_gp = gpar(fontsize = 18),
#         #     title_gp = gpar(fontsize = 20)
#         # )
#     ))
#     dev.off()
# }
# 
# 
# 
# ### Make heatmap
# gene_pos <- read.table("09_bulk/RNA_seq/gene_pos.txt")
# rownames(gene_pos) <- make.unique(as.character(gene_pos$V6))
# 
# 
# cell_color = c(
#     wes_palette("Darjeeling1"), 
#     wes_palette("Moonrise3"),
#     wes_palette("Darjeeling2"),
#     wes_palette("Moonrise1")
# )[1:length(unique(meta$cell_name1))]
# 
# names(cell_color) = unique(meta$cell_name1)
# 
# disease_colors = c(
#     "LUAD" = "#0084D1", 
#     "Normal (LUAD)"="#FECC1B", 
#     "LUAD_Normal"="#FECC1B", 
#     "Normal"="#73BDFF", 
#     "Normal (LUSC)"="#778793", 
#     "LUSC_Normal"="#778793", 
#     "LUSC"="#A0C807", 
#     "BSPN"="#6A0019", 
#     "BSPN_Normal"="#C5000B"
# )
# 
# 
# chrom_color <- c(
#     wes_palette("Darjeeling1"), 
#     wes_palette("Darjeeling2"), 
#     wes_palette("Moonrise1"), 
#     wes_palette("Moonrise2"), 
#     wes_palette("Moonrise3"),
#     wes_palette("Royal1"),
#     wes_palette("Royal2"),
#     wes_palette("Chevalier1"),
#     wes_palette("FantasticFox1"),
#     wes_palette("GrandBudapest1"),
#     wes_palette("GrandBudapest2")
# )
# chrom_color <- chrom_color[1:length(unique(gene_pos$V1))]
# names(chrom_color) <- unique(gene_pos$V1)
# 
# 
# files <- list.files("11_CNV/", pattern = "run.final.infercnv_obj", recursive = T, full.names = T)
# 
# infercnv_genes <- c()
# res = list()
# for (i in files) {
#     print(i)
#     
#     if (str_detect(i, "(Neuroendocrine|Alveolar_II|Ciliated)")) {
#         temp = read_infercnv_obj(i, FALSE)
#         res[[basename(dirname(i))]] <- temp
#         
#         infercnv_genes <- unique(c(infercnv_genes, rownames(temp)))
#     }
# }
# 
# temp_data = NULL
# for (i in names(res)) {
#     j <- res[[i]]
#     j <- j[, !colnames(j) %in% c("chr", "start", "stop", "cell")]
#     print(dim(j))
#     if (is.null(temp_data)) {
#         temp_data = j
#     } else {
#         genes = intersect(rownames(temp_data), rownames(j))
#         temp_data = cbind(temp_data[genes, ], j[genes, ])
#     }
# }
# 
# 
# 
# meta <- read.csv("02_rds/meta.csv", row.names = 1, stringsAsFactors = F)
# temp_meta <- meta[as.character(colnames(temp_data)), ]
# 
# 
# 
# 
# registerDoMC(10)
# heatmaps = foreach (cell = unique(temp_meta$cell_name1)) %dopar% {
#     
#     if (cell == "WGS") {
#         return()
#     }
#     
#     print(cell)
#     
#     temp_cells <- rownames(temp_meta)[temp_meta$cell_name1 == cell]
#     temp_cells <- c(temp_cells, rownames(temp_meta)[temp_meta$cell_name1 == "WGS"])
#     
#     la <- rowAnnotation(
#         # cell = temp_meta[temp_cells, "cell_name1"],
#         disease = temp_meta[temp_cells, "Disease"],
#         col = list(
#             # cell = cell_color,
#             disease = disease_colors
#         ),
#         show_annotation_name = F
#     )
#     
#     temp_pos = data.frame(
#         gene = rownames(temp_data),
#         chrom = gene_pos[rownames(temp_data), "V1"]
#     )
#     
#     temp_pos <- temp_pos[order(temp_pos$chrom), ]
#     
#     temp_anno_genes = which(rownames(temp_data) %in% head(unique(wgs$gene_name), 20))
#     ba <- HeatmapAnnotation(
#         # chrom = temp_pos$chrom,
#         # col = list(
#         #     chrom = chrom_color
#         # )
#         foo = anno_mark(
#             at=temp_anno_genes,
#             labels=rownames(temp_data)[temp_anno_genes],
#             which = "column",
#             side = "bottom"
#         ),
#         show_annotation_name = F
#     )
#     
#     p <- Heatmap(
#         t(as.matrix(temp_data[temp_pos$gene, temp_cells])),
#         name = "",
#         column_title = cell,
#         column_title_gp = gpar(fontsize = 25),
#         show_row_names = F,
#         left_annotation = la,
#         row_title_rot = 90,
#         border = T,
#         show_column_names = F,
#         bottom_annotation = ba,
#         cluster_columns = T,
#         column_split = as.character(temp_pos$chrom),
#         row_title_gp = gpar(fontsize = 15)
#     )
#     
#     
#     pdf(paste0("11_CNV/compare_with_bulk/", cell, ".pdf"), width = 15, height = 8)
#     p <- draw(p)
#     dev.off()
#     
#     return(p)
# }
# 
# 