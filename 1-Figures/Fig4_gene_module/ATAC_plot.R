#!/usr/bin/env python3
# Created at 2019.10.11
options(stringAsfactors = F)
library(stringr)
library(ggplot2)
library(doMC)
library(reshape2)
library(dplyr)
library(openxlsx)
library(ComplexHeatmap)
library(circlize)
library(ggrastr)


files = list.files("09_bulk/ATAC/computedMatrix/", pattern = ".csv", full.names = T)

registerDoMC(10)
cnv = foreach(i = files, .combine="rbind") %dopar% {
    # print(i)
    temp = read.csv(i)
    temp$ident = str_split(basename(i), "_")[[1]][2]
    
    temp$chrom <- sapply(temp$X, function(x) {
        str_split(x, ":")[[1]][1]
    })
    
    temp$start <- sapply(temp$X, function(x) {
        str_split(str_split(x, ":")[[1]][2], "-")[[1]][1]
    })
    
    temp$end <- sapply(temp$X, function(x) {
        str_split(str_split(x, ":")[[1]][2], "-")[[1]][2]
    })
    
    temp$strand <- sapply(temp$X, function(x) {
        str_split(str_split(x, ":")[[1]][3], "\t")[[1]][1]
    })
    
    temp$gene_id <- sapply(temp$X, function(x) {
        str_split(x, "\t")[[1]][2]
    })
    
    temp
}



gene_pos <- read.table("09_bulk/RNA_seq/gene_pos.txt")
rownames(gene_pos) <- gene_pos$V5

gene_trans <- read.table("09_bulk/gene_iso.txt", header = T)
rownames(gene_trans) <- gene_trans$transcript_id

cnv$trans_id <- sapply(cnv$gene_id, function(x){str_split(x, "\\.")[[1]][1]})
cnv$gene_id <- gene_trans[as.character(cnv$trans_id), "gene_id"]
cnv$gene_name <- gene_pos[as.character(cnv$gene_id), "V6"]


saveRDS(cnv, "09_bulk/ATAC/computedMatrix/all.rds")


#### 
cnv <- readRDS("09_bulk/ATAC/computedMatrix/all.rds")

# registerDoMC(10)
# data = foreach(i = sort(unique(cnv$ident)), .combine="rbind") %dopar% {
#     temp = cnv[cnv$ident == i, str_detect(colnames(cnv), "X\\d+")]
#     temp[is.na(temp)] = 0
#     
#     temp <- as.data.frame(apply(temp, 2, sum))
#     colnames(temp) = "y"
#     temp$x = 1:nrow(temp)
#     temp$ident = i
# 
#     temp
# }


temp <- unique(cnv[, c("X", "gene_id", "trans_id", "gene_name")])

### select genes by clustering
#### the cluster should have better correlation with disease info

calculate_corr <- function(data, dist_method = "euclidean", cluster_method="complete") {
    data = melt(data[, str_detect(colnames(data), "(X\\d+|ident)")])
    data$variable = as.numeric(str_replace(data$variable, "X", ""))
    
    # data <- data %>% 
    #     group_by(ident, variable) %>%
    #     mutate(value=sum(value)) %>%
    #     unique()
    # data$value = as.numeric(data$value)
    
    data <- dcast(data, variable ~ ident, fun.aggregate = sum)
    data <- data[, colnames(data) != "variable"]
    
    data <- na.omit(t(scale(data)))
    hclust(dist(data, method = dist_method), method = cluster_method)
}



### make heatmap of different genes
cnv <- cnv[!is.na(cnv$gene_name), ]


# read_markers <- function(path) {
#     data = read.xlsx(path)
#     data = data[order(abs(data$avg_logFC), decreasing = T), ]
#     data[abs(data$avg_logFC) > 0.25 & data$p_val_adj < 0.05,]
# }
# 
# 
# atii <- read_markers("03_each_cells/total/Alveolar_II/disease/markers.xlsx")
# basal <- read_markers("03_each_cells/total/Alveolar_II/disease/markers.xlsx")
# 
# bulk <- read.xlsx("09_bulk/RNA_seq/Disease_DEGs.xlsx")
# bulk <- bulk[bulk$Type == "Not", ]
# 
# 
# genes <- intersect(atii$gene[1: 100], basal$gene[1: 100])
# 
# genes1 <- intersect(genes, bulk$symbol)

# 
# h <- readRDS("03_each_cells/select_features/ATII_importance_gene_heatmap.rds")
# temp_gene = h@ht_list$Expr@row_names_param$labels[h@ht_list$Expr@row_order]
# temp_gene_mark <- round(length(temp_gene) / 2)
# temp_gene_mark <- temp_gene[(temp_gene_mark - 120): (temp_gene_mark + 30)]
# 
# targets <- c(
#     "ALKBH3", "ALOX5AP", "BRAF", 
#     "CASP8", "CCDC130", "EGFR",
#     "EIF2S1", "ERP27", "G6PD", 
#     "INSR", "KRAS", "MET", 
#     "PSMB4", "RBPJ", "ROS1", 
#     "RPS6KA1", "SYT12", "TNFSF9",
#     "GZMK", "GZMA", "CCL4", "CXCL13",
#     "DNAJB1", "HSPB1", "CCL3L1", "CD69",
#     "MT2A", "XIST", "SPINK1"
# )
# 
# gene_drug <- read.xlsx("gene_drug.xlsx")
# gene_drug <- intersect(cnv$gene_name, gene_drug$Symbol)
# gene_drug <- gene_drug[!gene_drug %in% c("KDR", "EFCAB13")]
# 
# 
# gene_drug <- unique(c(gene_drug, targets))


bulk_meta <- read.xlsx("09_bulk/RNATAC50commonSample-1119.xlsx")
bulk_meta <- bulk_meta[!is.na(bulk_meta$SampleID), ]
rownames(bulk_meta) <- bulk_meta$SampleID
bulk_meta$Stage <- str_replace(bulk_meta$Stage, "[^IV]", "")
bulk_meta$Stage[is.na(bulk_meta$Stage)] <- "BSPN"


# gene_mark = c("ACSL4", "ALPL", "AQP4", "BMP2", "SFTPC", "C4BPA")
# basal_mark = c("VIM", "NAPSA", "SCGB3A1", "KRT5", "IGFBP2", "ALDH3A1")
# 
# 
# atii_clt_markers <- read.xlsx("03_each_cells/LUAD/Alveolar_II/markers_by_cluster.xlsx", rowNames = T)
# atii_clt_markers <- atii_clt_markers %>%
#     # filter(p_val_adj < 0.05 & avg_logFC > 0.25) %>%
#     group_by(ident) %>%
#     top_n(5, wt = avg_logFC) %>%
#     as.data.frame()
# 
# 
# basal_clt_markers <- read.xlsx("03_each_cells/LUSC/Basal/markers_by_cluster.xlsx", rowNames = T)
# basal_clt_markers <- basal_clt_markers %>%
#     # filter(p_val_adj < 0.05 & avg_logFC > 0.25) %>%
#     group_by(ident) %>%
#     top_n(5, wt = avg_logFC) %>%
#     as.data.frame()

res_markers <- readRDS("Lung_cancer_10x/ATII/new_markers_multiple_seed.rds")

scores <- read.table("09_bulk/ATAC/computedMatrix/all.res", header = F)

gene_pos <- read.table("09_bulk/RNA_seq/gene_pos.txt")
rownames(gene_pos) <- gene_pos$V5

gene_trans <- read.table("09_bulk/gene_iso.txt", header = T)
rownames(gene_trans) <- gene_trans$transcript_id
scores$gene_id <- gene_trans[as.character(scores$V2), "gene_id"]
scores$gene_name = gene_pos[as.character(scores$gene_id), "V6"]

intersect(as.character(res_markers$gene)[res_markers$ident == 'LUAD-I'], as.character(scores$gene_name[scores$V1 > 0.6]))



make_atac_heatmap_by_scores <- function(cnv, ident, title) {
    temp = cnv[as.character(cnv$trans_id) == as.character(ident), ]
    
    rownames(temp) <- temp$ident
    temp <- temp[, str_detect(colnames(temp), "^X\\d+")]
    temp <- t(na.omit(scale(t(temp))))
    
    ra <- rowAnnotation(
        Disease = sapply(rownames(temp), function(x) {
            str_replace(x, "\\d+", "")
        }),
        Stage = bulk_meta[rownames(temp), "Stage"],
        col = list(
            Disease=c(
                "LUAD" = "#0084D1", 
                "Normal"="#73BDFF", 
                "LUSC"="#A0C807", 
                "BSPN"="#6A0019", 
                "BSPN_Normal"="#C5000B"
            ),
            Stage = c(
                "I"="#65A9A3", 
                "II"="#4A933E", 
                "III"="#EC7A21", 
                "IV"="#D73F47", 
                "BSPN"="#6A0019", 
                "LUAD_Normal"="#FECC1B", 
                "LUSC_Normal"="#778793"
                
            )
        ),
        show_annotation_name = FALSE
    )
    
    draw(
        Heatmap(
            temp,
            name="Expr",
            cluster_rows = T,
            cluster_columns = F,
            right_annotation = ra,
            show_column_names = F,
            column_title = title,
            column_title_gp = gpar(fontsize = 30)
        )
    )
}


make_atac_heatmap_by_scores(cnv, "ENST00000494263", "ENST00000494263")


expr <- cnv[cnv$trans_id == "ENST00000540698",]
rownames(expr) <- expr$ident
expr <- as.matrix(expr[, str_detect(colnames(expr), "X\\d+")])
expr <- t(na.omit(scale(t(expr))))
expr <- melt(expr)
expr$Var1 <- as.character(expr$Var1)
expr$Var2 <- as.numeric(as.character(str_replace_all(expr$Var2, "^X", "")))
expr$Disease <- as.character(str_replace_all(expr$Var1, "\\d+", ""))


expr <- expr %>%
    group_by(Disease, Var2) %>%
    mutate(value = mean(value)) %>%
    dplyr::select(Var2, value, Disease) %>%
    unique()


ggplot(expr[expr$Var1 %in% c("LUSC1", "LUSC14", "LUSC3", "LUAD7", "LUAD12", "LUAD25"), ], aes(x=Var2, y = value, color = Disease)) +
    geom_line() +
    facet_grid(Var1~.)


registerDoMC(5)
foreach(i = gene_names$trans_id[gene_names$gene_name %in% c("DHX9", "NBPF1", "RBM42", "KLHDC9", "GCSH", "ZNF704", "UTP15")], .errorhandling = "pass") %dopar% {
    
    target.dir = "09_bulk/ATAC/scores/"
    dir.create(target.dir, showWarnings = F, recursive = T)
    
    gn = gene_names[i, "gene_name"]
    if (gn == "") {
        output = paste0(i, ".pdf")
        title = i
    } else {
        output = paste0(gn, "_", i, ".pdf")
        title = paste0(gn, " (", i, ")")
    }

    pdf(paste0(target.dir, output), width = 10, height = 8)
    make_atac_heatmap_by_scores(cnv, ident = i, title = title)
    dev.off()
}

# 
# res_markers <- readRDS("Lung_cancer_10x/ATII/new_markers_multiple_seed.rds")
# 
# intersect(as.character(res_markers$gene[res_markers$ident == "LUAD-I"]), as.character(scores$gene_name[scores$V1 > 0.6]))
# 
# intersect(as.character(new_markers$gene[new_markers$ident == "LUAD-I"]), as.character(scores$gene_name[scores$V1 > 0.6]))


registerDoMC(5)
foreach(i = temp_markers$gene[temp_markers$ident == "LUAD-I"], .errorhandling = "pass") %dopar% {
    
    target.dir = "09_bulk/ATAC/atii_markers/"
    dir.create(target.dir, showWarnings = F, recursive = T)
    
    print(i)
    temp = cnv[as.character(cnv$gene_name) == i, ]
    
    temp = melt(temp[, str_detect(colnames(temp), "(X\\d+|ident)")])
    temp$variable = as.numeric(str_replace(temp$variable, "X", ""))
    
    # temp <- temp %>% 
    #     group_by(ident, variable) %>%
    #     mutate(value=sum(value)) %>%
    #     unique()
    # temp$value = as.numeric(temp$value)
    
    temp <- dcast(temp, variable ~ ident, fun.aggregate = sum)
    temp <- temp[, colnames(temp) != "variable"]
    
    temp <- na.omit(t(scale(temp)))
    
    ra <- rowAnnotation(
        Disease = sapply(rownames(temp), function(x) {
            str_replace(x, "\\d+", "")
        }),
        Stage = bulk_meta[rownames(temp), "Stage"],
        col = list(
            Disease=c(
                "LUAD" = "#0084D1", 
                "Normal"="#73BDFF", 
                "LUSC"="#A0C807", 
                "BSPN"="#6A0019", 
                "BSPN_Normal"="#C5000B"
            ),
            Stage = c(
                "I"="#65A9A3", 
                "II"="#4A933E", 
                "III"="#EC7A21", 
                "IV"="#D73F47", 
                "BSPN"="#6A0019", 
                "LUAD_Normal"="#FECC1B", 
                "LUSC_Normal"="#778793"
                
            )
        ),
        show_annotation_name = FALSE
    )
    
    pdf(paste0(target.dir, i, ".pdf"), width = 10, height = 8)
    draw(
        Heatmap(
            temp,
            name="Expr",
            cluster_rows = T,
            cluster_columns = F,
            right_annotation = ra,
            column_title = i,
            column_title_gp = gpar(fontsize = 30)
        )
    )
    dev.off()
}



temp <- merge(res_markers, scores, by.x = "gene", by.y = "gene_name")


ggplot(temp, aes(x=ident, y = V1)) + geom_violin()


atii_markers <- read.xlsx("03_each_cells/total/Alveolar_II/disease/markers.xlsx", rowNames = T)

atii_markers$ident <- atii_markers$p_val_adj < 0.05 & abs(atii_markers$avg_logFC) > 0.25
atii_markers <- merge(atii_markers, scores, by.x = "gene", by.y = "gene_name")

ggplot(atii_markers, aes(x=ident, y = V1)) + geom_violin()

####



make_voca_plot <- function(data, genes.label, genes.use = NULL, title = "") {
    
    if (is.null(genes.use)) {
        genes.use = unique(data$gene)
    } 
    
    data = data[data$gene %in% genes.use, ]
    
    data$ident[
        data$avg_logFC > 0.25 & data$p_val_adj < 0.05
        ] = "LUAD"
    
    data$ident[
        data$avg_logFC < -0.25 & data$p_val_adj < 0.05
        ] = "LUSC"
    
    data$ident[
        abs(data$avg_logFC) <= 0.25 | data$p_val_adj >= 0.05
        ] = "Not" 
    
    # data$ident <- apply(data, 1, function(row) {
    #     if (row[2] > 0.25 && row[5] < 0.05) {
    #         return("1")
    #     } else if (row[2] < -0.25 && row[5] < 0.05) {
    #         return("2")
    #     } 
    
    #     return("3")
    # })
    
    data$p_val_adj <- -log10(data$p_val_adj)
    # data$gene <- rownames(data)
    data$ident = as.factor(data$ident)
    
    p <- ggplot(data = data, aes(x=avg_logFC, y=p_val_adj, color = ident)) + 
        geom_point_rast() + 
        scale_color_manual(
            values=c(
                "LUAD" = "#0084D1", 
                "LUSC"="#A0C807",
                "Not"="grey50"
            ),
            breaks = c(
                "LUAD",
                "LUSC"
            )
        ) +
        geom_text_repel(
            aes(label=ifelse(gene %in% genes.label & ident != "Not", as.character(gene),''))
        ) + theme_bw() +
        theme(
            legend.position = c(0.1, 0.95),
            axis.text = element_text(size = 20),
            axis.title = element_text(size = 20),
            legend.text = element_text(size = 15),
            legend.background = element_blank(),
            legend.key = element_blank(),
            title = element_text(size = 20)
        ) +
        labs(x="log2(Fold Change)", y="-log10(FDR)", title = title, color="") +
        guides(color = guide_legend(override.aes = list(size = 5)))
    
    return(p)
}

markers <- read.xlsx("03_each_cells/total/CD8/disease/markers.xlsx", rowNames = T)
markers$gene <- rownames(markers)
genes.label = c()
temp_markers = markers[markers$p_val_adj < 0.05, ]

temp_markers = temp_markers[order(temp_markers$avg_logFC, decreasing = T), ]
genes.label = temp_markers$gene[1:20]
temp_markers = temp_markers[order(temp_markers$avg_logFC, decreasing = F), ]
genes.label = c(genes.label, temp_markers$gene[1:20])


p <- make_voca_plot(
    markers, 
    genes.label = c(genes.label, "GNLY")
) +
    theme(legend.position = c(0.8, 0.9))

ggsave(
    filename = "03_each_cells/total/CD8/disease/volcano.pdf",
    plot = p,
    width = 6,
    height = 6
)



markers <- read.xlsx("03_each_cells/total/CD4/disease/markers.xlsx", rowNames = T)
markers$gene <- rownames(markers)
genes.label = c()
temp_markers = markers[markers$p_val_adj < 0.05, ]

temp_markers = temp_markers[order(temp_markers$avg_logFC, decreasing = T), ]
genes.label = temp_markers$gene[1:15]
temp_markers = temp_markers[order(temp_markers$avg_logFC, decreasing = F), ]
genes.label = c(genes.label, temp_markers$gene[1:15])


p <- make_voca_plot(
    markers, 
    genes.label = c(genes.label, "GNLY", "DNAJB1")
) +
    theme(legend.position = c(0.8, 0.9))

ggsave(
    filename = "03_each_cells/total/CD4/disease/volcano.pdf",
    plot = p,
    width = 6,
    height = 6
)


markers <- read.xlsx("03_each_cells/total/B/disease/markers.xlsx", rowNames = T)
markers$gene <- rownames(markers)
genes.label = c()
temp_markers = markers[markers$p_val_adj < 0.05, ]

temp_markers = temp_markers[order(temp_markers$avg_logFC, decreasing = T), ]
genes.label = temp_markers$gene[1:15]
temp_markers = temp_markers[order(temp_markers$avg_logFC, decreasing = F), ]
genes.label = c(genes.label, temp_markers$gene[1:15])


p <- make_voca_plot(
    markers, 
    genes.label = c(genes.label, "GNLY", "DNAJB1")
) +
    theme(legend.position = c(0.8, 0.9))

ggsave(
    filename = "03_each_cells/total/CD4/disease/volcano.pdf",
    plot = p,
    width = 6,
    height = 6
)