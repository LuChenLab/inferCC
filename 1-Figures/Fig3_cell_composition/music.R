library(MuSiC)
library(xbioc)
library(dplyr)
library(openxlsx)
library(stringr)
library(ComplexHeatmap)
library(circlize)
library("FactoMineR")
library("factoextra")
library(ggfortify)
library(umap)
library(Rtsne)
library(reshape2)

set.seed(1)
options(stringsAsFactor = F)


make_heatmap_music <- function(data, scale_dot_size=0.9) {
    col_fun <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
    
    colors = list(
        Stage = c("I"="#65A9A3", "II"="#4A933E", "III"="#EC7A21", "IV"="#D73F47", "LUAD_Normal"="#FECC1B", "LUSC_Normal"="#778793", "BSPN"="#73BDFF"),
        Disease = c("LUAD" = "#0084D1", "NL(AD)"="#FECC1B", "NL(AD"="#778793", "BSPN"="#73BDFF", "NL(SC)"="#778793", "LUSC"="#A0C807"),
        Batch = c("1"="#E9B0B7", "2"="#90D0E0")
    )
    
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
    
    Heatmap(
        data,
        name = "Corr",
        col = col_fun,
        right_annotation = ra,
        bottom_annotation = ba,
        show_row_names = F,
        show_column_names = F,
        rect_gp = gpar(type = "none"),
        cell_fun = function(j, i, x, y, width, height, fill) {
            s = min(unit.c(convertWidth(width, "cm"), convertHeight(height, "cm")))
            grid.circle(
                x = x,
                y = y, 
                r = s * (abs(data[i, j]) + 0.1) * scale_dot_size,
                gp = gpar(col = "white", fill = fill)
            )
        }
    )
}


make_heatmap_tcga <- function(
    data,
    cluster_method = "ward.D2",
    border = TRUE,
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
    
    col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
    
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
    )
}

axis_size = 18
axis_title_size = 20
title_size = 25
legend.position="right"
legend_text_size = 15
legend_title_size = 20

my_theme <- theme(
    aspect.ratio=1,
    plot.title = element_text(hjust = 0.5, face="bold", size = title_size),
    legend.position = legend.position,
    axis.text = element_text(size = axis_size),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = axis_size),
    axis.title = element_text(size = axis_title_size),
    legend.text = element_text(size = legend_text_size),
    legend.title = element_text(size = legend_title_size)
)


make_pca_analysis <- function(data, output_prefix) {
    res.pca <- prcomp(data)
    
    p <- fviz_eig(
        res.pca, 
        ggtheme = theme_classic() + my_theme
    )
    ggsave(
        filename = paste(output_prefix, "scree_plot.pdf", sep = "_"),
        plot = p,
        device = "pdf",
        width = 6,
        height = 6
    )
    
    ggsave(
        filename = paste(output_prefix, "contrib_dim1.pdf", sep = "_"),
        device = "pdf",
        plot = fviz_contrib(
            res.pca, 
            choice = "var", 
            axes = 1, 
            top = 10, 
            ggtheme = theme_classic() + theme(
                plot.title = element_text(hjust = 0.5, face="bold", size = title_size),
                legend.position = legend.position,
                axis.text = element_text(size = axis_size),
                axis.text.x = element_text(size = axis_size),
                axis.title = element_text(size = axis_title_size),
                legend.text = element_text(size = legend_text_size),
                legend.title = element_text(size = legend_title_size)
            )
        ) ,
        width = 8,
        height = 6
    )
    
    
    ggsave(
        filename = paste(output_prefix, "contrib_dim2.pdf", sep = "_"),
        plot = fviz_contrib(
            res.pca, 
            choice = "var",
            axes = 2, 
            top = 10,
            ggtheme = theme_classic() + theme(
                plot.title = element_text(hjust = 0.5, face="bold", size = title_size),
                legend.position = legend.position,
                axis.text = element_text(size = axis_size),
                axis.text.x = element_text(size = axis_size),
                axis.title = element_text(size = axis_title_size),
                legend.text = element_text(size = legend_text_size),
                legend.title = element_text(size = legend_title_size)
            )
        ),
        device = "pdf",
        width = 8,
        height = 6
    )
    
    
    p <- fviz_pca_var(
        res.pca,
        col.var = "contrib", # Color by contributions to the PC
        gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
        repel = TRUE,     # Avoid text overlapping
        ggtheme = theme_classic() + my_theme,
        labelsize = 8
    )
    ggsave(
        filename = paste(output_prefix, "contrib_pca.pdf", sep = "_"),
        device = "pdf",
        plot = p,
        width = 6,
        height = 6
    )
}

#### Function to make umap of pca 
make_umap <- function(x, dot.size = 2) {
    data = as.data.frame(x$layout[, 1:2])
    
    colnames(data) <- c("UMAP1", "UMAP2")

    data$Disease <- sapply(rownames(data), function(x) {
        if (str_detect(x, "^TCGA")) {
            x = str_split(x, "\\.")[[1]][8]
        } else {
            x = str_replace_all(x, "\\d+", "")
        }
        
        return(x)
    })
    data$Stage <- sapply(rownames(data), function(x) {
        if (str_detect(x, "^TCGA")) {
            x = str_split(x, "\\.")[[1]][9]
        } else {
            x = bulk_meta[x, "Stage"]
        }
        
        return(x)
    })

    p1 <- ggplot(data, aes(x = UMAP1, y = UMAP2, color  = Disease)) +
        geom_point(size = dot.size) +
        theme_bw() +
        my_theme + 
        theme(legend.position = "bottom") +
        labs(color = "") +
        scale_color_manual(values = c(
            "LUAD" = "#0084D1", 
            "Normal (LUAD)"="#FECC1B", 
            "NL(AD"="#778793", 
            "BSPN"="#6A0019", 
            "Normal (LUSC)"="#778793", 
            "LUSC"="#A0C807",
            "LC"="#91822B",
            "NL(LC)"="#EDCAB0",
            "UPS"="#EBAAA4",
            "NL(UPS)"="#B43018"
        )) +
        guides(
            color = guide_legend(override.aes = list(size = 4))
        )
    
    p2 <- ggplot(data, aes(x = UMAP1, y = UMAP2, color  = Stage)) +
        geom_point(size = dot.size) +
        theme_bw() +
        my_theme + 
        theme(legend.position = "bottom") +
        labs(color = "") +
        scale_color_manual(values = c(
            "I"="#65A9A3", "II"="#4A933E", 
            "III"="#EC7A21", "IV"="#D73F47", 
            "LUAD_Normal"="#FECC1B", 
            "LUSC_Normal"="#778793", 
            "BSPN"="#73BDFF"
        )) +
        guides(
            color = guide_legend(override.aes = list(size = 4))
        )
    return (cowplot::plot_grid(p1, p2, ncol = 2))
}




make_tsne <- function(data, dot.size = 2) {
    
    colnames(data) <- c("tSNE_1", "tSNE_2")
    
    data$Disease <- sapply(rownames(data), function(x) {
        if (str_detect(x, "^TCGA")) {
            x = str_split(x, "\\.")[[1]][8]
        } else {
            x = str_replace_all(x, "\\d+", "")
        }
        
        return(x)
    })
    data$Stage <- sapply(rownames(data), function(x) {
        if (str_detect(x, "^TCGA")) {
            x = str_split(x, "\\.")[[1]][9]
        } else {
            x = bulk_meta[x, "Stage"]
        }
        
        return(x)
    })
    
    p1 <- ggplot(data, aes(x = tSNE_1, y = tSNE_2, color  = Disease)) +
        geom_point(size = dot.size) +
        theme_bw() +
        my_theme + 
        theme(legend.position = "bottom") +
        labs(color = "", title = "Disease") +
        scale_color_manual(values = c(
            "LUAD" = "#0084D1", 
            "Normal (LUAD)"="#FECC1B", 
            "NL(AD"="#778793", 
            "Normal"="#73BDFF", 
            "Normal (LUSC)"="#778793", 
            "LUSC"="#A0C807",
            "LC"="#91822B",
            "NL(LC)"="#EDCAB0",
            "UPS"="#EBAAA4",
            "NL(UPS)"="#B43018"
        )) +
        guides(
            color = guide_legend(override.aes = list(size = 4))
        )
    
    p2 <- ggplot(data, aes(x = tSNE_1, y = tSNE_2, color = Stage)) +
        geom_point(size = dot.size) +
        theme_bw() +
        my_theme + 
        theme(legend.position = "bottom") +
        labs(color = "", title = "Stage") +
        scale_color_manual(values = c(
            "I"="#65A9A3", "II"="#4A933E", 
            "III"="#EC7A21", "IV"="#D73F47", 
            "LUAD_Normal"="#FECC1B", 
            "LUSC_Normal"="#778793", 
            "BSPN"="#73BDFF"
        )) +
        guides(
            color = guide_legend(override.aes = list(size = 4))
        )

    return (cowplot::plot_grid(p1, p2, ncol = 2))
}


make_pca_reduction_plots <- function(data, output_prefix, n.pcs = 10, dot.size = 2) {
    library(umap)
    
    res.pca = res.pca <- prcomp(data)
    
    p <- ggplot(
        data.frame(x = 1:ncol(res.pca$x), y = apply(res.pca$x, 2, sd)), 
        aes(x = x, y = y)
    ) + geom_point()
    
    print(p)
    
    colors = data.frame(
        Disease = sapply(rownames(data), function(x) {
            if (str_detect(x, "^TCGA")) {
                x = str_split(x, "\\.")[[1]][8]
            } else {
                x = str_replace_all(x, "\\d+", "")
            }
            
            return(x)
        }),
        Stage = sapply(rownames(data), function(x) {
            if (str_detect(x, "^TCGA")) {
                x = str_split(x, "\\.")[[1]][9]
            } else {
                x = bulk_meta[x, "Stage"]
            }
            
            return(x)
        })
    )
    
    p1 <- autoplot(res.pca, data = colors,  colour = 'Disease', frame = TRUE) + 
        my_theme + theme(aspect.ratio=1, legend.position = "bottom") +
        labs(title = "Disease") +
        scale_color_manual(values = c(
            "LUAD" = "#0084D1", 
            "Normal (LUAD)"="#FECC1B", 
            "NL(AD"="#778793", 
            "Normal"="#73BDFF", 
            "Normal (LUSC)"="#778793", 
            "LUSC"="#A0C807",
            "LC"="#91822B",
            "NL(LC)"="#EDCAB0",
            "UPS"="#EBAAA4",
            "NL(UPS)"="#B43018"
        )) +
        guides(
            color = guide_legend(override.aes = list(size = 4))
        )
    
    p2 <- autoplot(res.pca, data = colors,  colour = 'Stage', frame = TRUE) + 
        my_theme + theme(aspect.ratio=1, legend.position = "bottom")  +
        labs(title = "Stage") +
        scale_color_manual(values = c(
            "I"="#65A9A3", "II"="#4A933E", 
            "III"="#EC7A21", "IV"="#D73F47", 
            "LUAD_Normal"="#FECC1B", 
            "LUSC_Normal"="#778793", 
            "BSPN"="#73BDFF"
        )) +
        guides(
            color = guide_legend(override.aes = list(size = 4))
        )
    
    ggsave(
        filename = paste(output_prefix, "pca.pdf", sep = "_"),
        plot = cowplot::plot_grid(p1, p2, ncol = 2),
        device = "pdf",
        width = 12,
        height = 6
    )
    
    res.umap <- umap(res.pca$x[, 1:n.pcs])
    ggsave(
        filename = paste(output_prefix, "umap.pdf", sep = "_"),
        plot = make_umap(res.umap, dot.size),
        device = "pdf",
        width = 12,
        height = 6
    )
    
    res.tsne <- Rtsne(
        res.pca$x[, 1:n.pcs], 
        dims = 2, 
        perplexity=min(30, n.pcs - 1), 
        verbose=F, 
        max_iter = 500,
        check_duplicates = F
    )
    
    data = as.data.frame(res.tsne$Y)
    rownames(data) <- rownames(res.pca$x)
    ggsave(
        filename = paste(output_prefix, "tsne.pdf", sep = "_"),
        plot = make_tsne(data, dot.size),
        device = "pdf",
        width = 12,
        height = 6
    )
}


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



### read bulk
bulk <- read.xlsx("09_bulk/RNA_seq/RSEM.xlsx", rowNames = T)
gene_pos <- read.table("09_bulk/RNA_seq/gene_pos.txt", header = F, row.names = 5)
gene_pos$V6 <- make.unique(as.character(gene_pos$V6))
rownames(bulk) <- gene_pos[rownames(bulk), "V6"]


bulk_meta <- read.xlsx("09_bulk/RNATAC50commonSample-1119.xlsx")
bulk_meta <- bulk_meta[!is.na(bulk_meta$SampleID), ]
rownames(bulk_meta) <- bulk_meta$SampleID
bulk_meta$Stage <- str_replace(bulk_meta$Stage, "[^IV]", "")

bulk_set <- ExpressionSet(as.matrix(bulk), phenoData = AnnotatedDataFrame(bulk_meta[colnames(bulk), ]))



### read single cell
meta <- read.csv("02_rds/meta_after_singleR.csv", row.names = 1, stringsAsFactors = F)
meta <- meta[meta$Batch != 3 & meta$cell_name != "", ]
meta$cell_short[meta$cell_short == "Mφ"] <- "Mo"

expr <- readRDS("02_rds/all_cell_expr.rds")

sc_set <- ExpressionSet(as.matrix(expr[, rownames(meta)]), phenoData = AnnotatedDataFrame(meta))



dec = music_prop(
    bulk.eset = bulk_set, 
    sc.eset = sc_set, 
    clusters = 'cell_short',
    samples = 'SampleID', 
    markers = markers$Markers,
    verbose = F
)

saveRDS(dec, "09_bulk/RNA_seq/MuSic/bulk.rds")




## TCGA
tcga <- read.csv("TCGA/Lung.csv", row.names = 1, stringsAsFactors = F)


### read bulk
tcga_meta <- data.frame(
    Disease = sapply(colnames(tcga), function(x) { str_split(x, "\\.")[[1]][8] }),
    Stage = sapply(colnames(tcga), function(x) { str_split(x, "\\.")[[1]][9] })
)

rownames(tcga_meta) <- colnames(tcga)
tcga_meta <- tcga_meta[tcga_meta$Stage != "", ]


bulk_set <- ExpressionSet(
    as.matrix(tcga[, rownames(tcga_meta)]), 
    phenoData = AnnotatedDataFrame(tcga_meta)
)


dec = music_prop(
    bulk.eset = bulk_set, 
    sc.eset = sc_set, 
    clusters = 'cell_short',
    samples = 'SampleID', 
    markers = markers$Markers,
    verbose = T
)


saveRDS(dec, "09_bulk/RNA_seq/MuSic/tcga.rds")


## East asian

### East Asian
read_east_asian <- function() {
    tumor <- read.table("09_bulk/East_Asians/GIS031/GSK_RSEM_expCounts/GSK_RSEM_rerun_expCounts_172_Tumor.tsv", header = T)
    tumor <- tumor[!is.na(tumor$Gene.symbol), ]
    rownames(tumor) <- make.unique(as.character(tumor$Gene.symbol))
    tumor <- tumor[, str_detect(colnames(tumor), "A\\d+")]
    colnames(tumor) <- paste(colnames(tumor), "Tumor", sep = "_")
    
    east <- read.table("09_bulk/East_Asians/GIS031/GSK_RSEM_expCounts/GSK_RSEM_rerun_expCounts_88_Normal.tsv", header = T)
    east <- east[!is.na(east$Gene.symbol), ]
    rownames(east) <- make.unique(as.character(east$Gene.symbol))
    east <- east[, str_detect(colnames(east), "A\\d+")]
    colnames(east) <- paste(colnames(east), "Normal", sep = "_")
    
    for (i in colnames(tumor)) {
        east[, i] <- tumor[rownames(east), i]
    }
    
    east
}


read_east_meta <- function() {
    east_meta <- read.table("09_bulk/East_Asians/GIS031/GIS031.clinical.patient.txt", header = T, sep = "\t")
    temp <- east_meta
    temp$PATIENT_ID <- paste(temp$PATIENT_ID, "Tumor", sep = "_")
    
    east_meta$PATIENT_ID <- paste(east_meta$PATIENT_ID, "Normal", sep = "_")
    
    east_meta <- rbind(east_meta, temp)
    rownames(east_meta) <- east_meta$PATIENT_ID
    east_meta
}

east <- read_east_asian()
east_meta <- read_east_meta()

east_set <- ExpressionSet(as.matrix(east), phenoData = AnnotatedDataFrame(east_meta[colnames(east), ]))


dec = music_prop(
    bulk.eset = east_set, 
    sc.eset = sc_set, 
    clusters = 'cell_short',
    samples = 'SampleID', 
    markers = markers$Markers,
    verbose = F
)


saveRDS(dec, "09_bulk/RNA_seq/MuSic/east.rds")




## GTEx

gtex <- read.table("09_bulk/GTEx/All_Tissue_Site_Details.combined.reads.gct", sep = "\t", header = T)
row.names(gtex) <- make.unique(as.character(gtex$Description))

gtex_meta <- read.table("09_bulk/GTEx/meta.txt", sep = "\t")
rownames(gtex_meta) <- gtex_meta$V1


sel_cells = sapply(gtex_meta$V1, function(x) {
    str_replace_all(x, "-", ".")
})

gtex <- gtex[, intersect(sel_cells, colnames(gtex))]
colnames(gtex) <- sapply(colnames(gtex), function(x) {
    str_replace_all(x, "\\.", "-")
})

saveRDS(gtex, "09_bulk/GTEx/gtex_lung.rds")


gtex_set <- ExpressionSet(as.matrix(gtex), phenoData = AnnotatedDataFrame(gtex_meta[colnames(gtex), ]))


dec = music_prop(
    bulk.eset = gtex_set, 
    sc.eset = sc_set, 
    clusters = 'cell_short',
    samples = 'SampleID', 
    markers = markers$Markers,
    verbose = F
)

saveRDS(dec, "09_bulk/RNA_seq/MuSic/gtex.rds")



# Plots
## Bulk

dec <- readRDS("09_bulk/RNA_seq/MuSic/bulk.rds")
data <- cor(t(dec$Est.prop.weighted))


est <- dec$Est.prop.weighted
est <- est[str_detect(rownames(est), "^LU"), ]

pdf("09_bulk/RNA_seq/MuSic/wo_BSPN_heatmap_all_sample.pdf", width = 5, height = 4)
make_heatmap_music(cor(t(est)), scale_dot_size = 0.4)
dev.off()

pdf("09_bulk/RNA_seq/MuSic/wo_BSPN_heatmap_all_sample_spearman.pdf", width = 5, height = 4)
make_heatmap_music(cor(t(est), method = "spearman"), scale_dot_size = 0.4)
dev.off()



make_pca_analysis(est, "02_rds/pca/Bulk/all_sample")

make_pca_reduction_plots(est, "02_rds/pca/Bulk/all_sample", n.pcs = 5, dot.size = 4)


###
immu <- c(
    "B",
    "CD4",
    "CD8",
    "DC",
    "Mast",
    "Mo",
    "NK",
    "Tregs",
    "Gran"
)


temp <- melt(as.matrix(est))
temp$Disease <- str_replace_all(temp$Var1, "\\d+", "")
temp$Type <- ifelse(temp$Var2 %in% immu, "Immune", "Non-immune") 
temp$Var2 <- factor(temp$Var2, levels = sort(unique(as.character(temp$Var2))))

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
    scale_color_manual(values =c(
        "LUAD" = "#0084D1", 
        "NL(AD)"="#FECC1B", 
        "NL(AD"="#778793", 
        "BSPN"="#FB290F", 
        "NL(SC)"="#778793", 
        "LUSC"="#A0C807")
    )


ggsave(
    filename = "09_bulk/RNA_seq/MuSic/w0_BSPN_boxplot.pdf",
    plot = p,
    width = 12,
    height = 8
)



est <- dec$Est.prop.weighted
est <- est[str_detect(rownames(est), "^LU"), ]


pdf("09_bulk/RNA_seq/MuSic/wo_BSPN_heatmap_immune.pdf", width = 5, height = 4)
temp_est <- cor(t(est[, colnames(est) %in% immu]))
temp_est[is.na(temp_est)] <- 0
make_heatmap_music(temp_est, scale_dot_size = 0.4)
dev.off()


pdf("09_bulk/RNA_seq/MuSic/wo_BSPN_heatmap_immune_spearman.pdf", width = 5, height = 4)
temp_est <- cor(t(est[, colnames(est) %in% immu]), method = "spearman")
temp_est[is.na(temp_est)] <- 0
make_heatmap_music(temp_est, scale_dot_size = 0.4)
dev.off()


pdf("09_bulk/RNA_seq/MuSic/wo_BSPN_heatmap_non_immune.pdf", width = 5, height = 4)
make_heatmap_music(cor(t(est[, !colnames(est) %in% immu])), scale_dot_size = 0.4)
dev.off()


pdf("09_bulk/RNA_seq/MuSic/wo_BSPN_heatmap_non_immune_spearman.pdf", width = 5, height = 4)
make_heatmap_music(cor(t(est[, !colnames(est) %in% immu]), method = "spearman"), scale_dot_size = 0.4)
dev.off()



### Stage
for (i in unique(bulk_meta$Stage)) {
    
    if (!is.na(i)) {
        print(i)
        
        pdf(paste0("09_bulk/RNA_seq/MuSic/heatmap_", i, ".pdf"), width = 10, height = 8)
        draw(
            make_heatmap_music(
                cor(t(est[intersect(rownames(est), bulk_meta$SampleID[bulk_meta$Stage == i]), ])), 
                scale_dot_size = 0.4
            )
        )
        dev.off()
        
        
        pdf(paste0("09_bulk/RNA_seq/MuSic/heatmap_", i, "_non_immune.pdf"), width = 10, height = 8)
        draw(
            make_heatmap_music(
                cor(t(est[intersect(rownames(est), bulk_meta$SampleID[bulk_meta$Stage == i]), !colnames(est) %in% immu])), 
                scale_dot_size = 0.4
            )
        )
        dev.off()
        
        
        pdf(paste0("09_bulk/RNA_seq/MuSic/heatmap_", i, "_immune.pdf"), width = 10, height = 8)
        draw(
            make_heatmap_music(
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

# est <- dec$Est.prop.weighted
# est <- dec$Est.prop.allgene
# est <- dec$Var.prop

for (i in unique(bulk_meta$Disease)) {
    
    if (!is.na(i)) {
        print(i)
        
        # pearson
        pdf(paste0("09_bulk/RNA_seq/MuSic/heatmap_", i, ".pdf"), width = 10, height = 8)
        draw(
            make_heatmap_music(
                cor(t(est[intersect(rownames(est), bulk_meta$SampleID[bulk_meta$Disease == i]), ])), 
                scale_dot_size = 0.4
            )
        )
        dev.off()
        
        
        pdf(paste0("09_bulk/RNA_seq/MuSic/heatmap_", i, "_non_immune.pdf"), width = 10, height = 8)
        draw(
            make_heatmap_music(
                cor(t(est[intersect(rownames(est), bulk_meta$SampleID[bulk_meta$Disease == i]), !colnames(est) %in% immu])), 
                scale_dot_size = 0.4
            )
        )
        dev.off()
        
        
        pdf(paste0("09_bulk/RNA_seq/MuSic/heatmap_", i, "_immune.pdf"), width = 10, height = 8)
        temp_est <- cor(t(est[intersect(rownames(est), bulk_meta$SampleID[bulk_meta$Disease == i]), colnames(est) %in% immu]))
        temp_est[is.na(temp_est)] <- 0
        draw(
            make_heatmap_music(
                temp_est,
                scale_dot_size = 0.4
            )
        )
        dev.off()
        
        # spearman
        pdf(paste0("09_bulk/RNA_seq/MuSic/heatmap_", i, "_spearman.pdf"), width = 10, height = 8)
        draw(
            make_heatmap_music(
                cor(t(est[intersect(rownames(est), bulk_meta$SampleID[bulk_meta$Disease == i]), ]), method = "spearman"), 
                scale_dot_size = 0.4
            )
        )
        dev.off()
        
        
        pdf(paste0("09_bulk/RNA_seq/MuSic/heatmap_", i, "_non_immune_spearman.pdf"), width = 10, height = 8)
        draw(
            make_heatmap_music(
                cor(t(est[intersect(rownames(est), bulk_meta$SampleID[bulk_meta$Disease == i]), !colnames(est) %in% immu]), method = "spearman"), 
                scale_dot_size = 0.4
            )
        )
        dev.off()
        
        
        pdf(paste0("09_bulk/RNA_seq/MuSic/heatmap_", i, "_immune_spearman.pdf"), width = 10, height = 8)
        temp_est <- cor(t(est[intersect(rownames(est), bulk_meta$SampleID[bulk_meta$Disease == i]), colnames(est) %in% immu]), method = "spearman")
        temp_est[is.na(temp_est)] <- 0
        draw(
            make_heatmap_music(
                temp_est,
                scale_dot_size = 0.4
            )
        )
        dev.off()
    }
}



## TCGA

dec <- readRDS("09_bulk/RNA_seq/MuSic/tcga.rds")

est <- dec$Est.prop.weighted


pdf("09_bulk/RNA_seq/MuSic/tcga_heatmap_all_sample.pdf", width = 5, height = 4)
make_heatmap_tcga(cor(t(est)))
dev.off()


pdf("09_bulk/RNA_seq/MuSic/tcga_heatmap_all_sample_spearman.pdf", width = 5, height = 4)
make_heatmap_tcga(cor(t(est), method = "spearman"))
dev.off()


make_pca_analysis(est, "02_rds/pca/TCGA/all_sample")

make_pca_reduction_plots(est, "02_rds/pca/TCGA/all_sample", n.pcs = 5, dot.size = 4)



immu <- c(
    "B",
    "CD4",
    "CD8",
    "DC",
    "Mast",
    "Mo",
    "NK",
    "Tregs",
    "Granu"
)


# dec <- readRDS("09_bulk/RNA_seq/MuSic/bulk.rds")

temp <- melt(as.matrix(dec$Est.prop.weighted))
temp$Disease <- sapply(temp$Var1, function(x) { str_split(x, "\\.")[[1]][8] })
temp$Type <- ifelse(temp$Var2 %in% immu, "Immune", "Non-immune") 

temp$Var2 <- factor(temp$Var2, levels = sort(unique(as.character(temp$Var2))))

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
    scale_color_manual(values =c(
        "LUAD" = "#0084D1", 
        "NL(AD)"="#FECC1B", 
        "NL(AD"="#778793", 
        "BSPN"="#FB290F", 
        "NL(SC)"="#778793", 
        "LUSC"="#A0C807")
    )


ggsave(
    filename = "09_bulk/RNA_seq/MuSic/tcga_boxplot.pdf",
    plot = p,
    width = 12,
    height = 8
)


pdf("09_bulk/RNA_seq/MuSic/tcga_heatmap_immune.pdf", width = 5, height = 4)
temp_est <- cor(t(est[, colnames(est) %in% immu]))
temp_est[is.na(temp_est)] <- 0
make_heatmap_tcga(temp_est)
dev.off()


pdf("09_bulk/RNA_seq/MuSic/tcga_heatmap_immune_spearman.pdf", width = 5, height = 4)
temp_est <- cor(t(est[, colnames(est) %in% immu]), method = "spearman")
temp_est[is.na(temp_est)] <- 0
make_heatmap_tcga(temp_est)
dev.off()


pdf("09_bulk/RNA_seq/MuSic/tcga_heatmap_non_immune.pdf", width = 5, height = 4)
make_heatmap_tcga(cor(t(est[, !colnames(est) %in% immu])))
dev.off()


pdf("09_bulk/RNA_seq/MuSic/tcga_heatmap_non_immune_spearman.pdf", width = 5, height = 4)
make_heatmap_tcga(cor(t(est[, !colnames(est) %in% immu]), method = "spearman"))
dev.off()


set.seed(1)
for (i in unique(tcga_meta$Disease)) {
    temp_tcga_meta <- tcga_meta[tcga_meta$Disease == i, ]
    temp_tcga_meta$id <- rownames(temp_tcga_meta)
    
    temp_tcga_meta <- temp_tcga_meta %>%
        group_by(Stage) %>%
        sample_n(min(table(temp_tcga_meta$Stage)))

    
    if (!is.na(i)) {
        print(i)

        pdf(paste0("09_bulk/RNA_seq/MuSic/tcga_heatmap_", i, ".pdf"), width = 10, height = 8)
        draw(
            make_heatmap_tcga(
                cor(t(est[intersect(rownames(est), temp_tcga_meta$id), ])),
                scale_dot_size = 0.4
            )
        )
        dev.off()

        print("non-immune")
        pdf(paste0("09_bulk/RNA_seq/MuSic/tcga_heatmap_", i, "_non_immune.pdf"), width = 10, height = 8)
        temp_est <- cor(t(est[intersect(rownames(est), temp_tcga_meta$id), !colnames(est) %in% immu]))
        temp_est[is.na(temp_est)] <- 0
        draw(
            make_heatmap_tcga(
                temp_est,
                scale_dot_size = 0.4
            )
        )
        dev.off()

        print("immune")
        pdf(paste0("09_bulk/RNA_seq/MuSic/heatmap_", i, "_immune.pdf"), width = 10, height = 8)
        temp_est <- cor(t(est[intersect(rownames(est), temp_tcga_meta$id), colnames(est) %in% immu]))
        temp_est[is.na(temp_est)] <- 0
        draw(
            make_heatmap_tcga(
                temp_est,
                scale_dot_size = 0.4
            )
        )
        dev.off()
    }
}



sample_meta = data.frame(
    Disease = sapply(rownames(est), function(x) {
        str_split(x, "\\.")[[1]][8]
    }),
    Stage = sapply(rownames(est), function(x) {
        str_split(x, "\\.")[[1]][9]
    }),
    name = rownames(est),
    row.names = rownames(est)
)


for(i in 1:20) {
    set.seed(i)
    
    temp_meta <- sample_meta %>%
        group_by(Disease, Stage) %>%
        sample_n(8) %>%
        as.data.frame()
    
    pdf(paste0("09_bulk/RNA_seq/MuSic/TCGA_sample/", i,"_non_immune_spearman.pdf"), width = 10, height = 8)
    draw( make_heatmap_tcga(cor(t(est[as.character(temp_meta$name), !colnames(est) %in% immu]), method = "spearman"), scale_dot_size = 0.4))
    dev.off()
    
    pdf(paste0("09_bulk/RNA_seq/MuSic/TCGA_sample/", i,"_non_immune.pdf"), width = 10, height = 8)
    draw(make_heatmap_tcga(cor(t(est[as.character(temp_meta$name), !colnames(est) %in% immu])), scale_dot_size = 0.4))
    dev.off()
}

set.seed(1)


#### Test gene module on deconvolution data
library(reshape2)
expr <- readRDS("02_rds/all_cell_expr.rds")

meta <- read.csv("02_rds/meta_after_singleR.csv", row.names = 1, stringsAsFactors = F)
meta$cell_short = meta$cell_name1
meta <- meta[meta$cell_short != "" & !is.na(meta$cell_short), ]
meta$cell_short[meta$cell_short == "Mφ"] = "Mo"

set.seed(1)
meta$Cells = rownames(meta)

cells <- meta %>%
    sample_n(3000) %>%
    dplyr::select(Cells, cell_short) %>%
    as.data.frame()


## prepare gene expression matrix
temp <- melt(as.matrix(expr[,  which(colnames(expr) %in% cells$Cells)]))

temp$Disease = as.character(meta[as.character(temp$Var2), "Disease"])
temp$Cell = as.character(meta[as.character(temp$Var2), "cell_short"])
temp$Stage = as.character(meta[as.character(temp$Var2), "Stage"])

temp <- temp[temp$value > 0, ]

## normalize gene expression by cell type, disease and stage
temp <- temp %>%
    group_by(Disease, Cell, Stage) %>%
    add_tally() %>%
    group_by(Disease, Cell, Stage, Var1) %>%
    mutate(value = sum(value) / n)

temp <- temp[, colnames(temp) != "Var2"]
temp <- unique(temp)

## calculate gene expression composition
temp <- temp %>%
    group_by(Var1, Disease, Stage) %>%
    mutate(perc = value / sum(value)) %>%
    as.data.frame()


## Read data
bulk <- read.xlsx("09_bulk/RNA_seq/RSEM.xlsx", rowNames = T)
gene_pos <- read.table("09_bulk/RNA_seq/gene_pos.txt", header = F, row.names = 5)
gene_pos$V6 <- make.unique(as.character(gene_pos$V6))
rownames(bulk) <- gene_pos[rownames(bulk), "V6"]

bulk_meta <- read.xlsx("09_bulk/RNATAC50commonSample-1119.xlsx")
bulk_meta <- bulk_meta[!is.na(bulk_meta$SampleID), ]
rownames(bulk_meta) <- bulk_meta$SampleID
bulk_meta$Stage <- sapply(bulk_meta$Stage, function(x) {
    str_replace_all(x, "[^IV]", "")
})


dec <- readRDS("09_bulk/RNA_seq/MuSic/bulk.rds")
est <- dec$Est.prop.weighted


atii <- read.csv("03_each_cells/ATII_Basal/ATII_mfuzz.csv", row.names = 1, stringsAsFactors = F)
basal <- read.csv("03_each_cells/ATII_Basal/Basal_mfuzz.csv", row.names = 1, stringsAsFactors = F)

## prepare bulk expression
prepare_bulk_expression_matrix <- function(bulk, genes, disease, cell) {
    genes = as.character(genes)
    temp_expr <- bulk[genes, str_detect(colnames(bulk), disease)]
    temp_expr <- melt(as.matrix(temp_expr))
    temp_expr$Stage <- bulk_meta[as.character(temp_expr$Var2), "Stage"]
    temp_expr$Cell_perc <- est[as.character(temp_expr$Var2), cell]
    
    
    temp_perc = temp[as.character(temp$Var1) %in% genes, ]
    temp_perc = temp_perc[as.character(temp_perc$Disease) == disease & as.character(temp_perc$Cell) == cell, ]
    temp_perc = temp_perc[!is.na(temp_perc$Var1), ]
    rownames(temp_perc) <- paste(temp_perc$Var1, temp_perc$Stage)
    temp_expr$Gene_perc <- temp_perc[paste(temp_expr$Var1, temp_expr$Stage), "perc"]
    
    temp_expr$value <- temp_expr$value / temp_expr$Cell_perc / temp_expr$Gene_perc
    
    temp_expr <- dcast(temp_expr, Var1~Var2, fun.aggregate = mean, value.var = "value")
    temp_expr <- temp_expr[!is.na(temp_expr$Var1), ]
    rownames(temp_expr) <- as.character(temp_expr$Var1)
    
    temp_expr[is.na(temp_expr)] <- 0
    rownames(temp_expr) <- as.character(temp_expr$Var1)
    
    temp_expr[, colnames(temp_expr) != "Var1"]
}


temp_expr <- prepare_bulk_expression_matrix(bulk, basal$gene, "LUSC", "Basal")
temp_expr <- temp_expr[, apply(temp_expr, 2, function(col) {
    sum(is.infinite(col)) == 0
})]


pdf("03_each_cells/ATII_Basal/Basal_mfuzz_dec.pdf", width = 6, height = 4)
Heatmap(
    t(scale(t(temp_expr[as.character(basal$gene), ]))), 
    name = "Expr",
    col = colorRamp2(c(-1, 0, 1), c("purple", "black", "yellow")),
    show_row_names = F,
    cluster_rows = F,
    row_split = basal$Clt,
    column_split = bulk_meta[colnames(temp_expr), "Stage"],
    border = T
)
dev.off()


temp_expr <- prepare_bulk_expression_matrix(bulk, atii$gene, "LUAD", "ATII")
temp_expr <- temp_expr[, apply(temp_expr, 2, function(col) {
    sum(is.infinite(col)) == 0 && sd(col) > 0
})]


Heatmap(
    t(scale(t(temp_expr[as.character(atii$gene), ]))), 
    name = "Expr",
    col = colorRamp2(c(-1, 0, 1), c("purple", "black", "yellow")),
    show_row_names = F,
    cluster_rows = F,
    row_split = atii$Clt,
    column_split = sapply(colnames(temp_expr), function(x) { str_split(x, "\\s+")[[1]][2] }),
    border = T
)



### TCGA survival by cell composition
library(TCGAbiolinks)

luad_cli <- readRDS("TCGA/LUAD_clinical.rds")
rownames(luad_cli) <- luad_cli$submitter_id

lusc_cli <- readRDS("TCGA/LUSC_clinical.rds")
rownames(lusc_cli) <- lusc_cli$submitter_id

dec <- readRDS("09_bulk/RNA_seq/MuSic/tcga.rds")
est <- as.data.frame(dec$Est.prop.weighted)

est$ID <- sapply(rownames(est), function(x) {
    temp = str_split(x, "\\.")[[1]]
    paste(temp[1], temp[2], temp[3], sep = "-")
})


for (i in colnames(est)) {
    if (i == "ID") {
        next
    }
    
    for (j in c("LUAD", "LUSC")) {
        temp_est <- est[str_detect(rownames(est), j), ]
        temp_est$group = ifelse(temp_est[, i] > median(temp_est[, i]), "high", "low")
        
        if (j == "LUAD") {
            luad_cli[temp_est$ID, "Cell"] = temp_est$group
            clinical = luad_cli
        } else {
            lusc_cli[temp_est$ID, "Cell"] = temp_est$group
            clinical = lusc_cli
        }
        
        p <- TCGAanalyze_survival(
            clinical,
            "Cell",
            filename = NULL,
            risk.table = F,
            conf.int = F,
            ncensor.plot = TRUE,
            labels = c("high", "low"),
            color = c("red", "blue"),
            conf.int.style = "step",  # customize style of confidence intervals
        )
        p <- p$plot
        p <- p + 
            labs(title = paste0(i, " (", j, ")")) +
            ggplot2::theme(
                aspect.ratio = 1,
                axis.text = element_text(size = 15),
                axis.title = element_text(size = 15),
                legend.position = c(0.75, 0.9),
                legend.title = element_blank(),
                legend.text = element_text(size = 15),
                plot.title = element_text(size = 15, hjust = 0.5)
            )
        
        ggsave(
            filename = paste0("09_bulk/RNA_seq/MuSic/Surv/median/", i, "_", j, ".pdf"),
            plot = p,
            width = 4,
            height = 4
        )
        
        
        # 3/4
        temp_est$group = ifelse(temp_est[, i] > summary(temp_est[, i])[5], "high", "low")
        
        if (j == "LUAD") {
            luad_cli[temp_est$ID, "Cell"] = temp_est$group
            clinical = luad_cli
        } else {
            lusc_cli[temp_est$ID, "Cell"] = temp_est$group
            clinical = lusc_cli
        }
        
        p <- TCGAanalyze_survival(
            clinical,
            "Cell",
            filename = NULL,
            risk.table = F,
            conf.int = F,
            ncensor.plot = TRUE,
            labels = c("high", "low"),
            color = c("red", "blue"),
            conf.int.style = "step",  # customize style of confidence intervals
        )
        p <- p$plot
        p <- p + 
            labs(title = paste0(i, " (", j, ")")) +
            ggplot2::theme(
                aspect.ratio = 1,
                axis.text = element_text(size = 15),
                axis.title = element_text(size = 15),
                legend.position = c(0.75, 0.9),
                legend.title = element_blank(),
                legend.text = element_text(size = 15),
                plot.title = element_text(size = 15, hjust = 0.5)
            )
        
        ggsave(
            filename = paste0("09_bulk/RNA_seq/MuSic/Surv/75/", i, "_", j, ".pdf"),
            plot = p,
            width = 4,
            height = 4
        )
    }
}



### east

make_heatmap_east <- function(data, meta, scale_dot_size=0.9) {
    col_fun <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
    
    colors = list(
        Stage = c("I"="#65A9A3", "II"="#4A933E", "III"="#EC7A21", "IV"="#D73F47", "LUAD_Normal"="#FECC1B", "LUSC_Normal"="#778793", "BSPN"="#73BDFF"),
        Disease = c("LUAD" = "#0084D1", "NL(AD)"="#FECC1B", "NL(AD"="#778793", "BSPN"="#73BDFF", "NL(SC)"="#778793", "LUSC"="#A0C807"),
        Batch = c("1"="#E9B0B7", "2"="#90D0E0")
    )
    
    ra <- rowAnnotation(
        Disease = sapply(rownames(data), function(x) {
            temp = str_split(x, "_")[[1]]
            ifelse(temp[2] == "Normal", "NL(AD)", "LUAD")
        }),
        Stage = meta[rownames(data), "STAGE"],
        col = colors,
        show_annotation_name = F,
        show_legend = T
    )
    
    ba <- HeatmapAnnotation(
        Disease = sapply(colnames(data), function(x) {
            temp = str_split(x, "_")[[1]]
            ifelse(temp[2] == "Normal", "NL(AD)", "LUAD")
        }),
        Stage = meta[colnames(data), "STAGE"],
        col = colors,
        show_annotation_name = T,
        show_legend = F
    )
    
    Heatmap(
        data,
        name = "Corr",
        col = col_fun,
        right_annotation = ra,
        bottom_annotation = ba,
        show_row_names = F,
        show_column_names = F,
        # rect_gp = gpar(type = "none"),
        # cell_fun = function(j, i, x, y, width, height, fill) {
        #     s = min(unit.c(convertWidth(width, "cm"), convertHeight(height, "cm")))
        #     grid.circle(
        #         x = x,
        #         y = y,
        #         r = s * (abs(data[i, j]) + 0.1) * scale_dot_size,
        #         gp = gpar(col = "white", fill = fill)
        #     )
        # }
    )
}


dec <- readRDS("09_bulk/RNA_seq/MuSic/east.rds")
est <- dec$Est.prop.weighted


immu <- c(
    "B",
    "CD4",
    "CD8",
    "DC",
    "Mast",
    "Mo",
    "NK",
    "Tregs",
    "Gran"
)


temp <- melt(as.matrix(dec$Est.prop.weighted))
temp$Disease <- sapply(temp$Var1, function(x) {
    temp = str_split(x, "_")[[1]]
    ifelse(temp[2] == "Normal", "NL(AD)", "LUAD")
})
temp$Type <- ifelse(temp$Var2 %in% immu, "Immune", "Non-immune") 
temp$Var2 <- factor(temp$Var2, levels = sort(unique(as.character(temp$Var2))))

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
    scale_color_manual(values =c(
        "LUAD" = "#0084D1", 
        "NL(AD)"="#FECC1B", 
        "NL(AD"="#778793", 
        "BSPN"="#FB290F", 
        "NL(SC)"="#778793", 
        "LUSC"="#A0C807")
    )


ggsave(
    filename = "09_bulk/RNA_seq/MuSic/east_boxplot.pdf",
    plot = p,
    width = 12,
    height = 8
)



pdf("09_bulk/RNA_seq/MuSic/east_heatmap_all_sample.pdf", width = 5, height = 4)
make_heatmap_east(cor(t(est)), east_meta)
dev.off()

pdf("09_bulk/RNA_seq/MuSic/east_heatmap_all_sample_spearman.pdf", width = 5, height = 4)
make_heatmap_east(cor(t(est), method = "spearman"), east_meta)
dev.off()


pdf("09_bulk/RNA_seq/MuSic/east_heatmap_immune.pdf", width = 5, height = 4)
temp_est <- cor(t(est[, colnames(est) %in% immu]))
temp_est[is.na(temp_est)] <- 0
make_heatmap_east(temp_est, east_meta)
dev.off()


pdf("09_bulk/RNA_seq/MuSic/east_heatmap_immune_spearman.pdf", width = 5, height = 4)
temp_est <- cor(t(est[, colnames(est) %in% immu]), method = "spearman")
temp_est[is.na(temp_est)] <- 0
make_heatmap_east(temp_est, east_meta)
dev.off()


pdf("09_bulk/RNA_seq/MuSic/east_heatmap_non_immune.pdf", width = 5, height = 4)
make_heatmap_east(cor(t(est[, !colnames(est) %in% immu])), east_meta)
dev.off()


pdf("09_bulk/RNA_seq/MuSic/east_heatmap_non_immune_spearman.pdf", width = 5, height = 4)
make_heatmap_east(cor(t(est[, !colnames(est) %in% immu]), method = "spearman"), east_meta)
dev.off()


### Bulk && East

#### East
dec <- readRDS("09_bulk/RNA_seq/MuSic/east.rds")
est <- dec$Est.prop.weighted


#### GTEx
gtex <- readRDS("09_bulk/RNA_seq/MuSic/gtex.rds")
gtex_est <- gtex$Est.prop.weighted

#### bulk
bulk <- readRDS("09_bulk/RNA_seq/MuSic/bulk.rds")
bulk_est <- bulk$Est.prop.weighted
# bulk_est <- bulk_est[str_detect(rownames(bulk_est), "^LU"), ]


### TCGA
tcga <- readRDS("09_bulk/RNA_seq/MuSic/tcga.rds")
tcga_est <- tcga$Est.prop.weighted

### merge est
rnames <- intersect(sort(colnames(bulk_est)), sort(colnames(est)))
rnames <- intersect(colnames(gtex_est), rnames)
rnames <- intersect(colnames(tcga_est), rnames)

merged_est <- cbind(t(bulk_est)[rnames, ], t(est)[rnames, ])
merged_est <- cbind(merged_est[rnames, ], t(gtex_est)[rnames, ])
merged_est <- cbind(merged_est[rnames, ], t(tcga_est)[rnames, ])


### merge meta
meta <- bulk_meta[, "Stage", drop=F]
meta$Disease <- sapply(rownames(meta), function(x) { str_replace_all(x, "\\d+", "") })
meta$Source <- "In-house"

temp_meta <- east_meta[, "STAGE", drop = F]
temp_meta$Disease <- sapply(rownames(temp_meta), function(x) {
    temp = str_split(x, "_")[[1]]
    ifelse(temp[2] == "Normal", "NL(AD)", "LUAD")
})
colnames(temp_meta)[1] <- "Stage"
temp_meta$Source = "East Asians"

meta <- rbind(meta, temp_meta)

temp_meta <- data.frame(
    Stage = rep("Normal", nrow(gtex_est)),
    Disease = rep("NL", nrow(gtex_est)),
    Source = rep("GTEx", nrow(gtex_est)),
    row.names = rownames(gtex_est)
)
### read bulk
tcga_meta <- data.frame(
    Disease = sapply(rownames(tcga_est), function(x) { str_split(x, "\\.")[[1]][8] }),
    Stage = sapply(rownames(tcga_est), function(x) { str_split(x, "\\.")[[1]][9] }),
    row.names = rownames(tcga_est)
)

# tcga_meta <- tcga_meta[tcga_meta$Stage != "", ]
tcga_meta$Source = "TCGA"

meta <- rbind(meta, temp_meta)
meta <- rbind(meta, tcga_meta)
meta <- meta[colnames(merged_est), ]


make_heatmap_merged <- function(
    data, meta, 
    scale_dot_size=0.9, 
    row_split=NULL, column_split=NULL,
    cluster_rows = T, cluster_columns = T
    ) {
    col_fun <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
    
    colors = list(
        Stage = c("I"="#65A9A3", "II"="#4A933E", "III"="#EC7A21", "IV"="#D73F47", "Normal"="#73BDFF"),
        Disease = c(
            "LUAD" = "#0084D1", 
            "NL(AD)"="#FECC1B", 
            "NL"="#778793", 
            "BSPN"="#73BDFF", 
            "NL(SC)"="#778793", 
            "LUSC"="#A0C807"
        ),
        Source = c(
            "In-house"="#F8AFA8", 
            "East Asians"="#FDDDA0",
            "GTEx"="#73A089",
            "TCGA"="#877E45"
        )
    )
    
    ra <- rowAnnotation(
        Disease = meta[rownames(data), "Disease"],
        Stage = meta[rownames(data), "Stage"],
        Source = meta[rownames(data), "Source"],
        col = colors,
        show_annotation_name = F,
        show_legend = T
    )
    
    ba <- HeatmapAnnotation(
        Disease = meta[colnames(data), "Disease"],
        Stage = meta[colnames(data), "Stage"],
        Source = meta[colnames(data), "Source"],
        col = colors,
        show_annotation_name = T,
        show_legend = F
    )
    
    Heatmap(
        data,
        name = "Corr",
        col = col_fun,
        right_annotation = ra,
        bottom_annotation = ba,
        show_row_names = F,
        show_column_names = F,
        row_split = row_split,
        column_split = column_split,
        cluster_rows = cluster_rows,
        cluster_columns = cluster_columns,
        border = T
    )
}


pdf("09_bulk/RNA_seq/MuSic/merged_heatmap_all_sample.pdf", width = 5, height = 4)
make_heatmap_merged(cor(merged_est), meta)
dev.off()


pdf("09_bulk/RNA_seq/MuSic/merged_heatmap_immune.pdf", width = 5, height = 4)
temp_cor = cor(merged_est[rownames(merged_est) %in% immu, ])
temp_cor[is.na(temp_cor)] <- 0
make_heatmap_merged(temp_cor, meta)
dev.off()


pdf("09_bulk/RNA_seq/MuSic/merged_heatmap_non_immune.pdf", width = 5, height = 4)
temp_cor = cor(merged_est[!rownames(merged_est) %in% immu, ])
temp_cor[is.na(temp_cor)] <- 0
make_heatmap_merged(temp_cor, meta)
dev.off()


#### boxplot

temp_est <- melt(merged_est)
colnames(temp_est) <- c("CellType", "Sample", "value")
temp_est$Source <- as.character(meta[as.character(temp_est$Sample), "Source"])
temp_est$Disease <- as.character(meta[as.character(temp_est$Sample), "Disease"])
temp_est$Type <- ifelse(temp_est$CellType %in% immu, "Immune", "Non-immune") 

temp_est$CellType <- factor(temp_est$CellType, levels = sort(unique(as.character(temp_est$CellType))))

p <- ggplot(temp_est, aes(x=CellType, y = value, color = Disease)) +
    geom_boxplot() +
    facet_grid(Source~Type, scales = "free_x", space = "free") +
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
    scale_color_manual(values =c(
        "LUAD" = "#0084D1", 
        "NL(AD)"="#FECC1B", 
        "NL"="#778793", 
        "BSPN"="#F8AFA8", 
        "NL(SC)"="#778793", 
        "LUSC"="#A0C807"
    ))

ggsave(
    filename = "09_bulk/RNA_seq/MuSic/merged_boxplot_with_BSPN.pdf",
    plot = p,
    width = 12,
    height = 8
)


#### PCA

res.pca <- prcomp(t(merged_est))

p <- ggplot(
    data.frame(x = 1:ncol(res.pca$x), y = apply(res.pca$x, 2, sd)), 
    aes(x = x, y = y)
) + geom_point()

print(p)


p1 <- autoplot(res.pca, data = meta[rownames(res.pca$x), ],  colour = 'Disease', frame = TRUE) + 
    my_theme + theme(aspect.ratio=1, legend.position = "bottom") +
    labs(title = "Disease") +
    scale_color_manual(values = c(
        "LUAD" = "#0084D1", 
        "NL(AD)"="#FECC1B", 
        "NL"="#778793", 
        "BSPN"="#73BDFF", 
        "NL(SC)"="#778793", 
        "LUSC"="#A0C807"
    )) +
    guides(
        color = guide_legend(override.aes = list(size = 4))
    )

p2 <- autoplot(res.pca, data = meta[rownames(res.pca$x), ],  colour = 'Stage', frame = TRUE) + 
    my_theme + theme(aspect.ratio=1, legend.position = "bottom")  +
    labs(title = "Stage") +
    scale_color_manual(values = c(
        "I"="#65A9A3", "II"="#4A933E", 
        "III"="#EC7A21", "IV"="#D73F47", 
        "Normal"="#FECC1B", 
        "BSPN"="#73BDFF"
    )) +
    guides(
        color = guide_legend(override.aes = list(size = 4))
    )


p3 <- autoplot(res.pca, data = meta[rownames(res.pca$x), ],  colour = 'Source', frame = TRUE) + 
    my_theme + theme(aspect.ratio=1, legend.position = "bottom")  +
    labs(title = "Stage") +
    scale_color_manual(values = c(
        "In-house"="#F8AFA8", 
        "East Asians"="#FDDDA0",
        "GTEx"="#73A089"
    )) +
    guides(
        color = guide_legend(override.aes = list(size = 4))
    )

ggsave(
    filename = "09_bulk/RNA_seq/MuSic/merged_pca_gtex.pdf",
    plot = cowplot::plot_grid(p1, p2, p3, ncol = 3),
    device = "pdf",
    width = 18,
    height = 6
)




n.pcs = 10
dot.size = 0.6
res.umap <- umap(res.pca$x[, 1:n.pcs])

data = as.data.frame(res.umap$layout[, 1:2])
colnames(data) <- c("UMAP1", "UMAP2")

data$Disease <- meta[rownames(data), "Disease"]
data$Source <- meta[rownames(data), "Source"]
data$Stage <- meta[rownames(data), "Stage"]


p1 <- ggplot(data, aes(x = UMAP1, y = UMAP2, color  = Disease)) +
    geom_point(size = dot.size) +
    theme_bw() +
    my_theme + 
    theme(legend.position = "bottom") +
    labs(color = "") +
    scale_color_manual(values = c(
        "LUAD" = "#0084D1", 
        "NL(AD)"="#FECC1B", 
        "NL"="#778793", 
        "BSPN"="#73BDFF", 
        "NL(SC)"="#778793", 
        "LUSC"="#A0C807"
    )) +
    guides(
        color = guide_legend(override.aes = list(size = 4), ncol = 3)
    )

p2 <- ggplot(data, aes(x = UMAP1, y = UMAP2, color  = Stage)) +
    geom_point(size = dot.size) +
    theme_bw() +
    my_theme + 
    theme(legend.position = "bottom") +
    labs(color = "") +
    scale_color_manual(values = c(
        "I"="#65A9A3", "II"="#4A933E", 
        "III"="#EC7A21", "IV"="#D73F47", 
        "Normal"="#FECC1B", 
        "BSPN"="#73BDFF"
    )) +
    guides(
        color = guide_legend(override.aes = list(size = 4))
    )


p3 <- ggplot(data, aes(x = UMAP1, y = UMAP2, color  = Source)) +
    geom_point(size = dot.size) +
    theme_bw() +
    my_theme + 
    theme(legend.position = "bottom") +
    labs(color = "") +
    scale_color_manual(values = c(
        "In-house"="#F8AFA8", 
        "East Asians"="#FDDDA0",
        "GTEx"="#73A089"
    )) +
    guides(
        color = guide_legend(override.aes = list(size = 4), nrow=2)
    )


ggsave(
    filename = "09_bulk/RNA_seq/MuSic/merged_umap_gtex.pdf",
    plot = cowplot::plot_grid(p1, p2, p3, ncol = 3),
    device = "pdf",
    width = 18,
    height = 6
)



### using cutree and max number of samples as label
make_heatmap_number_merged <- function(
    data, meta, 
    scale_dot_size=0.9, 
    row_split=NULL, column_split=NULL,
    cluster_rows = T, cluster_columns = T,
    show_row_names = F, show_column_names = T
) {
    col_fun <- colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
    
    colors = list(
        Stage = c("I"="#65A9A3", "II"="#4A933E", "III"="#EC7A21", "IV"="#D73F47", "Normal"="#73BDFF"),
        Disease = c(
            "LUAD" = "#0084D1", 
            "NL(AD)"="#FECC1B", 
            "NL"="#778793", 
            "BSPN"="#73BDFF", 
            "NL(SC)"="#778793", 
            "LUSC"="#A0C807"
        ),
        Source = c(
            "In-house"="#F8AFA8", 
            "East Asians"="#FDDDA0",
            "GTEx"="#73A089",
            "TCGA"="#877E45"
        )
    )
    
    ra <- rowAnnotation(
        Disease = meta[rownames(data), "Disease"],
        Stage = meta[rownames(data), "Stage"],
        Source = meta[rownames(data), "Source"],
        col = colors,
        show_annotation_name = F,
        show_legend = T
    )
    
    # ba <- HeatmapAnnotation(
    #     Disease = meta[colnames(data), "Disease"],
    #     Stage = meta[colnames(data), "Stage"],
    #     Source = meta[colnames(data), "Source"],
    #     col = colors,
    #     show_annotation_name = T,
    #     show_legend = F
    # )
    
    Heatmap(
        data,
        name = "Weight",
        col = col_fun,
        right_annotation = ra,
        # bottom_annotation = ba,
        show_row_names = show_row_names,
        show_column_names = show_column_names,
        row_split = row_split,
        column_split = column_split,
        cluster_rows = cluster_rows,
        cluster_columns = cluster_columns,
        border = T
    )
}


### GTEx、EA、inhouse merged
temp_est <- t(merged_est)

set.seed(42)
temp_kmeans <- kmeans(temp_est, 9)
temp_tree <- as.data.frame(temp_kmeans$cluster)
colnames(temp_tree) <- "Clt"

temp_tree$Disease <- meta[rownames(temp_tree), "Disease"]
temp_tree <- temp_tree[order(temp_tree$Clt), , drop = F]

temp_tree$Disease[temp_tree$Disease == "NL"] <- "NL(AD)"

temp_score <- temp_tree %>%
    group_by(Clt, Disease) %>%
    add_tally() %>%
    unique() %>%
    group_by(Clt) %>%
    mutate(p = n / sum(n)) %>%
    group_by(Clt) %>%
    top_n(1, wt = p) %>%
    as.data.frame()

mean(temp_score$p)

dir.create("09_bulk/RNA_seq/MuSic/Est", showWarnings = F)

pdf("09_bulk/RNA_seq/MuSic/Est/merged_heatmap.pdf", width = 5, height = 4)
make_heatmap_number_merged(
    temp_est,
    meta, 
    row_split = temp_tree[rownames(temp_est), "Clt"], 
    # column_split = rownames(temp_tree),
    cluster_rows = F,
    cluster_columns = T,
    show_c
)
dev.off()



#### inhouse
temp_est <- t(merged_est[, str_detect(colnames(merged_est), "^LU")])

set.seed(42)
temp_kmeans <- kmeans(temp_est, 4)
temp_tree <- as.data.frame(temp_kmeans$cluster)
colnames(temp_tree) <- "Clt"

temp_tree$Disease <- meta[rownames(temp_tree), "Disease"]
temp_tree <- temp_tree[order(temp_tree$Clt), , drop = F]


temp_score <- temp_tree %>%
    group_by(Clt, Disease) %>%
    add_tally() %>%
    unique() %>%
    group_by(Clt) %>%
    mutate(p = n / sum(n)) %>%
    group_by(Clt) %>%
    top_n(1, wt = p) %>%
    as.data.frame()

mean(temp_score$p)


pdf("09_bulk/RNA_seq/MuSic/Est/inhouse.pdf", width = 5, height = 4)
make_heatmap_number_merged(
    temp_est,
    meta, 
    row_split = temp_tree[rownames(temp_est), "Clt"], 
    # column_split = rownames(temp_tree),
    cluster_rows = F,
    cluster_columns = T,
    show_c
)
dev.off()


##### EA
temp_est <- t(merged_est[, str_detect(colnames(merged_est), "^A\\d+")])

set.seed(42)
temp_kmeans <- kmeans(temp_est, 5)
temp_tree <- as.data.frame(temp_kmeans$cluster)
colnames(temp_tree) <- "Clt"

temp_tree$Disease <- meta[rownames(temp_tree), "Disease"]
temp_tree <- temp_tree[order(temp_tree$Clt), , drop = F]


temp_score <- temp_tree %>%
    group_by(Clt, Disease) %>%
    add_tally() %>%
    unique() %>%
    group_by(Clt) %>%
    mutate(p = n / sum(n)) %>%
    group_by(Clt) %>%
    top_n(1, wt = p) %>%
    as.data.frame()

mean(temp_score$p)


pdf("09_bulk/RNA_seq/MuSic/Est/East_Asians.pdf", width = 5, height = 4)
make_heatmap_number_merged(
    temp_est,
    meta, 
    row_split = temp_tree[rownames(temp_est), "Clt"], 
    # column_split = rownames(temp_tree),
    cluster_rows = F,
    cluster_columns = T,
    show_c
)
dev.off()



##### GTEx
temp_est <- t(merged_est[, str_detect(colnames(merged_est), "GTE")])

set.seed(42)
temp_kmeans <- kmeans(temp_est, 5)
temp_tree <- as.data.frame(temp_kmeans$cluster)
colnames(temp_tree) <- "Clt"

temp_tree$Disease <- meta[rownames(temp_tree), "Disease"]
temp_tree <- temp_tree[order(temp_tree$Clt), , drop = F]

temp_score <- temp_tree %>%
    group_by(Clt, Disease) %>%
    add_tally() %>%
    unique() %>%
    group_by(Clt) %>%
    mutate(p = n / sum(n)) %>%
    group_by(Clt) %>%
    top_n(1, wt = p) %>%
    as.data.frame()

mean(temp_score$p)


pdf("09_bulk/RNA_seq/MuSic/Est/GTEx.pdf", width = 5, height = 4)
make_heatmap_number_merged(
    temp_est,
    meta, 
    row_split = temp_tree[rownames(temp_est), "Clt"], 
    # column_split = rownames(temp_tree),
    cluster_rows = F,
    cluster_columns = T,
    show_c
)
dev.off()



#### Inhouse LUSC
temp_est <- t(merged_est[, str_detect(colnames(merged_est), "LUSC")])

fviz_nbclust(temp_est, kmeans, method = "wss") +
    geom_vline(xintercept = 4, linetype = 2)+
    labs(subtitle = "Elbow method")

set.seed(42)
temp_kmeans <- kmeans(temp_est, 4)
temp_tree <- as.data.frame(temp_kmeans$cluster)
colnames(temp_tree) <- "Clt"

temp_tree$Disease <- meta[rownames(temp_tree), "Disease"]
temp_tree <- temp_tree[order(temp_tree$Clt), , drop = F]


pdf("09_bulk/RNA_seq/MuSic/Est/inhouse_LUSC.pdf", width = 5, height = 4)
make_heatmap_number_merged(
    temp_est,
    meta, 
    row_split = temp_tree[rownames(temp_est), "Clt"], 
    # column_split = rownames(temp_tree),
    cluster_rows = F,
    cluster_columns = T,
    show_column_names = T
)
dev.off()



##### Inhouse LUAD
temp_est <- t(merged_est[, str_detect(colnames(merged_est), "LUAD")])

fviz_nbclust(temp_est, kmeans, method = "wss") +
    geom_vline(xintercept = 4, linetype = 2)+
    labs(subtitle = "Elbow method")

set.seed(42)
temp_kmeans <- kmeans(temp_est, 4)
temp_tree <- as.data.frame(temp_kmeans$cluster)
colnames(temp_tree) <- "Clt"

temp_tree$Disease <- meta[rownames(temp_tree), "Disease"]
temp_tree <- temp_tree[order(temp_tree$Clt), , drop = F]


pdf("09_bulk/RNA_seq/MuSic/Est/inhouse_LUAD.pdf", width = 5, height = 4)
make_heatmap_number_merged(
    temp_est,
    meta, 
    row_split = temp_tree[rownames(temp_est), "Clt"], 
    # column_split = rownames(temp_tree),
    cluster_rows = F,
    cluster_columns = T,
    show_column_names = T
)
dev.off()



##### EA LUAD
temp_est <- t(merged_est[, str_detect(colnames(merged_est), "^A\\d+_Tumor")])

fviz_nbclust(temp_est, kmeans, method = "wss") +
    geom_vline(xintercept = 4, linetype = 2)+
    labs(subtitle = "Elbow method")

set.seed(42)
temp_kmeans <- kmeans(temp_est, 4)
temp_tree <- as.data.frame(temp_kmeans$cluster)
colnames(temp_tree) <- "Clt"

temp_tree$Disease <- meta[rownames(temp_tree), "Disease"]
temp_tree <- temp_tree[order(temp_tree$Clt), , drop = F]


pdf("09_bulk/RNA_seq/MuSic/Est/East_asians_LUAD.pdf", width = 5, height = 4)
make_heatmap_number_merged(
    temp_est,
    meta, 
    row_split = temp_tree[rownames(temp_est), "Clt"], 
    # column_split = rownames(temp_tree),
    cluster_rows = F,
    cluster_columns = T,
    show_column_names = T
)
dev.off()


##### Survival LUAD
library(TCGAbiolinks)

luad_cli <- readRDS("TCGA/LUAD_clinical.rds")
rownames(luad_cli) <- luad_cli$submitter_id

lusc_cli <- readRDS("TCGA/LUSC_clinical.rds")
rownames(lusc_cli) <- lusc_cli$submitter_id


temp_est <- as.data.frame(tcga_est[str_detect(rownames(tcga_est), "LUAD"), ])


if (file.exists("09_bulk/RNA_seq/MuSic/Est/TCGA_LUAD_tree.rds")) {
    temp_tree <- readRDS("09_bulk/RNA_seq/MuSic/Est/TCGA_LUAD_tree.rds")
} else {
    set.seed(42)
    temp_kmeans <- kmeans(temp_est, 4)
    temp_tree <- as.data.frame(temp_kmeans$cluster)
    colnames(temp_tree) <- "Clt"
    saveRDS(temp_tree, "09_bulk/RNA_seq/MuSic/Est/TCGA_LUAD_tree.rds")
}


pdf("09_bulk/RNA_seq/MuSic/Est/TCGA_LUAD.pdf", width = 5, height = 4)
make_heatmap_number_merged(
    temp_est,
    meta, 
    row_split = temp_tree[rownames(temp_est), "Clt"], 
    # column_split = rownames(temp_tree),
    cluster_rows = F,
    cluster_columns = T,
    show_column_names = T
)
dev.off()



temp_est$ID <- sapply(as.character(rownames(temp_est)), function(x) {
    temp = str_split(x, "\\.")[[1]]
    paste(temp[1], temp[2], temp[3], sep = "-")
})
temp_est$Clt <- temp_tree[rownames(temp_est), "Clt"]
luad_cli[temp_est$ID, "group"] <- temp_est$Clt


p <- TCGAanalyze_survival(
    luad_cli,
    "group",
    filename = NULL,
    risk.table = F,
    conf.int = F,
    ncensor.plot = TRUE,
    labels = c("1", "2", "3", "4"),
    # color = c("red", "blue"),
    conf.int.style = "step",  # customize style of confidence intervals
)

p <- p$plot
p <- p + 
    labs(title = "LUAD") +
    ggplot2::theme(
        aspect.ratio = 1,
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.position = c(0.75, 0.9),
        legend.title = element_blank(),
        legend.text = element_text(size = 15),
        plot.title = element_text(size = 15, hjust = 0.5)
    )

p

ggsave(
    plot = p,
    filename = "09_bulk/RNA_seq/MuSic/Est/TCGA_LUAD_surv.pdf",
    width = 4,
    height = 4
)


###### LUSC
temp_est <- as.data.frame(tcga_est[str_detect(rownames(tcga_est), "LUSC"), ])


if (file.exists("09_bulk/RNA_seq/MuSic/Est/TCGA_LUSC_tree.rds")) {
    temp_tree <- readRDS("09_bulk/RNA_seq/MuSic/Est/TCGA_LUSC_tree.rds")
} else {
    set.seed(42)
    temp_kmeans <- kmeans(temp_est, 4)
    temp_tree <- as.data.frame(temp_kmeans$cluster)
    colnames(temp_tree) <- "Clt"
    saveRDS(temp_tree, "09_bulk/RNA_seq/MuSic/Est/TCGA_LUSC_tree.rds")
}


pdf("09_bulk/RNA_seq/MuSic/Est/TCGA_LUSC.pdf", width = 5, height = 4)
make_heatmap_number_merged(
    temp_est,
    meta, 
    row_split = 4, #temp_tree[rownames(temp_est), "Clt"], 
    # column_split = rownames(temp_tree),
    cluster_rows = F,
    cluster_columns = T,
    show_column_names = T
)
dev.off()

temp_est$ID <- sapply(as.character(rownames(temp_est)), function(x) {
    temp = str_split(x, "\\.")[[1]]
    paste(temp[1], temp[2], temp[3], sep = "-")
})
temp_est$Clt <- temp_tree[rownames(temp_est), "Clt"]
lusc_cli[temp_est$ID, "group"] <- temp_est$Clt


p <- TCGAanalyze_survival(
    lusc_cli,
    "group",
    filename = NULL,
    risk.table = F,
    conf.int = F,
    ncensor.plot = TRUE,
    labels = c("1", "2", "3", "4"),
    # color = c("red", "blue"),
    conf.int.style = "step",  # customize style of confidence intervals
)

p <- p$plot
p <- p + 
    labs(title = "LUAD") +
    ggplot2::theme(
        aspect.ratio = 1,
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.position = c(0.75, 0.9),
        legend.title = element_blank(),
        legend.text = element_text(size = 15),
        plot.title = element_text(size = 15, hjust = 0.5)
    )

p

ggsave(
    plot = p,
    filename = "09_bulk/RNA_seq/MuSic/Est/TCGA_LUSC_surv.pdf",
    width = 4,
    height = 4
)


#### AUC and ROC test

temp_est <- t(as.data.frame(merged_est[, rownames(meta)[meta$Disease %in% c("LUAD", "NL", "NL(AD)")]]))

fviz_nbclust(temp_est, kmeans, method = "wss") +
    geom_vline(xintercept = 4, linetype = 2)+
    labs(subtitle = "Elbow method")

set.seed(42)
temp_kmeans <- kmeans(temp_est, 4)
temp_tree <- as.data.frame(temp_kmeans$cluster)
colnames(temp_tree) <- "Clt"


pdf("09_bulk/RNA_seq/MuSic/Est/LUAD_NL.pdf", width = 5, height = 4)
make_heatmap_number_merged(
    temp_est,
    meta, 
    row_split = temp_tree[rownames(temp_est), "Clt"], 
    # column_split = rownames(temp_tree),
    cluster_rows = F,
    cluster_columns = T,
    show_column_names = T
)
dev.off()


temp_est <- t(as.data.frame(merged_est[, rownames(meta)[meta$Disease %in% c("LUSC", "NL")]]))

fviz_nbclust(temp_est, kmeans, method = "wss") +
    geom_vline(xintercept = 4, linetype = 2)+
    labs(subtitle = "Elbow method")

set.seed(42)
temp_kmeans <- kmeans(temp_est, 4)
temp_tree <- as.data.frame(temp_kmeans$cluster)
colnames(temp_tree) <- "Clt"


pdf("09_bulk/RNA_seq/MuSic/Est/LUSC_NL.pdf", width = 5, height = 4)
make_heatmap_number_merged(
    temp_est,
    meta, 
    row_split = temp_tree[rownames(temp_est), "Clt"], 
    # column_split = rownames(temp_tree),
    cluster_rows = F,
    cluster_columns = T,
    show_column_names = T
)
dev.off()



temp_est <- t(as.data.frame(merged_est[, rownames(meta)[meta$Disease %in% c("LUSC", "LUAD")]]))

fviz_nbclust(temp_est, kmeans, method = "wss") +
    geom_vline(xintercept = 4, linetype = 2)+
    labs(subtitle = "Elbow method")

set.seed(42)
temp_kmeans <- kmeans(temp_est, 4)
temp_tree <- as.data.frame(temp_kmeans$cluster)
colnames(temp_tree) <- "Clt"


pdf("09_bulk/RNA_seq/MuSic/Est/LUSC_LUAD.pdf", width = 5, height = 4)
make_heatmap_number_merged(
    temp_est,
    meta, 
    row_split = temp_tree[rownames(temp_est), "Clt"], 
    # column_split = rownames(temp_tree),
    cluster_rows = F,
    cluster_columns = T,
    show_column_names = T
)
dev.off()


write.csv(merged_est, "09_bulk/RNA_seq/MuSic/All_est.csv")
write.csv(meta, "09_bulk/RNA_seq/MuSic/All_est_meta.csv")