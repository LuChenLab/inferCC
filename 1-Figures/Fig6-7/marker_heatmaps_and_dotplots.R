options(stringsAsFactors = F)
library(Seurat)
library(ComplexHeatmap)
library(circlize)
library(wesanderson)
library(openxlsx)
library(stringr)
library(reshape2)
library(ggthemes)
library(dplyr)
library(tidyr)

extrafont::loadfonts()

root.dir = "LungCancer10x"


here <- function(path) {
    paste(root.dir, path, sep = "/")
}


obj <- readRDS(here("02_rds/seurat_obj_multi_reduction.rds"))

## Format meta info
# nm_meta = read.xlsx("08_pseudobulk/NM_meta.xlsx")
# 
# nm_patients = list()
# for(i in 1:nrow(nm_meta)) {
#     nm_patients[[nm_meta[i, "SampleID"]]] = paste0(
#         "DL",
#         nm_meta[i, "PatientID"]
#     )
# }
# 
# obj@meta.data$PatientID[obj@meta.data$Batch == 3] = NA
# obj@meta.data$PatientID = as.character(obj@meta.data$PatientID)
# obj@meta.data$PatientID[is.na(obj@meta.data$PatientID)] = sapply(as.character(obj@meta.data$SampleID)[is.na(obj@meta.data$PatientID)], function(x) {
#     nm_patients[[x]]
# })
# 
# rownames(nm_meta) <- nm_meta$SampleID
# 
# obj@meta.data$Disease <- as.character(obj@meta.data$Disease)
# obj@meta.data$Disease[obj@meta.data$Batch == 3] <- nm_meta[
#     as.character(obj@meta.data$SampleID[obj@meta.data$Batch == 3]), 
#     "Disease"
# ]
# 
# 
short_names <- c(
    "APII"="AT2",
    "Alveolar II"="AT2",
    "AT1"="AT1",
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
# 
# 
# obj@meta.data$cell_name = as.character(obj@meta.data$cell_name)
# 
# obj@meta.data$cell_name1 <- short_names[obj@meta.data$cell_name, 1]
# 
# unique(obj@meta.data[, c("cell_name", "cell_name1")])
# 
# 
# markers = read.xlsx(here("20190701_gene_markers.xlsx"), sheet = 3)
# markers <- na.omit(markers[, 1:2])
# markers <- unique(markers)
#  
# 
# new_obj <- CreateSeuratObject(
#     obj@raw.data[as.character(c(markers$Markers, "ITGAM", "CSF1R")), ],
#     meta.data = obj@meta.data
# )
# 
# new_obj@scale.data <- obj@scale.data[as.character(c(markers$Markers, "CSF1R", "ITGAM")), ]
# 
# saveRDS(new_obj, here("02_rds/seurat_obj_markers.rds"))
# 
# ### Read cell markers
# new_obj <- readRDS(here("02_rds/seurat_obj_markers.rds"))


expr <- obj@scale.data
# meta <- read.csv(here("02_rds/meta_after_singleR.csv"), row.names = 1, stringsAsFactors = F)
# # meta$cell_short = meta$cell_name1
# meta <- meta[meta$cell_short != "" & !is.na(meta$cell_short), ]
# meta$cell_short[meta$cell_short == "Mφ"] = "Mo"

meta <- readRDS("11_CNV/meta.rds")
meta$cell_short = as.character(meta$cell_short)
meta$cell_short[meta$cell_short == "Mo"] = "Mφ"
meta$cell_short[meta$cell_short == "ATII"] = "AT2"
meta$cell_short[meta$cell_short == "Epi"] = "AT1"
# extract common cell types
cell_types = unique(meta$cell_short)
cell_types = cell_types[!is.na(cell_types)]

# order markers
markers = read.xlsx(here("20190701_gene_markers.xlsx"), sheet = 3)
markers <- na.omit(markers[, 1:2])
markers <- unique(markers)
temp_markers = markers[order(markers$Cells, markers$Markers), ]
temp_markers$Cells = short_names[temp_markers$Cells, 1]
temp_markers <- na.omit(temp_markers)

temp_markers <- temp_markers[temp_markers$Cells %in% c(cell_types, "Mo", "ATII"), ]

temp_markers$Cells[temp_markers$Cells == "Mo"] = "Mφ"
temp_markers$Cells[temp_markers$Cells == "ATII"] = "AT2"

non_immu = unique(temp_markers$Cells)[!unique(temp_markers$Cells) %in% immu]
cell_levels = c(sort(non_immu), sort(immu))

# temp_markers <- rbind(
#     temp_markers,
#     data.frame(
#         Markers = c("CSF1R", "ITGAM"), # "KLRB1", , "CD80"
#         Cells = c("Mφ", "Mφ")  #
#     )
# )


# order cells
set.seed(1)
meta$Cells = rownames(meta)

cells <- meta %>%
    group_by(cell_short) %>%
    sample_n(1000, replace = T) %>%
    dplyr::select(Cells, cell_short) %>%
    as.data.frame()


## prepare marker expr
temp <- melt(as.matrix(expr[,  which(colnames(expr) %in% cells$Cells)]))
temp$Disease <- as.character(meta[as.character(temp$Var2), "Disease"])
temp$Cell <- as.character(meta[as.character(temp$Var2), "cell_short"])
temp <- temp[temp$Var1 %in% markers$Markers, ]

temp <- temp %>%
    group_by(Var1, Disease, Cell) %>%
    mutate(value = mean(value)) %>%
    dplyr::select(Var1, Disease, Cell, value) %>%
    unique() %>%
    as.data.frame()

colnames(temp)[1] <- "gene"



immu = c(
    "B",
    "CD4",
    "CD8",
    "DC",
    #"EC",
    "Gran",
    "Mast",
    "Mo",
    "NK",
    "Tregs",
    "Mφ"
)


### sort cells
cells$isImmu <- cells$cell_short %in% immu
cells <- cells[order(cells$cell_short, cells$isImmu), ]

cells <- cells$Cells
cells <- cells[!is.na(cells)]


### sort markers
temp_markers$isImmu <- temp_markers$Cells %in% immu
temp_markers$Cells <- as.character(temp_markers$Cells)
temp_markers <- temp_markers[order(temp_markers$isImmu, temp_markers$Cells), ]
unique(temp_markers$Cells)

# cell_color = wes_palette("Moonrise3", length(cell_types), type = "continuous")
cell_color = c(
    wes_palette("Darjeeling1"), 
    wes_palette("Moonrise3"),
    wes_palette("Darjeeling2"),
    wes_palette("Moonrise1")
)[1:length(cell_types)]

names(cell_color) = cell_types

col_fun = colorRamp2(c(-3, 0, 3), c("blue", "white", "red"))


# left anno
ha = rowAnnotation(
    cell_left = anno_simple(meta[cells, "cell_short"]),
    col = list(
        cell_left=cell_color,
        Type=c(
            "Immune"="#E9B0B7",
            "Stomal"="#90D0E0"
        )
    ),
    show_legend = F,
    show_annotation_name = F
)


# top anno
la = HeatmapAnnotation(
    Type = ifelse(temp_markers$Cells %in% immu, "Immune", "Stromal"),
    col = list(
        cell=cell_color,
        Type=c(
            "Immune"="#E9B0B7",
            "Stromal"="#90D0E0"
        )
    ),
    show_legend = T,
    show_annotation_name = F,
    annotation_legend_param = list(
        labels_gp = gpar(fontsize = 15, fontfamily="Arial Unicode MS"),
        title_gp = gpar(fontsize = 20, fontfamily="Arial Unicode MS"),
        # grid_height = unit(2, "cm")
        direction = "horizontal"
    )
)


meta$cell_short = factor(meta$cell_short, levels = cell_levels)
temp_markers$Cells = factor(temp_markers$Cells, levels = cell_levels)


extrafont::loadfonts()

cairo_pdf(here("01_first_plots/marker_heatmap1.pdf"), height = 8, width = 16)
draw(
    Heatmap(
        t(as.matrix(expr[as.character(temp_markers$Markers), as.character(cells)])),
        top_annotation = la,
        # left_annotation = ha,
        cluster_rows = F,
        cluster_columns = F,
        name = "mat",
        show_row_names = F,
        show_column_names = T,
        show_heatmap_legend = T,
        column_split = temp_markers$Cells,
        row_split = meta[cells, "cell_short"],
        border = T,
        column_title_gp = gpar(fontsize=20, fontfamily="Arial Unicode MS"),
        column_title_rot = 0,
        row_title_gp = gpar(fontsize = 20, fontfamily="Arial Unicode MS"),
        row_title_rot = 0,
        heatmap_legend_param = list(
            title = "Expr",
            # legend_height = unit(10, "cm"),
            # legend_width = unit(10, "cm"),
            labels_gp = gpar(fontsize = 15, fontfamily="Arial Unicode MS"),
            title_gp = gpar(fontsize = 20, fontfamily="Arial Unicode MS"),
            direction = "horizontal"
        ),
        col = col_fun,
        column_names_gp = gpar(fontfamily="Arial Unicode MS"),
        row_names_gp = gpar(fontfamily="Arial Unicode MS"),
        use_raster = T
    ),
    merge_legend = TRUE,
    heatmap_legend_side = "bottom"
)
dev.off()



# bottom anno
ba = HeatmapAnnotation(
    cell = as.character(meta[cells, "cell_short"]),
    # cel = anno_mark(
    #     at = temp$m,
    #     labels = temp$cell_name1,
    #     which = "column",
    #     side="bottom",
    #     labels_gp = gpar(fontsize = 25)
    # ),
    col = list(
        cell=cell_color
    ),
    show_legend = F,
    show_annotation_name = F
)


# left anno
la = rowAnnotation(
    # cell = temp_markers$Cells,
    Type = ifelse(temp_markers$Cells %in% immu, "Immune", "Stromal"),
    col = list(
        cell=cell_color,
        Type=c(
            "Immune"="#E9B0B7",
            "Stromal"="#90D0E0"
        )
    ),
    show_legend = T,
    show_annotation_name = F,
    annotation_legend_param = list(
        labels_gp = gpar(fontsize = 15, fontfamily="Arial Unicode MS"),
        title_gp = gpar(fontsize = 20, fontfamily="Arial Unicode MS"),
        grid_height = unit(2, "cm")
    )
)


cairo_pdf(here("01_first_plots/marker_heatmap_vertical.pdf"), height = 16, width = 8)
draw(
    Heatmap(
        as.matrix(expr[temp_markers$Markers, cells]),
        left_annotation = la,
        # top_annotation = ba,
        cluster_rows = F,
        cluster_columns = F,
        name = "Expr",
        show_row_names = T,
        show_column_names = F,
        show_heatmap_legend = T,
        column_split = meta[cells, "cell_short"],
        row_split = temp_markers$Cells,
        border = T,
        column_title_gp = gpar(fontsize=20, fontfamily="Arial Unicode MS"),
        column_title_rot = 90,
        column_title_side = "top",
        row_title_gp = gpar(fontsize = 20, fontfamily="Arial Unicode MS"),
        row_title_rot = 0,
        col = col_fun,
        heatmap_legend_param = list(
            title = "Expr",
            # direction = "vertical",
            legend_height = unit(10, "cm"),
            legend_width = unit(10, "cm"),
            labels_gp = gpar(fontsize = 25, fontfamily="Arial Unicode MS"),
            title_gp = gpar(fontsize = 20, fontfamily="Arial Unicode MS")
            # title_position = "leftcenter"
        ),
        column_names_gp = gpar(fontfamily="Arial Unicode MS"),
        row_names_gp = gpar(fontfamily="Arial Unicode MS"),
        use_raster = T
    ),
    heatmap_legend_side = "right"
)

dev.off()




######
# DotPlot
######
obj <- readRDS(here("02_rds/seurat_obj_markers.rds"))

## Convert cell name to short
immu = c(
    "B",
    "CD4",
    "CD8",
    "DC",
    # "EC",
    "Gran",
    "Mast",
    "Mo",
    "NK",
    "Tregs",
    "Mφ"
)



markers = read.xlsx(here("20190701_gene_markers.xlsx"), sheet = 3)
markers <- na.omit(markers[, 1:2])
markers <- unique(markers)
temp_markers = markers[order(markers$Cells, markers$Markers), ]
temp_markers$Cells = short_names[temp_markers$Cells, 1]
temp_markers <- na.omit(temp_markers)

temp_markers <- temp_markers[temp_markers$Cells %in% c(cell_types, "Mo", "ATII"), ]

temp_markers$Cells[temp_markers$Cells == "Mo"] = "Mφ"
temp_markers$Cells[temp_markers$Cells == "ATII"] = "AT2"


markers <- temp_markers[temp_markers$Cells %in% immu, ]


PercentAbove <- function(x, threshold){
    return(length(x = x[x > threshold]) / length(x = x))
}


ModifiedDotData <- function(
    object,
    group.by,
    markers,
    cols.use = c("lightgrey", "blue"),
    col.min = -2.5,
    col.max = 2.5,
    dot.min = 0,
    dot.scale = 6,
    scale.by = 'radius'
) {
    scale.func <- switch(
        EXPR = scale.by,
        'size' = scale_size,
        'radius' = scale_radius,
        stop("'scale.by' must be either 'size' or 'radius'")
    )
    if (!missing(x = group.by)) {
        object <- SetAllIdent(object = object, id = group.by)
    }
    
    genes.plot <- intersect(row.names(object@raw.data), markers$Markers)
    
    data.to.plot <- data.frame(FetchData(object = object, vars.all = genes.plot))
    colnames(x = data.to.plot) <- genes.plot
    data.to.plot$cell <- rownames(x = data.to.plot)
    data.to.plot$id <- object@ident
    
    data.to.plot %>% gather(
        key = genes.plot,
        value = expression,
        -c(cell, id)
    ) -> data.to.plot
    data.to.plot %>%
        group_by(id, genes.plot) %>%
        summarize(
            avg.exp = mean(expm1(x = expression)),
            pct.exp = PercentAbove(x = expression, threshold = 0)
        ) -> data.to.plot
    data.to.plot %>%
        ungroup() %>%
        group_by(genes.plot) %>%
        mutate(avg.exp.scale = scale(x = avg.exp)) %>%
        mutate(avg.exp.scale = MinMax(
            data = avg.exp.scale,
            max = col.max,
            min = col.min
        )) ->  data.to.plot
    data.to.plot$genes.plot <- factor(
        x = data.to.plot$genes.plot,
        levels = rev(x = genes.plot)
    )
    # data.to.plot$genes.plot <- factor(
    #   x = data.to.plot$genes.plot,
    #   levels = rev(x = sub(pattern = "-", replacement = ".", x = genes.plot))
    # )
    data.to.plot$pct.exp[data.to.plot$pct.exp < dot.min] <- NA
    data.to.plot$pct.exp <- data.to.plot$pct.exp * 100
    data.to.plot <- merge(data.to.plot, markers, by.x = "genes.plot", by.y = "Markers")
    
    return(data.to.plot)
}



#### make dot plot after singleR
# meta_after <- read.csv(here("02_rds/meta_after_singleR.csv"),row.names = 1, stringsAsFactors = F)
# obj@meta.data$cell_name_bk <- obj@meta.data$cell_name
# obj@meta.data$cell_name = as.character(obj@meta.data$cell_name)
# obj@meta.data$cell_name <- as.character(meta[rownames(obj@meta.data), "cell_name1"])

meta <- readRDS("11_CNV/meta.rds")
meta$cell_short[meta$cell_short == "ATII"] = "AT2"
meta$cell_short[meta$cell_short == "Mo"] = "Mφ"
obj@meta.data = meta

data = ModifiedDotData(obj, markers = markers, group.by = "cell_short")

required_immu <- c(
    "B", "CD4", "CD8",
    "DC", "Epi", 
    "Gran", "Mast",
    "Mo", "Mφ", "NK",
    "Tregs"  # "EC", 
)

data <- data[as.character(data$id) %in% required_immu, ]

cairo_pdf(here("01_first_plots/dotplot_immu_after_singleR.pdf"), width = 8, height = 10)
ggplot(data, aes(x=id, y=genes.plot, size = pct.exp, color=avg.exp.scale)) +
    geom_point() +
    scale_radius(range = c(0, 6), limits = c(NA, NA)) +
    scale_color_gradient(low = "lightgrey", high = "blue") +
    facet_grid(Cells~., scales = "free", space = "free") +
    theme_bw(base_family = "Arial Unicode MS") +
    theme(
        axis.text = element_text(size = 15),
        strip.text.y = element_text(size = 15),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 15),
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
    ) +
    labs(x="", y="", color = "avg.exp")
dev.off()



# markers = read.xlsx(here("20190701_gene_markers.xlsx"), sheet = 3)
# markers <- na.omit(markers[, 1:2])
# markers <- unique(markers)
# temp_markers = markers[order(markers$Cells, markers$Markers), ]
# temp_markers$Cells = short_names[temp_markers$Cells, 1]
# temp_markers <- na.omit(temp_markers)
# 
# temp_markers <- temp_markers[temp_markers$Cells %in% cell_types, ]
# 
# 
# temp_markers <- temp_markers[!temp_markers$Cells %in% immu, ]


# ggsave("LungCancer10x/01_first_plots/dotplot_08.pdf", plot = p, width = 12, height = 30, units = "in", dpi = 600)
obj@meta.data <- meta
data = ModifiedDotData(obj, markers = temp_markers[!temp_markers$Cells %in% immu, ], group.by = "cell_short")

required_immu <- c(
    "B", "CD4", "CD8",
    "DC", 
    "Gran", "Mast",
    "Mo", "Mφ", "NK",
    "Tregs"  # "EC", 
)

data <- data[!as.character(data$id) %in% required_immu, ]
data$id = as.character(data$id)

cairo_pdf(here("01_first_plots/dotplot_non_immu.pdf"), width = 8, height = 12)
ggplot(data, aes(x=as.character(id), y=genes.plot, size = pct.exp, color=avg.exp.scale)) +
    geom_point() +
    scale_radius(range = c(0, 6), limits = c(NA, NA)) +
    scale_color_gradient(low = "lightgrey", high = "blue") +
    facet_grid(Cells~., scales = "free", space = "free") +
    theme_bw(base_family = "Arial Unicode MS") +
    theme(
        axis.text = element_text(size = 15),
        strip.text.y = element_text(size = 15),
        panel.grid = element_blank(),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 15),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
    ) +
    labs(x="", y="", color = "avg.exp")
dev.off()


## Dotplot of non_immu
temp_markers = read.xlsx("LungCancer10x/markers.xlsx", sheet = 1)

colnames(temp_markers) <- c("Cells", "Markers", "Type")
data = ModifiedDotData(obj, markers = temp_markers, group.by = "cell_short")
data$Cells = factor(data$Cells, levels = unique(temp_markers$Cells))


required_immu <- c(
    "B", "CD4", "CD8",
    "DC", "Gran", "Mast",
    "Mo", "Mφ", "NK",
    "Tregs"  # "EC", 
)


for (i in unique(data$Type)) {

    temp = data[data$Type == i, ]
    print(i)
    print(nrow(temp))

    p <- ggplot(temp[!as.character(temp$id) %in% required_immu, ], aes(x=id, y=genes.plot, size = pct.exp, color=avg.exp.scale)) +
        geom_point() +
        scale_radius(range = c(0, 6), limits = c(NA, NA)) +
        scale_color_gradient(low = "lightgrey", high = "blue") +
        facet_grid(Cells~., scales = "free", space = "free") +
        theme_bw(base_family = "Arial Unicode MS") +
        theme(
            axis.text = element_text(size = 15),
            strip.text.y = element_text(size = 15, angle = 0),
            legend.title = element_text(size = 20),
            legend.text = element_text(size = 15),
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
        ) +
        labs(x="", y="", color = "avg.exp") +
        scale_size_continuous(range = c(0.1,4))
    
    
    
    ggsave(
        filename = here(paste0("01_first_plots/dotplot_", i, "_bioarx.pdf")),
        plot = p, width = 8, height = 1 + nrow(temp)  / 50, device = cairo_pdf
    )
    
}


data$id = as.character(data$id)
data$id[data$id == "Epi"] = "AT1"

p <- ggplot(data[!as.character(data$id) %in% required_immu, ], aes(y=id, x=genes.plot, size = pct.exp, color=avg.exp.scale)) +
    geom_point() +
    scale_radius(range = c(0, 6), limits = c(NA, NA)) +
    scale_color_gradient(low = "lightgrey", high = "blue") +
    facet_grid(.~Cells, scales = "free", space = "free") +
    theme_bw(base_family = "Arial Unicode MS") +
    theme(
        axis.text = element_text(size = 15),
        strip.text = element_text(size = 15, angle = 90),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 15),
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    ) +
    labs(x="", y="", color = "avg.exp") +
    scale_size_continuous(range = c(0.1,4))



pdf(here(paste0("01_first_plots/dotplot_non_immu_bioarx.pdf")), height = 6, width = 21)
g <- ggplot_gtable(ggplot_build(p))
stripr <- which(grepl('strip-t', g$layout$name))
fills <- c(rep("#0084D1", 10), 
           rep("grey75", 4), 
           rep("#A0C807", 7))
k <- 1
for (i in stripr) {
    j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
    g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
    k <- k+1
}
grid::grid.draw(g)
dev.off()



temp_markers = read.xlsx("LungCancer10x/markers.xlsx", sheet = 2)

colnames(temp_markers) <- c("Cells", "Markers")
data = ModifiedDotData(obj, markers = temp_markers, group.by = "cell_short")

required_immu <- c(
    "B", "CD4", "CD8",
    "DC", "Gran", "Mast",
    "Mo", "Mφ", "NK",
    "Tregs"  # "EC", 
)

data$Cells = factor(data$Cells, levels = unique(temp_markers$Cells))
data$id = factor(data$id, levels = c("B", "CD4", "CD8", "Tregs", "NK", "DC", "Mast", "Gran", "Mo", "Mφ"))


p <- ggplot(data[as.character(data$id) %in% required_immu, ], aes(y=id, x=genes.plot, size = pct.exp, color=avg.exp.scale)) +
    geom_point() +
    scale_radius(range = c(0, 6), limits = c(NA, NA)) +
    scale_color_gradient(low = "lightgrey", high = "blue") +
    facet_grid(.~Cells, scales = "free", space = "free") +
    theme_bw(base_family = "Arial Unicode MS") +
    theme(
        axis.text = element_text(size = 15),
        strip.text.y = element_text(size = 15, angle = 0),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 15),
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        strip.text = element_text(angle = 90, size = 15)
    ) +
    labs(x="", y="", color = "avg.exp")

p

ggsave(
    filename = here("01_first_plots/dotplot_immu_bioarx.pdf"),
    plot = p, height = 6, width = 16, device = cairo_pdf
)