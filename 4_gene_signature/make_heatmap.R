# scripts to make heatmap
args = commandArgs(trailingOnly = TRUE)


rds = args[1]
genes = strsplit(args[2], ",")
outfile = args[3]
title = args[4]
group = args[5]


suppressMessages(function() {
    library(Seurat)
    library(ggplot2)
})


obj <- readRDS(rds)


p <- DoHeatmap(
    obj, 
    data.use = NULL, 
    genes.use = genes, 
    disp.min = -2.5, 
    disp.max = 2.5, 
    group.by = group,
    group.order = NULL, 
    draw.line = TRUE, 
    col.low = "#FF00FF",
    col.mid = "#000000", 
    col.high = "#FFFF00", 
    slim.col.label = TRUE,
    remove.key = TRUE, 
    rotate.key = FALSE, 
    title = title, 
    cex.col = 10,
    cex.row = 10, 
    group.label.loc = "bottom", 
    group.label.rot = FALSE,
    group.cex = 15, 
    group.spacing = 0.15, 
    assay.type = "RNA",
    do.plot = FALSE
)

height = length(genes) / 8

if (height < 5) {
    height = 5
} else if (height > 40) {
    height = 40
}

ggsave(filename = p, plot = p, width = 6, height = height, units = "in", dpi = 600)

