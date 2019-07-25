library(Mfuzz)
library(openxlsx)
library(ComplexHeatmap)
library(reticulate)
library("wesanderson")
# if(!require("Cairo")) {
#     install.packages("Cairo")
#     require("Cairo")
# }


args = commandArgs(trailingOnly = T)


input_dir = args[1]
output_dir = args[2]

print(args)

dir.create(output_dir, showWarnings = F, recursive = T)

## read markers
markers <- read.xlsx(paste(input_dir, "markers_by_stage.xlsx", sep = "/"), rowNames = T)
markers <- markers[markers$p_val_adj < 0.05 & markers$avg_logFC > 0.5, ]

nms = unique(markers$gene)
ig_genes = c(grep("^IGJ", nms, v=T), 
             grep("^IGH",nms,v=T),
             grep("^IGK", nms, v=T), 
             grep("^IGL", nms, v=T))

bad_genes = unique(c(grep("^MT-", nms, v=T), grep("^MTMR", nms, v=T), grep("^MTND", nms, v=T),"NEAT1","TMSB4X", "TMSB10", ig_genes))

markers <- markers[!markers$gene %in% bad_genes, ]

## read expression data
obj <- readRDS(paste(input_dir, "seurat.rds", sep = "/"))


## run mfuzz
if(file.exists(paste(output_dir, "results.csv", sep = "/"))) {
# if (FALSE) {
    res = read.csv(paste(output_dir, "results.csv", sep = "/"))
} else {
    expr = ExpressionSet(
        as.matrix(obj@scale.data[unique(markers$gene),])
    )
    
    m = mestimate(expr)
    
    temp = Dmin(expr, m=m, crange=seq(2,10,1), repeats=3, visu=FALSE)
    
    if (sum(!is.na(temp)) == 0) {
        c = length(unique(obj@meta.data$Stage))
    } else {
        new_temp = temp[!is.na(temp)]
        
        # using kneed to selected best cluster
        kee = import("kneed")
        c = kee$KneeLocator(
            1:length(new_temp),
            new_temp, 
            curve='convex', 
            direction='decreasing'
        )
        c = c$knee
        
        c = which(temp == new_temp[c])
    }
 
    
    cl <- mfuzz(expr, c=c, m=m)
    
    
    # extract results
    res = as.data.frame(cl$cluster)
    colnames(res) <- "Clt"
    res$gene = rownames(res)
    
    res = res[order(res$Clt), ]
    
    write.csv(res, paste(output_dir, "results.csv", sep = "/"))
}

print(head(res))

# make heatmap

meta = obj@meta.data
meta = meta[order(meta$Stage), ]

gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
}

patient_colors = gg_color_hue(35)
names(patient_colors) <- c(
    sapply(1:23, function(x){sprintf("PA%02d", x)}),
    sapply(1:12, function(x){sprintf("PS%02d", x)})
)



ra_col = wes_palette("Zissou1", length(unique(res$Clt)), type = "continuous")
print(length(unique(res$Clt)))
names(ra_col) <- 1:length(unique(res$Clt))

ra = rowAnnotation(
    Mfuzz=res$Clt,
    col = list(
        Mfuzz=ra_col
    )
)

ha = HeatmapAnnotation(
    Stage=meta$Stage,
    Patient=meta$PatientID,
    col=list(
        Stage=c("I"="#EC7A21", "II"="#D73F47", "III"="#65A9A3", "IV"="#4A933E", "LUAD_Normal"="#FECC1B", "LUSC_Normal"="#DA6906"),
        Patient=patient_colors
    ),
    border = TRUE
)

pdf(paste(output_dir, "mfuzz.pdf", sep = "/"), width = 15, height = 15)
Heatmap(
    obj@scale.data[res$gene, rownames(meta)], 
    cluster_rows = F, 
    cluster_columns = F, 
    show_row_names = F, 
    show_column_names = F,
    name="Expr",
    bottom_annotation = ha,
    left_annotation = ra
)
dev.off()
