library(openxlsx)
library(stringr)
library(fmsb)
library(Seurat)
library(reshape2)
library(dplyr)
library("wesanderson")

set.seed(1)




format_radar_matrix <- function(object, meta) {
    meta$genes = as.character(meta$genes)

    mat = MinMax(object@scale.data[as.character(meta$genes), , drop = FALSE], -2.5, 2.5)
    mat = melt(as.matrix(mat))
    colnames(mat) <- c("genes", "cell", "value")

    # print(head(mat)) 
    # print(head(meta))

    expr <- merge(meta, mat, by = "genes")
    expr$PatientID = object@meta.data[expr$cell, "PatientID"]
    
    # print(head(expr))
    
    data = as.data.frame(
        expr %>% 
            dplyr::select(stage, PatientID, value) %>% 
            dplyr::group_by(stage, PatientID) %>%
            mutate(val=mean(value)) %>%
            dplyr::select(stage, PatientID, val) %>%
            unique()
        )
    
    data = dcast(data, stage~PatientID, value.var = "val")
    rownames(data) <- data$stage
    data <- data[, colnames(data) != "stage"]
    
    return(as.data.frame(t(data)))
}


make_radar_plots <- function(object, meta, group.by = "Stage", output_prefix=NULL) {
    
    current_stages = sapply(meta$stage, function(x) {
        return(str_split(x, "\\.")[[1]][2])
    })
    current_stages = intersect(as.character(object@meta.data[, group.by]), current_stages)

    for(i in current_stages) {
        print(i)
        temp_meta = object@meta.data[as.character(object@meta.data[, group.by]) == i, ]
        obj <- CreateSeuratObject(
            object@raw.data[, rownames(temp_meta), drop = FALSE],
            meta = temp_meta
        )
        
        obj@scale.data = object@scale.data[, rownames(temp_meta), drop = FALSE]
        
        expr = format_radar_matrix(obj, meta)    
        expr <- rbind(
            rep(ceiling(max(expr)), ncol(expr)),
            rep(floor(min(expr)), ncol(expr)),
            expr
        )
    
        # print(expr)
        if (!is.null(output_prefix)) {
            png(paste(output_prefix, "_", i, "_radar.png", sep = ""), width = 12, height = 6, res = 600, units = "in")
        }
        
        par(
            mfrow=c(1, 2), 
            bty='n', 
            oma = c(0.5,0.5,0,0) + 0.1,
            mar = c(1,0,0,1) + 0.1
        )
    
        legend_labels = rownames(expr)
        legend_labels = legend_labels[3: length(legend_labels)]
        colors_border = wes_palette("Zissou1", length(legend_labels), type = "continuous")
    
        radarchart(expr, axistype=1, 
            #custom polygon
            pcol=colors_border, # pfcol=colors_in, plwd=4 , plty=1, 
            #custom the grid
            cglcol="grey", cglty=1, axislabcol="grey", 
            caxislabels=seq(expr[2, 1], expr[1, 1], (expr[1, 1] - expr[2, 1]) / 5), cglwd=0.8,
            #custom labels
            vlcex=0.8 
        )
        plot(0, col="white", cex=0, axes=F, ann=FALSE)
    
        legend(
            "left",
            legend = legend_labels, 
            bty = "n", pch=20 , 
            col=colors_border , 
            text.col = "grey", 
            cex=1.2, 
            pt.cex=3,
            ncol = 2
        )
        
        if (!is.null(output_prefix)) {
            dev.off()
        }
    }
}



args = commandArgs(trailingOnly = T)
module = read.xlsx(args[1], sheet = 2)
obj = readRDS(args[2])


temp_module = module[module$Cell_name == args[4], ]

# print(head(temp_module))

# format meta
meta = NULL
for(j in 1:nrow(temp_module)) {
    stage = paste("M", temp_module[j, "Stage"], sep = ".")
    genes = str_split(temp_module[j, "Genes"], "\\|")[[1]]
    
    meta = rbind(meta, data.frame(stage=stage, genes=genes))
}

# print(meta)

# two random select groups
for(i in 1:2) {
    meta = rbind(meta, data.frame(
        stage = paste("R", i, sep = ""), 
        genes = sample(
            rownames(obj@raw.data)[!rownames(obj@raw.data) %in% meta$genes],
            min(100, sum(!rownames(obj@raw.data) %in% meta$genes))
        )
    ))

}

# format and plot
make_radar_plots(object = obj, meta = meta, group.by = "Stage", output_prefix = args[3])
