library(Seurat)
library(openxlsx)
library(clusterProfiler)
library(DOSE)
library(KEGG.db)
library(GO.db)
library(org.Hs.eg.db)
library(wesanderson)
library(ggrepel)
library(dplyr)
library(stringr)





make_voca_plot <- function(data, genes.label, genes.use = NULL) {

    if (is.null(genes.use)) {
        genes.use = unique(data$gene)
    } 

    data = data[data$gene %in% genes.use, ]

    data$ident[
        data$avg_logFC > 0.5 & data$p_val_adj < 0.05
    ] = "1"

    data$ident[
        data$avg_logFC < -0.5 & data$p_val_adj < 0.05
    ] = "2"

    data$ident[
        abs(data$avg_logFC) <= 0.5 | data$p_val_adj >= 0.05
    ] = "3" 

    # data$ident <- apply(data, 1, function(row) {
    #     if (row[2] > 0.25 && row[5] < 0.05) {
    #         return("1")
    #     } else if (row[2] < -0.25 && row[5] < 0.05) {
    #         return("2")
    #     } 
        
    #     return("3")
    # })
    
    print(unique(data$ident))
    
    data$p_val_adj <- -log10(data$p_val_adj)
    # data$gene <- rownames(data)
    data$ident = as.factor(data$ident)
    print(sum(data$gene %in% genes.label))
    p <- ggplot(data = data, aes(x=avg_logFC, y=p_val_adj, color = ident)) + 
        geom_point() + 
        scale_color_manual(values=c(
            "1"="red", 
            "2"="green", 
            "3"="grey50"
        )) +
        geom_text_repel(
            aes(label=ifelse(gene %in% genes.label, as.character(gene),''))
        ) + theme_bw() +
        theme(legend.position = "none") +
        labs(x="log2(FC)", y="-log10(FDR)")

    return(p)
}



make_compare_dotplot_between_cluster <- function(object, genes.label, scale=FALSE, group.by = "res.0.8", cluster = 1, genes.use=NULL) {
    meta = object@meta.data

    if (is.null(genes.use)) {
        genes.use = rownames(object@raw.data)
    }
    
    if(scale) {
        data = as.matrix(object@scale.data[ ,rownames(meta)])
    } else {
        data = as.matrix(object@raw.data[, rownames(meta)])
    }
    
    temp = as.data.frame(matrix(NA, nrow = nrow(data), ncol = 3))
    
    colnames(temp) <- c("gene", cluster, "Rest")
    
    temp[, 1] = rownames(data)
    temp[, cluster] = apply(data[, rownames(meta[meta[, group.by] == cluster,])], 1, function(x) {
        return(mean(as.numeric(x)))
    })
    
    temp[, "Rest"] = apply(data[, rownames(meta[meta[, group.by] != cluster,])], 1, function(x) {
        return(mean(as.numeric(x)))
    })
    
    
    temp[,2] = log2(temp[, 2] + 1)
    temp[,3] = log2(temp[, 3] + 1)
    
    max_ = ceiling(max(temp[, 2:3])) + 1
    min_ = floor(min(temp[, 2:3]))
    
    temp$gene[!temp$gene %in% genes.label] = ""
    
    p <- eval(
        parse(text = paste("ggplot(data = temp, aes(x =", cluster, ", y = Rest, label = gene))"))   
    )
    
    p = p + 
        geom_point(color=ifelse(temp$gene %in% genes.label, "red", 'grey50')) +
        xlim(min_, max_) +
        ylim(min_, max_) +
        geom_abline(
            color="red", 
            linetype="dashed"
        ) +
        # geom_text(
        #     aes(label=ifelse(gene %in% genes.label, as.character(gene),''), color = "red"),
        #     hjust=-0.2,
        #     vjust=0.5,
        #     angle = 0
        # ) +
        geom_text_repel(
            aes(label=ifelse(gene %in% genes.label, as.character(gene),''))
        ) +
        theme(legend.position = "none") +
        labs(
            x = paste0("-log2(", cluster, "); n=", sum(meta[, group.by] == cluster)),
            y = paste0("-log2(Rest); n=", sum(meta[, group.by] != cluster))
        )
    return(p)
    
}



args = commandArgs(trailingOnly = T)

input_rds = args[1]
output_dir = args[2]

obj <- readRDS(input_rds)

if (length(unique(obj@meta.data$res.0.8)) < 2) {
    quit("no")
}

dir.create(output_dir, showWarnings = F, recursive = T)

output_xlsx = paste(dirname(input_rds), "markers_by_stage.xlsx", sep = "/")

markers = read.xlsx(output_xlsx, sheet = 1, rowNames = T)



nms = rownames(obj@raw.data)
ig_genes = c(grep("^IGJ", nms, v=T), 
                grep("^IGH",nms,v=T),
                grep("^IGK", nms, v=T), 
                grep("^IGL", nms, v=T))

bad_genes = unique(c(grep("^MT-", nms, v=T), grep("^MTMR", nms, v=T), grep("^MTND", nms, v=T),"NEAT1","TMSB4X", "TMSB10", ig_genes))

markers = markers[!rownames(markers) %in% bad_genes,]
markers_de <- markers[markers$p_val_adj < 0.05 & abs(markers$avg_logFC) > 0.5, ]


genes.use = rownames(obj@raw.data)[!rownames(obj@raw.data) %in% bad_genes]


for (i in unique(markers$ident)) {
    
    temp_genes_use = markers[markers$ident == i & markers$gene %in% genes.use & markers$p_val_adj < 0.05, ]

    genes.label = c()
    temp_genes_use <- temp_genes_use[order(temp_genes_use$avg_logFC, decreasing = T), ]
    genes.label = c(temp_genes_use$gene[1:10])

    temp_genes_use <- temp_genes_use[order(temp_genes_use$avg_logFC, decreasing = F), ]
    genes.label = c(genes.label, temp_genes_use$gene[1:10])

    p <- make_compare_dotplot_between_cluster(
        object = obj, 
        scale = F, 
        genes.label = genes.label,
        genes.use = genes.use,
        cluster = i,
        group.by = "Stage"
    )

    ggsave(
        filename = paste(output_dir, paste0(i, "_gene.pdf"), sep = "/"),
        plot = p,
        width = 6,
        height = 6,
        dpi = 600,
        units = "in"
    ) 


    p <- make_voca_plot(markers[markers$ident == i & markers$gene %in% genes.use, ], genes.label = genes.label, genes.use = genes.use)

    ggsave(
        filename = paste(output_dir, paste0(i, "_volca.pdf"), sep = "/"),
        plot = p,
        width = 6,
        height = 6,
        dpi = 600,
        units = "in"
    ) 
}

