library(wesanderson)
library(Seurat)
library(ggrepel)
library(openxlsx)
library(ggplot2)


make_compare_dotplot_between_cluster <- function(object, genes.label, scale=FALSE, group.by = "res.0.6", cluster = 1) {
    meta = object@meta.data
    
    if(scale) {
        data = as.matrix(object@scale.data[ ,rownames(meta)])
    } else {
        data = as.matrix(object@raw.data[, rownames(meta)])
    }
    
    temp = as.data.frame(matrix(NA, nrow = nrow(data), ncol = 3))
    
    colnames(temp) <- c("gene", paste0("Cluster", cluster), "Rest")
    
    temp[, 1] = rownames(data)
    temp[, paste0("Cluster", cluster)] = apply(data[, rownames(meta[meta[, group.by] == cluster,])], 1, function(x) {
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
        parse(text = paste("ggplot(data = temp, aes(x =", paste0("Cluster", cluster), ", y = Rest, label = gene))"))   
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
        geom_text_repel() +
        theme(legend.position = "none") +
        labs(
            x = paste0("-log2(", paste0("Cluster", cluster), ")"),
            y = paste0("-log2(Rest)")
        )
    return(p)
    
}


args = commandArgs(trailingOnly = T)

xlsx = paste(dirname(args[1]), "annotation_results_by_cluster.xlsx", sep = "/")
obj <- readRDS(args[1])

if (file.exists(xlsx) && length(unique(obj@meta.data$res.0.6 > 1))) {
    
    output = args[2]
    
    dir.create(output, showWarnings = F, recursive = T)
    
    clt = read.xlsx(xlsx)
    
    
    extract_gene_from_cluster <- function(markers, topn = 10, ident = 1, two.sides = TRUE, pval = 0.05) {
        markers = markers[markers$p_val_adj < pval & markers$ident == ident, ]
        
        markers = markers[order(markers$avg_logFC, decreasing = T), ]
        res = markers$gene[1:topn]
        
        if(two.sides) {
            markers = markers[order(markers$avg_logFC, decreasing = F), ]
            res = c(res, markers$gene[1:topn])
        }
        
        return(res)
    }
    
    
    for(i in unique(obj@meta.data$res.0.6)) {
        p <- make_compare_dotplot_between_cluster(
            obj, 
            genes.label = extract_gene_from_cluster(clt, ident = i, two.sides = F), 
            group.by = "res.0.6", 
            cluster = i
        )
        
        ggsave(
            filename = paste(output, paste0(i, ".png"), sep = "/"),
            plot = p,
            width = 6,
            height = 6,
            units = "in",
            dpi = 600
        )
    }
    
    
}

