#!/usr/bin/env Rscript
# create at 2019.05.29
# author zhangyiming

load <- function() {
    library(Seurat)
    library(Mfuzz)
    library(openxlsx)
    library(mosaic)
    library(ggplot2)
}

suppressPackageStartupMessages(load())

args = commandArgs(trailingOnly = T)

root.dir = args[1]
rds = args[2]


root.dir = "/mnt/raid62/Lung_cancer_10x/02_figures_each_cell/Basal_SCC/"
rds = "/mnt/raid62/Lung_cancer_10x/02_figures_each_cell/Basal_SCC/Basal.rds"

set.seed(1)
setwd(root.dir)


cluster_markers = "annotation_results_by_cluster.xlsx"
cluster_markers = read.xlsx(paste(root.dir, cluster_markers, sep = "/"))

obj <- readRDS(rds)

# function to read and format the normalized_counts from sctransform
read_sctransform <- function(path="paga/normalized_counts.csv.gz") {
    r = gzfile(path)
    data = read.csv(r, row.names = 1)
    
    colnames(data) = gsub("\\.", "-", colnames(data), perl=F)
    colnames(data) = gsub("^X", "", colnames(data), perl=T)
    return(data)
}


# function to select cluster_markers by qvalue and logfc
# :param cluster_markers cluster markers from Seurat FindMarkers
# :param qvalue: q value threshold
# :param logfc: logFC threshold
# :return vector of gene names
select_cluster_markers <- function(cluster_markers, qvalue=0.05, logfc=0.5) {
    return(cluster_markers[cluster_markers$p_val_adj < qvalue & cluster_markers$avg_logFC, "gene"])
}


# function to make heatmaps
# :param obj Seurat obj
# :param cluster results from perform_mfuzz
# :param output_prefix: the output file prefix
make_heatmap <- function(obj, cluster, output_prefix) {
    for (i in sort(unique(cluster))) {
        temp_genes = names(cluster[cluster == i])
        
        p <- DoHeatmap(
            obj, 
            group.by = "Stage", 
            genes.use = temp_genes, 
            cex.col=0,
            slim.col.label = TRUE, 
            remove.key = TRUE,
            do.plot = F,
            title = i
        )
        
        height = length(temp_genes) / 8
        
        if (height > 40) {
            height = 40
        }
        
        ggsave(
            paste0(output_prefix, i, ".png"),
            p,
            width = 12,
            height = height,
            dpi = 300,
            units = "in"
        )
    }
}


# function to perform mfuzz on expression based on cluster selection
# :param expr dataframe, the expression matrix
# :param cluster_markers qvlaue see select_cluster_markers
# :param logfc see select_cluster_markers
# :param init_cluster init cluster
# :param m see mfuzz
perform_mfuzz <- function(expr, cluster_markers, qvalue=0.05, logfc=0.5, init_cluster=9, m=NULL ) {
    
    selected_markers = select_cluster_markers(cluster_markers, qvalue, logfc)
    
    # filter genes and create ExpressionSet
    expr = ExpressionSet(
        as.matrix(expr[intersect(rownames(expr), selected_markers),])
    )
    expr = standardise(expr)
    
    # best m, if m is null, then estimate based on expr
    if (is.null(m)) {
        m = mestimate(expr)
    }
    cl = mfuzz(expr, c=init_cluster, m=m)
    
    return(cl)
}



# function to construct matrix of mean zscore
construct_stage_group_zscore <- function(obj, expr, cluster) {
    groups = sort(unique(cluster$cluster))
    stages = sort(unique(obj@meta.data$Stage))
    
    res = as.data.frame(matrix(0, nrow = length(groups), ncol = length(stages)))
    rownames(res) = groups
    colnames(res) = stages
    
    for (i in groups) {
        temp_gene = names(cluster$cluster[cluster$cluster == i])
        
        for (j in stages) {
            # print(paste0(i, j))
            temp_cells = rownames(obj@meta.data[obj@meta.data$Stage == j, ])
            
            temp_expr = expr[intersect(temp_gene, rownames(expr)), temp_cells]
            # temp_expr = scale(temp_expr)
            
            res[i, j] = mean(apply(temp_expr, 1, sum))
        }
    }
    
    res[is.na(res)] = 0
    return(res)
}




expr = read_sctransform()

cl = perform_mfuzz(expr, cluster_markers)

mtx = construct_stage_group_zscore(obj, expr, cluster=cl)

temp = hclust(dist(mtx))

pdf("test.pdf")
plot(temp)
dev.off()

