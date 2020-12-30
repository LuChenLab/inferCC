rm(list=ls())

options(stringsAsFactors = F)

set.seed(1)

#!/usr/bin/env Rscript
# create at 2019.05.29
# author zhangyiming

load <- function() {
    library(Seurat)
    library(Mfuzz)
    library(openxlsx)
    library(ggplot2)
    library(clusterProfiler)
    library(DOSE)
    library(GO.db)
    library(KEGG.db)
    library(dplyr)
    library(reshape2)
}

suppressPackageStartupMessages(load())

root.dir = "/mnt/raid62/Lung_cancer_10x/02_figures_each_cell/Mast_SCC"
rds = "/mnt/raid62/Lung_cancer_10x/02_figures_each_cell/Mast_SCC/Mast.rds"

dir.create(paste(root.dir, "cluster_gene_module", sep = "/"), showWarnings = F)

setwd(root.dir)
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
    return(cluster_markers[cluster_markers$p_val_adj < qvalue & cluster_markers$avg_logFC > logfc, "gene"])
}


# function to make heatmaps
# :param obj Seurat obj
# :param cluster results from perform_mfuzz
# :param output_prefix: the output file prefix
make_heatmap <- function(obj, cluster, output_prefix=NULL) {
    for (i in sort(unique(cluster))) {
        temp_genes = names(cluster[cluster == i])
        
        p <- DoHeatmap(
            obj, 
            group.by = "res.0.6", 
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
        } else if (height < 5) {
            height = 5
        }
        
        if (is.null(output_prefix)) {
            print(p)
        } else {
            ggsave(
                paste(output_prefix, i, ".png", sep = ""),
                p,
                width = 12,
                height = height,
                dpi = 300,
                units = "in"
            )
        }
    }
}



# function to make dotplot
# :param obj Seurat obj
# :param cluster results from perform_mfuzz
# :param output_prefix: the output file prefix
make_dotplot <- function(obj, cluster, output_prefix=NULL) {
    for (i in sort(unique(cluster))) {
        temp_genes = names(cluster[cluster == i])
        
        p <- DotPlot(
            obj, 
            group.by = "res.0.6", 
            genes.plot = temp_genes, 
            x.lab.rot = TRUE, 
            do.return = T
        ) + labs(title = i) + coord_flip()
        
        height = length(temp_genes) / 8
        
        if (height > 40) {
            height = 40
        } else if (height < 5) {
            height = 5
        }
        
        if (is.null(output_prefix)) {
            print(p)
        } else {
            ggsave(
                paste(output_prefix, i, ".png", sep = ""),
                p,
                width = 12,
                height = height,
                dpi = 300,
                units = "in"
            )
        }
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
    print(dim(expr))
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
    stages = sort(unique(obj@meta.data$res.0.6))
    
    res = as.data.frame(matrix(0, nrow = length(groups), ncol = length(stages)))
    rownames(res) = groups
    colnames(res) = stages
    
    for (i in groups) {
        temp_gene = names(cluster$cluster[cluster$cluster == i])
        
        for (j in stages) {
            # print(paste0(i, j))
            temp_cells = rownames(obj@meta.data[obj@meta.data$res.0.6 == j, ])
            
            temp_expr = expr[intersect(temp_gene, rownames(expr)), temp_cells]
            # temp_expr = scale(temp_expr)
            
            res[i, j] = mean(apply(temp_expr, 1, sum))
        }
    }
    
    res[is.na(res)] = 0
    return(res)
}




do_kegg <- function(eg, cluster=NA, pvalueCutoff = 0.05) {
    kk <- enrichKEGG(gene     = eg$ENTREZID,
                     organism     = 'hsa',
                     pvalueCutoff = pvalueCutoff)
    kk <- as.data.frame(kk)
    
    if (nrow(kk) == 0) {
        return(kk)
    }
    kk$cluster <- cluster
    return(kk)
}


do_go <- function(eg, cluster = NA, pvalueCutoff = 0.01, qvalueCutoff = 0.05, cutoff=0.7) {
    
    res = NULL
    # for(i in c("BP", "CC", "MF")) {
    #     ego <- enrichGO(gene      = eg$ENTREZID,
    #                     keyType       = 'ENTREZID',
    #                     OrgDb         = org.Hs.eg.db,
    #                     ont           = i,
    #                     pAdjustMethod = "BH",
    #                     pvalueCutoff  = pvalueCutoff,
    #                     qvalueCutoff  = qvalueCutoff,
    #                     readable      = FALSE)
    #     
    #     if(is.null(ego)) {
    #         return(ego)
    #     }
    #     ego <- clusterProfiler::simplify(ego, cutoff=cutoff, by="p.adjust", select_fun=min)
    #     
    #     ego = as.data.frame(ego)
    #     
    #     if (nrow(ego) > 0) {
    #         ego$ONTOLOGY = i
    #     }
    #     
    #     if (is.null(res)) {
    #         res = ego
    #     } else {
    #         res = rbind(res, ego)
    #     }
    # }
    # 
    # if (nrow(res) == 0) {
    #     return(res)
    # }
    
    ego <- enrichGO(gene      = eg$ENTREZID,
                    keyType       = 'ENTREZID',
                    OrgDb         = org.Hs.eg.db,
                    ont           = "ALL",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = pvalueCutoff,
                    qvalueCutoff  = qvalueCutoff,
                    readable      = FALSE)
    res = as.data.frame(ego)
    
    if (is.null(res) || nrow(res) <= 0) {
        return(res)
    }
    
    res$cluster <- cluster
    
    return(res)
}


do_do <- function(eg, cluster = NA, pvalueCutoff = 0.05, qvalueCutoff = 0.05, minGSSize = 5, maxGSSize = 500) {
    do <- enrichDO(gene     = eg$ENTREZID,
                   ont           = "DO",
                   pvalueCutoff  = pvalueCutoff,
                   pAdjustMethod = "BH",
                   minGSSize     = minGSSize,
                   maxGSSize     = maxGSSize,
                   qvalueCutoff  = qvalueCutoff,
                   readable      = TRUE)
    
    do = as.data.frame(do)
    
    if (nrow(do) == 0) {
        return(do)
    }
    do$cluster = cluster
    
    return(do)
} 


get_entrzid <- function(markers) {
    res = NULL
    for(i in unique(markers)) {
        print(i)
        
        temp <- names(markers[markers == i])
        eg <- bitr(unique(temp), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
        
        if(!is.null(eg) && nrow(eg) > 0) {
            eg$ident = i
            
            if(is.null(res)) {
                res = eg
            } else {
                res = rbind(res, eg)
            }
        }
    }
    
    return(res)
}

print(getwd())
print(rds)

cluster_markers = "annotation_results_by_cluster.xlsx"
cluster_markers = read.xlsx(paste(root.dir, cluster_markers, sep = "/"))

expr = read_sctransform()
obj <- readRDS(rds)
qvalue = 0.01
logfc = 0

cl = perform_mfuzz(expr, cluster_markers, logfc=logfc, qvalue=qvalue)

mtx = construct_stage_group_zscore(obj, expr, cluster=cl)

temp = hclust(dist(mtx))
# make_heatmap(obj, cl$cluster)
plot(temp)
# new_cl = c()

new_cl = perform_mfuzz(expr, cluster_markers, logfc=logfc, qvalue=qvalue, init_cluster = 3)
new_cl = new_cl$cluster
# make_heatmap(obj, new_cl)
make_heatmap(obj, new_cl, output_prefix = paste(root.dir, "cluster_gene_module/heatmap_", sep = "/"))
make_dotplot(obj, new_cl, output_prefix = paste(root.dir, "cluster_gene_module/dotplot_", sep = "/"))
eg = get_entrzid(new_cl)
kk <- NULL

for(i in unique(eg$ident)) {
    print(i)
    temp = do_kegg(eg[eg$ident == i, ])
    
    if(is.null(temp) || nrow(temp) == 0) {
        next
    }
    
    temp$ident = i
    if(is.null(kk)) {
        kk = temp
    } else {
        kk = rbind(kk, temp)
    }
}
go <- NULL

for(i in unique(eg$ident)) {
    print(i)
    temp = do_go(eg[eg$ident == i, ])
    
    if(is.null(temp) || nrow(temp) == 0) {
        next
    }
    
    temp$ident = i
    if(is.null(go)) {
        go = temp
    } else {
        go = rbind(go, temp)
    }
}
do <- NULL

for(i in unique(eg$ident)) {
    print(i)
    temp = do_do(eg[eg$ident == i, ])
    
    if(is.null(temp) || nrow(temp) == 0) {
        next
    }
    
    temp$ident = i
    if(is.null(do)) {
        do = temp
    } else {
        do = rbind(do, temp)
    }
}
print(getwd())

temp = as.data.frame(new_cl)
colnames(temp) <- "gene_module_id"
temp$gene = rownames(temp)

wb = createWorkbook()
addWorksheet(wb, "gene_module")
writeData(wb, 1, temp)

addWorksheet(wb, "kegg")
writeData(wb, 2, kk)

addWorksheet(wb, "go")
writeData(wb, 3, go)

addWorksheet(wb, "do")
writeData(wb, 4, do)

saveWorkbook(wb, file = paste(root.dir, "cluster_gene_module/results.xlsx", sep = "/"), overwrite = T)
