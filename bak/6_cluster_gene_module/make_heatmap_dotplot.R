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
    library("wesanderson")
    library(RColorBrewer)
    library(ComplexHeatmap)
}

suppressPackageStartupMessages(load())


args = commandArgs(trailingOnly = TRUE)

root.dir = args[1]
rds = args[2]
topN = as.numeric(args[3])
group.by = args[4]


setwd(root.dir)

dir.create(paste(root.dir, "cluster_gene_module", sep = "/"), showWarnings = F)

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


gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
}


make_complex_heatmap <- function(obj, module, module_id=NULL, group.by = "res.0.6", order.by = "Cluster") {

    if (!is.null(module_id)) {
        module = module[module$ident == module_id, ]
    }

    rownames(module) <- module$gene
    temp_module <- module[, 1, drop=FALSE]
    temp_module[,1] <- as.character(temp_module[, 1])

    # print(colnames(obj@meta.data))
    meta <- obj@meta.data[, c("Stage", "PatientID", group.by)]
    colnames(meta) <- c("Stage", "Patient", "Cluster")
    meta <- meta[order(meta[, order.by]), ]

    ### extract temp information
    stage_colors = c("Benign"="#3E6596", "I"="#EC7A21", "II"="#D73F47", "III"="#65A9A3", "IV"="#4A933E")
    disease_colors = c("ADC" = "#063373", "LCC"="#FECC1B", "Normal"="#73BDFF", "Other"="#DA6906", "SCC"="#A0C807", "LELC"="#6A0019", "CPI"="#C5000B", "LBT"="#0084D1")
    tissue_colors = c("Bronichi"="#FB290F", "Lymph node"="#488F17", "Tumor"="#FC8210")
    patient_colors = gg_color_hue(50)

    mat = MinMax(as.matrix(obj@scale.data)[rownames(module), rownames(meta)], -2.5, 2.5)

    cols_stage = stage_colors[unique(meta$Stage)]
    # print(cols_stage)

    cols_patients = patient_colors[as.numeric(gsub("P", "", sort(unique(meta$Patient))))]
    names(cols_patients) <- sort(unique(meta$Patient))
    # print(cols_patients)

    cols_cluster = wes_palette("Darjeeling2", length(unique(meta$Cluster)), type = "continuous")
    names(cols_cluster) <- unique(meta$Cluster)
    # print(cols_cluster)

    ann_colors = list(
        Stage = cols_stage,
        Patient = cols_patients,
        Cluster = cols_cluster
    )


    ha = HeatmapAnnotation(
        df = meta,
        col = ann_colors
    )

    ht = Heatmap(
        mat,
        name = "mat",
        cluster_rows = F,
        cluster_columns = F,
        show_column_names = F,
        top_annotation = ha,
        column_title = module_id
    )
    draw(ht)
}



# function to make heatmaps
# :param obj Seurat obj
# :param cluster results from perform_mfuzz
# :param output_prefix: the output file prefix
make_heatmap <- function(obj, cluster, output_prefix=NULL, group.by = "res.0.6") {


    for (i in sort(unique(cluster$ident))) {
        temp_genes = cluster[cluster$ident == i, "gene"]

        height = length(temp_genes) / 6

        if (height > 45) {
            height = 45
        } else if (height < 5) {
            height = 5
        }

        if (is.null(output_prefix)) {
            make_complex_heatmap(obj, cluster, i, group.by = group.by)
        } else {
            # ggsave(
            #     paste(output_prefix, i, ".png", sep = ""),
            #     p,
            #     width = 12,
            #     height = height,
            #     dpi = 300,
            #     units = "in"
            # )

            png(paste(output_prefix, i, ".png", sep = ""), width = 12, height = height, res = 600, units = "in")
            make_complex_heatmap(obj, cluster, i, group.by = group.by)
            dev.off()
        }
    }
}




# function to make dotplot
# :param obj Seurat obj
# :param cluster results from perform_mfuzz
# :param output_prefix: the output file prefix
make_dotplot <- function(obj, cluster, output_prefix=NULL, group.by="res.0.6") {
    for (i in sort(unique(cluster$ident))) {

        temp_genes = cluster[cluster$ident == i, "gene"]

        p <- DotPlot(
            obj, 
            temp_genes,
            group.by = group.by, 
            x.lab.rot = TRUE, 
            do.return = T
        ) + labs(title = i) + coord_flip()
        
        height = length(temp_genes) / 6
        
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
    stages = sort(unique(obj@meta.data[group.by]))
    
    res = as.data.frame(matrix(0, nrow = length(groups), ncol = length(stages)))
    rownames(res) = groups
    colnames(res) = stages
    
    for (i in groups) {
        temp_gene = names(cluster$cluster[cluster$cluster == i])
        
        for (j in stages) {
            # print(paste0(i, j))
            temp_cells = rownames(obj@meta.data[obj@meta.data[group.by] == j, ])
            
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
    for(i in unique(markers$ident)) {
        print(i)
        
        temp <- markers[markers$ident == i, "gene"]
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


# expr = read_sctransform()
obj <- readRDS(rds)

cluster_markers = "annotation_results_by_cluster.xlsx"
cluster_markers = read.xlsx(paste(root.dir, cluster_markers, sep = "/"), rowNames = T)


print(unique(cluster_markers$ident))

genes = as.data.frame(cluster_markers %>% 
                          filter(p_val_adj < 0.05 & avg_logFC > 0) %>% 
                          group_by(ident) %>% 
                          top_n(topN, wt=avg_logFC))

print("heatmap")
make_heatmap(
    obj, 
    genes,
    group.by = group.by,
    output_prefix = paste(root.dir, "cluster_gene_module/heatmap_", sep = "/")
    )

print("dotplot")
make_dotplot(
    obj, 
    genes,
    group.by = group.by,
    output_prefix = paste(root.dir, "cluster_gene_module/dotplot_", sep = "/")
)



## KEGG, GO and DOSE
#
# eg = get_entrzid(genes)
#
# kk <- NULL
#
# for(i in unique(eg$ident)) {
#     print(i)
#     temp = do_kegg(eg[eg$ident == i, ])
#
#     if(is.null(temp) || nrow(temp) == 0) {
#         next
#     }
#
#     temp$ident = i
#     if(is.null(kk)) {
#         kk = temp
#     } else {
#         kk = rbind(kk, temp)
#     }
# }
# go <- NULL
#
# for(i in unique(eg$ident)) {
#     print(i)
#     temp = do_go(eg[eg$ident == i, ])
#
#     if(is.null(temp) || nrow(temp) == 0) {
#         next
#     }
#
#     temp$ident = i
#     if(is.null(go)) {
#         go = temp
#     } else {
#         go = rbind(go, temp)
#     }
# }
# do <- NULL
#
# for(i in unique(eg$ident)) {
#     print(i)
#     temp = do_do(eg[eg$ident == i, ])
#
#     if(is.null(temp) || nrow(temp) == 0) {
#         next
#     }
#
#     temp$ident = i
#     if(is.null(do)) {
#         do = temp
#     } else {
#         do = rbind(do, temp)
#     }
# }
#
#
# wb = createWorkbook()
# addWorksheet(wb, "gene_module")
# writeData(wb, 1, genes)
#
# addWorksheet(wb, "kegg")
# writeData(wb, 2, kk)
#
# addWorksheet(wb, "go")
# writeData(wb, 3, go)
#
# addWorksheet(wb, "do")
# writeData(wb, 4, do)
#
# saveWorkbook(wb, file = paste(root.dir, "cluster_gene_module/results.xlsx", sep = "/"), overwrite = T)
