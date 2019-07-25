#!/usr/bin/env Rscript
# -*- coding:utf-8 -*-

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
        geom_text_repel(
            aes(label=ifelse(gene %in% genes.label, as.character(gene),''))
        ) +
        theme(legend.position = "none") +
        labs(
            x = paste0("-log2(", paste0("Cluster", cluster), "); n=", sum(meta[, group.by] == cluster)),
            y = paste0("-log2(Rest); n=", sum(meta[, group.by] != cluster))
        )
    return(p)
    
}



find_markers <- function(object, group.by = "res.0.8") {
    
    groupsby <- sort(unique(object@meta.data[, group.by]))
    res = NULL
    for(i in groupsby) {
        groups = rownames(object@meta.data[object@meta.data[, group.by] == i, ])
        
        if(length(groups) >= 2) {
            temp <- FindMarkers(
                object = object, 
                ident.1 = groups,
                group.by = "res.0.8"
            )
            
            if (!is.null(temp) && nrow(temp) > 0) {
                temp$ident = i
                temp$gene = rownames(temp)
                res = rbind(res, temp)
            }
        }
    }
    
    return(res)
}


make_clusterprofiler_plot <- function(data, group.by = NULL, topn = 20, p.val = 0.05, title = NULL) {
    data = data[data$p.adjust < 0.05, ]
    
    if (!is.null(group.by)) {
        res = NULL
        for (i in unique(data[, group.by])) {
            temp = data[data[, group.by] == i, ]
            temp = as.data.frame(temp %>% top_n(-topn, p.adjust))
            
            res = rbind(res, temp)
        }
    } else {
        res = as.data.frame(data %>% top_n(-topn, p.adjust))
    }
    
    res$GeneRatio <- sapply(res$GeneRatio, function(x) {
        temp = str_split(x, "/")[[1]]
        return(as.numeric(temp[1]) / as.numeric(temp[2]))
    })
    
    p <- ggplot(data = res, aes(x=GeneRatio, y = reorder(Description, -p.adjust), color = p.adjust, size = Count)) + geom_point()
    
    if (!is.null(group.by)) {
        p <- eval(
            parse(
                text = paste0("p + facet_grid(.~", group.by, ", scales = \"free_y\")")
            )
        )
    }
    
    p = p + 
        scale_color_gradientn(colors=rev(wes_palette("Zissou1", 100, type = "continuous"))) +
        labs(title = title, y = "")
    return(p)
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
    stage_colors = c("I"="#EC7A21", "II"="#D73F47", "III"="#65A9A3", "IV"="#4A933E", "LUAD_Normal"="#FECC1B", "LUSC_Normal"="#DA6906")
    disease_colors = c("LUAD" = "#0084D1", "LUAD_Normal"="#FECC1B", "Normal"="#73BDFF", "LUSC_Normal"="#DA6906", "LUSC"="#A0C807", "LELC"="#6A0019", "CPI"="#C5000B")
    tissue_colors = c("Bronichi"="#FB290F", "Lymph node"="#488F17", "Tumor"="#FC8210")
    patient_colors = gg_color_hue(35)
    names(patient_colors) <- c(
        sapply(1:23, function(x){sprintf("PA%02d", x)}),
        sapply(1:12, function(x){sprintf("PS%02d", x)})
    )

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


args = commandArgs(trailingOnly = T)

input_rds = args[1]
output_dir = args[2]

obj <- readRDS(input_rds)

if (length(unique(obj@meta.data$res.0.8)) < 2) {
    quit("no")
}

dir.create(output_dir, showWarnings = F, recursive = T)

output_xlsx = paste(dirname(input_rds), "markers_by_cluster.xlsx", sep = "/")

if(file.exists(output_xlsx)) {
    print("Reading results instead of calculating")
    markers = read.xlsx(output_xlsx, sheet = 1, rowNames = T)
} else {
    markers = find_markers(obj)

    write.xlsx(markers, file = output_xlsx, rowNames = T)
}


nms = rownames(obj@raw.data)
ig_genes = c(grep("^IGJ", nms, v=T), 
                grep("^IGH",nms,v=T),
                grep("^IGK", nms, v=T), 
                grep("^IGL", nms, v=T))

bad_genes = unique(c(grep("^MT-", nms, v=T), grep("^MTMR", nms, v=T), grep("^MTND", nms, v=T),"NEAT1","TMSB4X", "TMSB10", ig_genes))

markers = markers[!rownames(markers) %in% bad_genes,]
markers_de <- markers[markers$p_val_adj < 0.05 & markers$avg_logFC > 0.5, ]


genes.use = rownames(obj@raw.data)[!rownames(obj@raw.data) %in% bad_genes]


kegg = NULL
go = NULL
do = NULL
for (i in unique(markers$ident)) {
    
    temp_genes_use = markers[markers$ident == i & markers$gene %in% genes.use & markers$p_val_adj < 0.05 & abs(markers$avg_logFC) > 0.5, ]

    genes.label = c()
    temp_genes_use <- temp_genes_use[order(temp_genes_use$avg_logFC, decreasing = T), ]
    genes.label = c(temp_genes_use$gene[1:10])

    temp_genes_use <- temp_genes_use[order(temp_genes_use$avg_logFC, decreasing = F), ]
    genes.label = c(genes.label, temp_genes_use$gene[1:10])

    # p <- make_compare_dotplot_between_cluster(
    #     object = obj, 
    #     scale = F, 
    #     genes.label = genes.label,
    #     genes.use = genes.use,
    #     cluster = i
    # )

    # ggsave(
    #     filename = paste(output_dir, paste0(i, "_gene.pdf"), sep = "/"),
    #     plot = p,
    #     width = 6,
    #     height = 6,
    #     dpi = 600,
    #     units = "in"
    # ) 


    p <- make_voca_plot(markers[markers$ident == i & markers$gene %in% genes.use, ], genes.label = genes.label, genes.use = genes.use)

    ggsave(
        filename = paste(output_dir, paste0(i, "_volca.pdf"), sep = "/"),
        plot = p,
        width = 6,
        height = 6,
        dpi = 600,
        units = "in"
    )


    ## cluster modules
    eg = get_entrzid(temp_genes_use)

    temp = do_kegg(eg, cluster = i)
    kegg = rbind(kegg, temp)

    if(!is.null(temp) && nrow(temp) > 0) {
        p = make_clusterprofiler_plot(temp, title = paste0("KEGG (", i, ")"))
        ggsave(
            filename = paste(output_dir, paste0(i, "_kegg.pdf"), sep = "/"),
            plot = p,
            width = 12,
            height = 12,
            dpi = 600,
            units = "in"
        )
    }


    temp = do_go(eg, cluster = i)
    go = rbind(go, temp)

    if(!is.null(temp) && nrow(temp) > 0) {
        p = make_clusterprofiler_plot(temp, group.by="ONTOLOGY", title = paste0("GO (", i, ")"))
        ggsave(
            filename = paste(output_dir, paste0(i, "_go.pdf"), sep = "/"),
            plot = p,
            width = 16,
            height = 12,
            dpi = 600,
            units = "in"
        )
    }

    temp = do_do(eg, cluster = i)
    do = rbind(do, temp)

    if(!is.null(temp) && nrow(temp) > 0) {
        p = make_clusterprofiler_plot(temp, title = paste0("DO (", i, ")"))
        ggsave(
            filename = paste(output_dir, paste0(i, "_do.pdf"), sep = "/"),
            plot = p,
            width = 12,
            height = 12,
            dpi = 600,
            units = "in"
        )
    }
}


wb = createWorkbook()
addWorksheet(wb, "KEGG")
writeData(wb, 1, kegg)

addWorksheet(wb, "GO")
writeData(wb, 2, go)

addWorksheet(wb, "DO")
writeData(wb, 3, do)

saveWorkbook(wb, file = paste(output_dir, "clusterprofiler.xlsx", sep = "/"), overwrite = T)
