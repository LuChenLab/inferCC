library(openxlsx)
library(ggplot2)
library(Seurat)
library(cowplot)
library(stringr)
library(dplyr)
library(reshape2)
library(tidyr)
library(scatterpie)
library(ggrepel)
library(ggrastr)
library(wesanderson)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)


args = commandArgs(trailingOnly=TRUE)

obj = readRDS(args[1])

output_dir = args[2]

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
}


patient_colors = gg_color_hue(35)
names(patient_colors) <- c(
    sapply(1:23, function(x){sprintf("PA%02d", x)}),
    sapply(1:12, function(x){sprintf("PS%02d", x)})
)

make_plot <- function(
    obj, 
    pt.size = 0.1, 
    alpha = 0.5,
    reduction.use = "tSNE",
    group.by = "batch",
    title = "",
    legend.position = "bottom",
    patients = FALSE
) {
    data = obj@meta.data
    
    data$x_c = obj@dr[[reduction.use]]@cell.embeddings[rownames(data), 1]
    data$y_c = obj@dr[[reduction.use]]@cell.embeddings[rownames(data), 2]
    
    data$group = data[, group.by]
    
    data = rbind(
        data[is.na(data$group), ],
        data[!is.na(data$group), ]
    )
    
    p <- ggplot(data = data, aes(x=x_c, y=y_c, color = group)) +
        geom_point_rast(size = pt.size, alpha = alpha) +
        theme_bw() +
        labs(
            title = str_to_title(title, locale = "en"), 
            x=paste(reduction.use, "1", sep = "_"),
            y=paste(reduction.use, "2", sep = "_"),
            color = ""
        ) +
        theme(
            plot.title = element_text(hjust = 0.5, face="bold", size = 20), 
            legend.position = legend.position,
            axis.text = element_text(size = 12),
            axis.title = element_text(size = 15)
        ) +
        guides(color = guide_legend(override.aes = list(size = 5)))
    
    if(patients) {
        p = p + scale_color_manual(
            values = c(
                patient_colors,
                c(
                    "LUAD" = "#0084D1", 
                    "LUAD_Normal"="#FECC1B", 
                    "Normal"="#73BDFF", 
                    "LUSC_Normal"="#DA6906", 
                    "LUSC"="#A0C807", 
                    "LELC"="#6A0019", 
                    "CPI"="#C5000B",
                    "Unknown"="lightgrey"
                ),
                c(
                    "I"="#EC7A21", 
                    "II"="#D73F47", 
                    "III"="#65A9A3", 
                    "IV"="#4A933E"
                )
            )
        ) 
    }

    p
}



p <- make_plot(obj, reduction.use = "tsne", group.by = "Disease")

ggsave(
    paste(output_dir, "disease_tsne.pdf", sep = "/"),
    plot = p,
    width = 6,
    height = 6
)


p <- make_plot(obj, reduction.use = "tsne", group.by = "Stage")
ggsave(
    paste(output_dir, "stage_tsne.pdf", sep = "/"),
    plot = p,
    width = 6,
    height = 6
)

luad_cells = rownames(obj@meta.data)[obj@meta.data$Disease == "LUAD"]
lusc_cells = rownames(obj@meta.data)[obj@meta.data$Disease == "LUSC"]

if(length(luad_cells) < 3 || length(lusc_cells) < 3) {
    q("no")
}


do_kegg <- function(eg, cluster=NA, pvalueCutoff = 0.05) {
    kk <- enrichKEGG(gene     = eg$ENTREZID,
                     organism     = 'hsa',
                     pvalueCutoff = pvalueCutoff)
    kk <- as.data.frame(kk)
    
    if (nrow(kk) == 0) {
        return(NULL)
    }
    # kk$cluster <- cluster
    return(kk)
}


do_go <- function(eg, cluster = NA, pvalueCutoff = 0.01, qvalueCutoff = 0.05, cutoff=0.7) {
    
    res = NULL
    for(i in c("BP", "CC", "MF")) {
        ego <- enrichGO(gene      = eg$ENTREZID,
                        keyType       = 'ENTREZID',
                        OrgDb         = org.Hs.eg.db,
                        ont           = i,
                        pAdjustMethod = "BH",
                        pvalueCutoff  = pvalueCutoff,
                        qvalueCutoff  = qvalueCutoff,
                        readable      = FALSE)
        
        if(is.null(ego)) {
            return(ego)
        }
        ego <- clusterProfiler::simplify(ego, cutoff=cutoff, by="p.adjust", select_fun=min)
        
        ego = as.data.frame(ego)
        
        if (nrow(ego) > 0) {
            ego$ONTOLOGY = i
        }
        
        if (is.null(res)) {
            res = ego
        } else {
            res = rbind(res, ego)
        }
    }
    
    if (nrow(res) == 0) {
        return(NULL)
    }
    
    # res$cluster <- cluster
    
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
        return(NULL)
    }
    # do$cluster = cluster
    
    return(do)
} 


get_entrzid <- function(markers) {
    res = NULL

    temp <- markers[markers$p_val_adj < 0.05 & abs(markers$avg_logFC) > 0.5, ]
    eg <- bitr(unique(temp$gene), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
    
    if(!is.null(eg) && nrow(eg) > 0) {
        if(is.null(res)) {
            res = eg
        } else {
            res = rbind(res, eg)
        }
    }
    
    return(res)
}

output_xlsx = paste(output_dir, "markers.xlsx", sep = "/")

if (file.exists(output_xlsx)) {
    markers = read.xlsx(output_xlsx, sheet = 1)
    kegg = read.xlsx(output_xlsx, sheet = 2)
    go = read.xlsx(output_xlsx, sheet = 3)
    do = read.xlsx(output_xlsx, sheet = 4)

    genes.label = c()
    temp_markers = markers[markers$p_val_adj < 0.05, ]

    temp_markers = temp_markers[order(temp_markers$avg_logFC, decreasing = T), ]
    genes.label = temp_markers$gene[1:10]
    temp_markers = temp_markers[order(temp_markers$avg_logFC, decreasing = F), ]
    genes.label = c(genes.label, temp_markers$gene[1:10])
} else {
    markers = FindMarkers(obj, ident.1 = luad_cells, ident.2 = lusc_cells, logfc.threshold = 0)
    markers$gene = rownames(markers)

    genes.label = c()
    temp_markers = markers[markers$p_val_adj < 0.05, ]

    temp_markers = temp_markers[order(temp_markers$avg_logFC, decreasing = T), ]
    genes.label = temp_markers$gene[1:10]
    temp_markers = temp_markers[order(temp_markers$avg_logFC, decreasing = F), ]
    genes.label = c(genes.label, temp_markers$gene[1:10])

    egs = get_entrzid(markers)

    kegg = do_kegg(egs)
    do = do_do(egs)
    go = do_go(egs)


    wb = createWorkbook()
    addWorksheet(wb, "marker genes")
    writeData(wb, 1, markers, rowNames = TRUE)

    addWorksheet(wb, "KEGG")
    writeData(wb, 2, kegg)

    addWorksheet(wb, "GO")
    writeData(wb, 3, go)

    addWorksheet(wb, "DOSE")
    writeData(wb, 4, do)

    saveWorkbook(wb, file = output_xlsx, overwrite = T)
}




make_voca_plot <- function(data, genes.label, genes.use = NULL, title = "") {

    if (is.null(genes.use)) {
        genes.use = unique(data$gene)
    } 

    data = data[data$gene %in% genes.use, ]

    data$ident[
        data$avg_logFC > 0.5 & data$p_val_adj < 0.05
    ] = "LUAD"

    data$ident[
        data$avg_logFC < -0.5 & data$p_val_adj < 0.05
    ] = "LUSC"

    data$ident[
        abs(data$avg_logFC) <= 0.5 | data$p_val_adj >= 0.05
    ] = "Not" 

    # data$ident <- apply(data, 1, function(row) {
    #     if (row[2] > 0.25 && row[5] < 0.05) {
    #         return("1")
    #     } else if (row[2] < -0.25 && row[5] < 0.05) {
    #         return("2")
    #     } 
        
    #     return("3")
    # })

    data$p_val_adj <- -log10(data$p_val_adj)
    # data$gene <- rownames(data)
    data$ident = as.factor(data$ident)

    p <- ggplot(data = data, aes(x=avg_logFC, y=p_val_adj, color = ident)) + 
        geom_point_rast() + 
        scale_color_manual(
            values=c(
                "LUAD" = "#0084D1", 
                "LUSC"="#A0C807",
                "Not"="grey50"
            ),
            breaks = c(
                "LUAD",
                "LUSC"
            )
        ) +
        geom_text_repel(
            aes(label=ifelse(gene %in% genes.label & ident != "Not", as.character(gene),''))
        ) + theme_bw() +
        theme(
            legend.position = c(0.1, 0.95),
            axis.text = element_text(size = 20),
            axis.title = element_text(size = 20),
            legend.text = element_text(size = 15),
            legend.background = element_blank(),
            legend.key = element_blank(),
            title = element_text(size = 20)
        ) +
        labs(x="log2(Fold Change)", y="-log10(FDR)", title = title, color="") +
        guides(color = guide_legend(override.aes = list(size = 5)))

    return(p)
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


print(head(markers))
print(head(kegg))
print(head(go))
print(head(do))

p <- make_voca_plot(markers, genes.label = genes.label)
ggsave(
    paste(output_dir, "volcano.pdf", sep = "/"),
    plot = p,
    width = 6,
    height = 6
)

if (!is.null(kegg)) {
    p <- make_clusterprofiler_plot(kegg, title = "KEGG")
    ggsave(
        paste(output_dir, "kegg.pdf", sep = "/"),
        plot = p,
        width = 10,
        height = 8
    )
}


if (!is.null(go)) {
    p <- make_clusterprofiler_plot(go, title = "GO", group.by = "ONTOLOGY")
    ggsave(
        paste(output_dir, "go.pdf", sep = "/"),
        plot = p,
        width = 20,
        height = 10
    )
}


if (!is.null(do)) {
    p <- make_clusterprofiler_plot(do, title = "DO")
    ggsave(
        paste(output_dir, "do.pdf", sep = "/"),
        plot = p,
        width = 10,
        height = 8
    )
}
