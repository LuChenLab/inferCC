library(openxlsx)
library(ggplot2)
library(cowplot)
library(reshape2)
library(stringr)
library(clusterProfiler)
library(DOSE)
library(tmod)
library(org.Hs.eg.db)
library(dplyr)
library(wesanderson)
library(Seurat)
library(ggrepel)


stage_colors = c("Benign"="#3E6596", "I"="#EC7A21", "II"="#D73F47", "III"="#65A9A3", "IV"="#4A933E")
disease_colors = c("ADC" = "#063373", "LCC"="#FECC1B", "Normal"="#73BDFF", "Other"="#DA6906", "SCC"="#A0C807", "LELC"="#6A0019", "CPI"="#C5000B", "LBT"="#0084D1")
tissue_colors = c("Bronichi"="#FB290F", "Lymph node"="#488F17", "Tumor"="#FC8210")



make_tmod_cp_style_dotplot <- function(resC, resU, resZ, topn=10, auc = 0.5, ident = NULL) {
    data = NULL
    
    if (!is.null(ident)) {
        resC = resC[resC$ident == ident, ]
        resU = resU[resU$ident == ident, ]
        resZ = resZ[resZ$ident == ident, ]
    }
    
    resC = resC[resC$AUC > auc, ]
    resU = resU[resU$AUC > auc, ]
    resZ = resZ[resZ$AUC > auc, ]
    
    #     print(resC %>% top_n(10, wt=P.Value))
    # temp = as.data.frame(resC %>% arrange(P.Value) %>% select(Title, N1, AUC, P.Value) %>% head(topn))
    temp = as.data.frame(resC %>% top_n(-topn, P.Value) %>% dplyr::select(Title, N1, AUC, P.Value))
    
    if(nrow(temp) > 0) {
        # print(temp)
        temp$ident = "CERNO"
        data = rbind(data, temp)
        
    }
    
    temp = NULL
    # temp = as.data.frame(resU %>% arrange(P.Value) %>% select(Title, N1, AUC, P.Value) %>% head(topn))
    temp = as.data.frame(resU %>% top_n(-topn, P.Value) %>% dplyr::select(Title, N1, AUC, P.Value))
    # print(temp)
    if(!is.null(temp) && nrow(temp) > 0) {
        # print(temp)
        temp$ident = "Utest"
        data = rbind(data, temp)
        
    }
    
    temp = NULL
    # temp = as.data.frame(resZ %>% arrange(P.Value) %>% select(Title, N1, AUC, P.Value) %>% head(topn))
    temp = as.data.frame(resZ %>% top_n(-topn, P.Value) %>% dplyr::select(Title, N1, AUC, P.Value))
    # print(temp)
    if(!is.null(temp) && nrow(temp) > 0) {
        # print(temp)
        temp$ident = "Ztest"
        data = rbind(data, temp)
        
    }
    data = na.omit(data)
    
    p <- ggplot(
        data = data,
        aes(x = AUC, y=reorder(Title, AUC), color = P.Value, size = N1)) +
        geom_point() +
        facet_grid(.~ident, scales = "free") +
        scale_color_gradientn(colors=rev(wes_palette("Zissou1", 100, type = "continuous"))) +
        labs(y="")
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


make_compare_dotplot_between_disease <- function(object, genes.label, scale=FALSE, disease.1 = "ADC", disease.2 = "SCC") {
    meta = object@meta.data[object@meta.data$Disease %in% c(disease.1, disease.2), ]
    
    if(scale) {
        data = as.matrix(object@scale.data[ ,rownames(meta)])
    } else {
        data = as.matrix(object@raw.data[, rownames(meta)])
    }
    
    temp = as.data.frame(matrix(NA, nrow = nrow(data), ncol = 3))
    
    colnames(temp) <- c("gene", disease.1, disease.2)
    
    temp[, 1] = rownames(data)
    temp[, disease.1] = apply(data[, rownames(meta[meta$Disease == disease.1,])], 1, function(x) {
        return(mean(as.numeric(x)))
    })
    temp[, disease.2] = apply(data[, rownames(meta[meta$Disease == disease.2,])], 1, function(x) {
        return(mean(as.numeric(x)))
    })
    
    temp[,2] = log2(temp[, 2] + 1)
    temp[,3] = log2(temp[, 3] + 1)
    
    max_ = ceiling(max(temp[, 2:3])) + 1
    min_ = floor(min(temp[, 2:3]))
    
    temp$gene[!temp$gene %in% genes.label] = ""
    
    p <- eval(
        parse(text = paste("ggplot(data = temp, aes(x =", disease.1, ", y =", disease.2, ", label=gene))"))   
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
            x = paste0("-log2(", disease.1, ")"),
            y = paste0("-log2(", disease.2, ")")
        )
    return(p)
    
}



args = commandArgs(trailingOnly = TRUE)


output = args[2]

dir.create(output, showWarnings = F, recursive = T)

obj <- readRDS(args[1])


markers <- FindMarkers(
    obj, 
    ident.1 = rownames(obj@meta.data[obj@meta.data$Disease == "SCC", ]), 
    ident.2 = rownames(obj@meta.data[obj@meta.data$Disease == "ADC", ])
)
markers$gene = rownames(markers)

markers_de <- markers[markers$p_val_adj < 0.05, ]

write.xlsx(markers, paste(output, "DEGs_ADC_SCC.xlsx", sep = "/"))



genes.label = c()
markers_de <- markers_de[order(markers_de$avg_logFC, decreasing = T), ]
genes.label = c(markers_de$gene[1:10])

markers_de <- markers_de[order(markers_de$avg_logFC, decreasing = F), ]
genes.label = c(genes.label, markers_de$gene[1:10])


p <- make_compare_dotplot_between_disease(object = obj, scale = F, genes.label = genes.label)

ggsave(filename = paste(output, "most_DEGs_scatter.png", sep = "/"), plot = p, width = 6, height = 6, units = "in", dpi = 600)



eg <- bitr(markers_de$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)


kk <- enrichKEGG(
    eg$ENTREZID,
    organism = "hsa", 
    keyType = "kegg",
    pvalueCutoff = 0.05, 
    pAdjustMethod = "BH", 
    minGSSize = 10, 
    maxGSSize = 500, 
    qvalueCutoff = 0.2,
    use_internal_data = FALSE
)

kk = as.data.frame(kk)



do <- enrichDO(gene         = eg$ENTREZID,
               ont           = "DO",
               pvalueCutoff  = 0.05,
               pAdjustMethod = "BH",
               minGSSize     = 5,
               maxGSSize     = 500,
               qvalueCutoff  = 0.05,
               readable      = FALSE)

do <- as.data.frame(do)


ego <- enrichGO(gene          = eg$ENTREZID,
                OrgDb         = org.Hs.eg.db,
                ont           = "All",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)

ego <- as.data.frame(ego)


msig <- tmodImportMSigDB("/mnt/raid62/Lung_cancer_10x/00_data_ingest/00_raw_data/msigdb_v6.2.xml")
sel <- msig$MODULES$Category %in% c("H", "C5", "C2")


resU = tmodUtest(
    markers$gene[order(abs(markers$avg_logFC), decreasing = T)],
    mset = msig[sel],
    qval=0.05
)

resZ = tmodZtest(
    markers$gene[order(abs(markers$avg_logFC), decreasing = T)],
    mset = msig[sel],
    qval=0.05
)

resC = tmodCERNOtest(
    markers$gene[order(abs(markers$avg_logFC), decreasing = T)],
    mset = msig[sel],
    qval=0.05
)



wb = createWorkbook()
addWorksheet(wb, "KEGG")
writeData(wb, 1, kk)

addWorksheet(wb, "DO")
writeData(wb, 2, do)

addWorksheet(wb, "GO")
writeData(wb, 3, ego)

addWorksheet(wb, "Utest")
writeData(wb, 4, resU)

addWorksheet(wb, "Ztest")
writeData(wb, 5, resZ)

addWorksheet(wb, "CERNO")
writeData(wb, 6, resC)

saveWorkbook(wb, paste(output, "annotation.xlsx", sep = "/"), overwrite = T)



p <- make_tmod_cp_style_dotplot(resC = resC, resZ = resZ, resU = resU, topn = 20)

ggsave(
    filename = paste(output, "tmod.png", sep = "/"),
    plot = p,
    width = 12,
    height = 12,
    dpi = 600,
    units = "in"
)


p <- make_clusterprofiler_plot(ego, group.by = "ONTOLOGY", title = "GO")
ggsave(
    filename = paste(output, "go.png", sep = "/"),
    plot = p,
    width = 16,
    height = 12,
    dpi = 600,
    units = "in"
)


p <- make_clusterprofiler_plot(kk, topn = 40, title = "KEGG")
ggsave(
    filename = paste(output, "kegg.png", sep = "/"),
    plot = p,
    width = 12,
    height = 12,
    dpi = 600,
    units = "in"
)


p <- make_clusterprofiler_plot(do, topn = 40, title = "DO")
ggsave(
    filename = paste(output, "do.png", sep = "/"),
    plot = p,
    width = 12,
    height = 12,
    dpi = 600,
    units = "in"
)

