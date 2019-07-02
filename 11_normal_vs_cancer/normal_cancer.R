#!/usr/bin/env Rscript

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


stage_colors = c("Benign"="#3E6596", "I"="#EC7A21", "II"="#D73F47", "III"="#65A9A3", "IV"="#4A933E")
disease_colors = c("ADC" = "#063373", "LCC"="#FECC1B", "Normal"="#73BDFF", "Other"="#DA6906", "SCC"="#A0C807", "LELC"="#6A0019", "CPI"="#C5000B", "LBT"="#0084D1")
tissue_colors = c("Bronichi"="#FB290F", "Lymph node"="#488F17", "Tumor"="#FC8210")

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


make_compare_dotplot_between_disease <- function(object, genes.label, scale=FALSE, group.by = "Disease", disease.1 = "ADC", disease.2 = "SCC", genes.use = NULL) {
    meta = object@meta.data

    if (is.null(genes.use)) {
        genes.use = rownames(object@raw.data)
    }
    
    if(scale) {
        data = as.matrix(object@scale.data[genes.use ,rownames(meta)[object@meta.data[, group.by] %in% c(disease.1, disease.2)]])
    } else {
        data = as.matrix(object@raw.data[genes.use, rownames(meta)[object@meta.data[, group.by] %in% c(disease.1, disease.2)]])
    }
    
    if (is.null(disease.2)) {
        disease.2 = "Rest"
    }
    
    temp = as.data.frame(matrix(NA, nrow = nrow(data), ncol = 3))
    
    colnames(temp) <- c("gene", as.character(disease.1), disease.2)
    
    temp[, 1] = rownames(data)
    temp[, as.character(disease.1)] = apply(data[, rownames(meta[meta[, group.by] == disease.1,])], 1, function(x) {
        return(mean(as.numeric(x)))
    })
    
    if(disease.2 != "Rest") {
        temp[, disease.2] = apply(data[, rownames(meta[meta[, group.by] == disease.2,])], 1, function(x) {
            return(mean(as.numeric(x)))
        })
    } else {
        temp[, disease.2] = apply(data[, rownames(meta[meta[, group.by] != disease.1,])], 1, function(x) {
            return(mean(as.numeric(x)))
        })
    }
    
    temp[,2] = log2(temp[, 2] + 1)
    temp[,3] = log2(temp[, 3] + 1)
    
    max_ = ceiling(max(temp[, 2:3]))
    min_ = floor(min(temp[, 2:3]))
    
    p <- eval(
        parse(text = paste("ggplot(data = temp, aes(x =", disease.1, ", y =", disease.2, "))"))   
    )
      
    p = p + 
        geom_point() +
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
        geom_text_repel(aes(label=ifelse(gene %in% genes.label, as.character(gene),''), color = "red")) +
        theme(legend.position = "none") +
        labs(
            x = paste0("-log2(", disease.1, "); n=", sum(meta[, group.by] == disease.1)),
            y = paste0("-log2(", disease.2, "); n=", sum(meta[, group.by] == disease.2))
        )
    return(p)
    
}


make_voca_plot <- function(data, genes.label, genes.use = NULL) {

    if (is.null(genes.use)) {
        genes.use = rownames(data)
    } 

    data = data[genes.use, ]
    data = na.omit(data)
    print(dim(data))

    data$ident <- apply(data, 1, function(row) {
        if (row[2] > 0.25 && row[5] < 0.05) {
            return("1")
        } else if (row[2] < -0.25 && row[5] < 0.05) {
            return("2")
        } 
        
        return("3")
    })
    
    
    data$p_val_adj <- -log10(data$p_val_adj)
    data$gene <- rownames(data)
    data$ident = as.factor(data$ident)

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


args = commandArgs(trailingOnly = T)

input_rds = args[1]
output_dir = args[2]


obj <- readRDS(input_rds)


if (!"Normal" %in% obj@meta.data$Disease) {
    quit("no")
}

dir.create(output_dir, showWarnings = F, recursive = T)
normal_cells = rownames(obj@meta.data)[obj@meta.data$Disease == "Normal"]

print(length(normal_cells))

nms = rownames(obj@raw.data)
ig_genes = c(grep("^IGJ", nms, v=T), 
                grep("^IGH",nms,v=T),
                grep("^IGK", nms, v=T), 
                grep("^IGL", nms, v=T))

bad_genes = unique(c(grep("^MT-", nms, v=T), grep("^MTMR", nms, v=T), grep("^MTND", nms, v=T),"NEAT1","TMSB4X", "TMSB10", ig_genes))

genes.use = rownames(obj@raw.data)[!rownames(obj@raw.data) %in% bad_genes]


for (i in c("ADC", "SCC")) {
    disease_cells = rownames(obj@meta.data)[obj@meta.data$Disease == i]

    if (length(normal_cells) > 10 && length(disease_cells) > 10) {

        output_xlsx = paste(output_dir, paste0(i, ".xlsx"), sep = "/")
        if(file.exists(output_xlsx)) {
            print("Reading results instead of calculating")
            markers = read.xlsx(output_xlsx, sheet = 1, rowNames = T)
            markers = markers[!rownames(markers) %in% bad_genes,]

            markers_de <- markers[markers$p_val_adj < 0.05 & abs(markers$avg_logFC) > 0.5, ]

            kk <- read.xlsx(output_xlsx, sheet = 2)
            do <- read.xlsx(output_xlsx, sheet = 3)
            ego <- read.xlsx(output_xlsx, sheet = 4)
        } else {
            markers = FindMarkers(obj, ident.1 = normal_cells, ident.2 = disease_cells, logfc.threshold = 0)
            markers = markers[!rownames(markers) %in% bad_genes,]

            markers_de <- markers[markers$p_val_adj < 0.05 & abs(markers$avg_logFC) > 0.5, ]
            eg <- bitr(rownames(markers_de), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)


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


            do <- enrichDO(gene = eg$ENTREZID,
                ont           = "DO",
                pvalueCutoff  = 0.05,
                pAdjustMethod = "BH",
                minGSSize     = 5,
                maxGSSize     = 500,
                qvalueCutoff  = 0.05,
                readable      = FALSE)

            do <- as.data.frame(do)


            ego <- enrichGO(gene  = eg$ENTREZID,
                    OrgDb         = org.Hs.eg.db,
                    ont           = "All",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05,
                    readable      = TRUE)

            ego <- as.data.frame(ego)

            wb = createWorkbook()
            addWorksheet(wb, "Markers")
            writeData(wb, 1, markers, rowNames = T)

            addWorksheet(wb, "KEGG")
            writeData(wb, 2, kk)

            addWorksheet(wb, "DO")
            writeData(wb, 3, do)

            addWorksheet(wb, "GO")
            writeData(wb, 4, ego)

            saveWorkbook(wb, paste(output_dir, paste0(i, ".xlsx"), sep = "/"), overwrite=T)
        }
        
        # print("go")

        # if(!is.null(ego) && nrow(ego) > 0) {
        #     p <- make_clusterprofiler_plot(ego, group.by = "ONTOLOGY", title = "GO")
        #     ggsave(
        #         filename = paste(output_dir, paste0(i, "_go.png"), sep = "/"),
        #         plot = p,
        #         width = 16,
        #         height = 12,
        #         dpi = 600,
        #         units = "in"
        #     )
        # }

        # print("kegg")
        # if(!is.null(kk) && nrow(kk) > 0) {
        #     p <- make_clusterprofiler_plot(kk, topn = 40, title = "KEGG")
        #     ggsave(
        #         filename = paste(output_dir, paste0(i, "_kegg.png"), sep = "/"),
        #         plot = p,
        #         width = 12,
        #         height = 12,
        #         dpi = 600,
        #         units = "in"
        #     )
        # }

        # print("do")
        # if(!is.null(do) && nrow(do) > 0) {
        #     p <- make_clusterprofiler_plot(do, topn = 40, title = "DO")
        #     ggsave(
        #         filename = paste(output_dir, paste0(i, "_do.png"), sep = "/"),
        #         plot = p,
        #         width = 12,
        #         height = 12,
        #         dpi = 600,
        #         units = "in"
        #     )
        # }

        genes.label = c()
        markers_de <- markers_de[order(markers_de$avg_logFC, decreasing = T), ]
        genes.label = c(rownames(markers_de)[1:10])

        markers_de <- markers_de[order(markers_de$avg_logFC, decreasing = F), ]
        genes.label = c(genes.label, rownames(markers_de)[1:10])


        # p <- make_compare_dotplot_between_disease(
        #     object = obj, 
        #     scale = F, 
        #     genes.label = genes.label,
        #     disease.1 = "Normal",
        #     disease.2 = i,
        #     genes.use = genes.use
        # )

        # print("scatter")
        # ggsave(
        #     filename = paste(output_dir, paste0(i, "_gene.png"), sep = "/"),
        #     plot = p,
        #     width = 6,
        #     height = 6,
        #     dpi = 600,
        #     units = "in"
        # ) 

        markers$avg_logFC <- -1 * markers$avg_logFC

        print(dim(markers))

        p <- make_voca_plot(markers, genes.label)

        ggsave(
            filename = paste(output_dir, paste0(i, "_voca.png"), sep = "/"),
            plot = p,
            width = 6,
            height = 6,
            dpi = 600,
            units = "in"
        ) 
    }
}

