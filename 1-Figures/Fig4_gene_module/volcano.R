options(stringsAsFactors = F)

library(Seurat)
library(ggrastr)
library(openxlsx)
library(ggrepel)
library(ggpubr)
library(ggthemes)





make_voca_plot <- function(data, genes.label, genes.use = NULL, title = "") {
    
    if (is.null(genes.use)) {
        genes.use = unique(data$gene)
    } 
    
    data = data[data$gene %in% genes.use, ]
    
    if (!"ident" %in% colnames(data)) {
        print("assign ident")
        data$ident[
            data$avg_logFC > 0.5 & data$p_val_adj < 0.05
        ] = "LUAD"
        
        data$ident[
            data$avg_logFC < -0.5 & data$p_val_adj < 0.05
        ] = "LUSC"
        
        data$ident[
            abs(data$avg_logFC) <= 0.5 | data$p_val_adj >= 0.05
        ] = "Not" 
    }
    
    
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
                "Not"="grey50",
                "ATII" = "#0084D1", 
                "Basal"="#A0C807"
            )
            # breaks = c(
            #     "LUAD",
            #     "LUSC"
            # )
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


get_gene_to_label <- function(markers, n=10) {
    genes.label = c()
    temp_markers = markers[markers$p_val_adj < 0.05, ]
    
    temp_markers = temp_markers[order(temp_markers$avg_logFC, decreasing = T), ]
    genes.label = temp_markers$gene[1:n]
    temp_markers = temp_markers[order(temp_markers$avg_logFC, decreasing = F), ]
    genes.label = c(genes.label, temp_markers$gene[1:n])
}


atii <- read.xlsx("03_each_cells/total/Alveolar_II/disease/markers.xlsx", rowNames = T)

genes.label = get_gene_to_label(atii)


p <- make_voca_plot(
    atii, 
    genes.label = c(
        genes.label, 
        c("ACSL4", "ALPL", "AQP4", "BMP2", "SFTPC", "C4BPA", "VIM")
    )
)


ggsave(
    "03_each_cells/total/Alveolar_II/disease/volcano.pdf",
    plot = p,
    width = 6,
    height = 6
)



basal <- read.xlsx("03_each_cells/total/Basal/disease/markers.xlsx", rowNames = T)

genes.label = get_gene_to_label(basal)


p <- make_voca_plot(
    basal, 
    genes.label = c(
        genes.label, 
        c("ACSL4", "ALPL", "AQP4", "BMP2", "SFTPC", "C4BPA", "VIM")
    )
)
p


ggsave(
    "03_each_cells/total/Basal/disease/volcano.pdf",
    plot = p,
    width = 6,
    height = 6
)


di <- read.csv("03_each_cells/select_features/atii_basal_markers.csv", row.names = 1, stringsAsFactors = F)

di$ident[di$avg_logFC > 0.5 & di$p_val_adj < 0.05] = "Basal"

di$ident[di$avg_logFC < -0.5 & di$p_val_adj < 0.05] = "ATII"

di$ident[abs(di$avg_logFC) <= 0.5 | di$p_val_adj >= 0.05] = "Not" 

genes.label = get_gene_to_label(di)
genes.label = c(
    genes.label, 
    c("ACSL4", "ALPL", "AQP4", "BMP2", "SFTPC", "C4BPA", "VIM")
)



p <- make_voca_plot(di, genes.label) +
    theme(legend.position = c(0.1, 0.15))
p

ggsave(
    "03_each_cells/select_features/volcano.pdf",
    plot = p,
    width = 6,
    height = 6
)