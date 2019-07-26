library(bigSCale)
library(Seurat)
library(Matrix)
library(openxlsx)
library(ggplot2)
library(ggrepel)
library(ggrastr)


EXHAUST = c(
    "CD274",
    "PDCD1",
    "CTLA4",
    "TIGIT",
    "LAG3",
    "NBL1",
    "CXCL13",
    "LAG3",
    "HAVCR2",
    "PTMS",
    "FAM3C",
    "IFNG",
    "AKAP5",
    "CD7",
    "PHLDA1",
    "ENTPD1",
    "TNS3",
    "CXCL13",
    "RDH10",
    "DGKH",
    "KIR2DL4",
    "LYST",
    "MIR155HG",
    "RAB27A",
    "CSF1",
    "TNFRSF9",
    "CD27",
    "CCL3",
    "ITGAE",
    "PAG1",
    "TNFRSF1B",
    "GALNT1",
    "GBP2",
    "MYO7A"
)


# obj <- readRDS("03_each_cells/LUSC/CD4/seurat.rds")

# writeMM(obj@raw.data, "test.mtx")

# out=iCells(file.dir = "test.mtx", target.cells = 100, sample.conditions = obj@meta.data[colnames(obj@raw.data),])


args = commandArgs(trailingOnly=TRUE)

# luad = readRDS("03_each_cells/LUAD/CD8/seurat.rds")
# lusc = readRDS("03_each_cells/LUSC/CD8/seurat.rds")

luad = readRDS(args[1])
lusc = readRDS(args[2])

dir.create(args[3], showWarnings = FALSE, recursive = TRUE)

setwd(args[3])


nms = c(rownames(luad@raw.data), rownames(lusc@raw.data))
ig_genes = c(grep("^IGJ", nms, v=T), 
                grep("^IGH",nms,v=T),
                grep("^IGK", nms, v=T), 
                grep("^IGL", nms, v=T))

bad_genes = unique(c(grep("^MT-", nms, v=T), grep("^MTMR", nms, v=T), grep("^MTND", nms, v=T),"NEAT1","TMSB4X", "TMSB10", ig_genes))

gene.names = nms[!nms %in% bad_genes]

luad_network = compute.network(expr.data = luad@raw.data, gene.names = gene.names, clustering = 'direct')
lusc_network =compute.network(expr.data = lusc@raw.data, gene.names = gene.names, clustering = 'direct')

saveRDS(luad_network, "LUAD_network.rds")
saveRDS(lusc_network, "LUSC_network.rds")

comparison=compare.centrality(list(luad_network$centrality,lusc_network$centrality),c('LUAD','LUSC'))

saveRDS(comparison, "comparison.rds")


        
wb = createWorkbook()
addWorksheet(wb, "Degree")
writeData(wb, 1, comparison$Degree, rowNames = TRUE)

addWorksheet(wb, "Betweenness")
writeData(wb, 2, comparison$Betweenness, rowNames = TRUE)

addWorksheet(wb, "Closeness")
writeData(wb, 3, comparison$Closeness, rowNames = TRUE)

addWorksheet(wb, "PAGErank")
writeData(wb, 4, comparison$PAGErank, rowNames = TRUE)

saveWorkbook(wb, file = "comparison.xlsx", overwrite = T)


## make plots

for(i in names(comparison)) {
    data = comparison[[i]]
    
    colnames(data) = c("LUAD", "LUSC", "DELTA", "Ranking")

    data$gene = rownames(data)
    data$ident = "Not"
    data$ident[data$gene %in% EXHAUST] = "red"

    # data[, 1] = log10(data[, 1])
    # data[, 2] = log10(data[, 2])
    
    p <- ggplot(data = data, aes(x=LUAD, y=LUSC, color = ident)) + 
        geom_point_rast() + 
        scale_color_manual(
            values=c(
                "red" = "red", 
                "Not"="grey50"
            )
        ) +
        geom_text_repel(
            aes(label=ifelse(gene %in% EXHAUST, as.character(gene),''))
        ) + theme_bw() +
        theme(legend.position = "none") +
        labs(title=i, color = "", x="LUAD", y="LUSC") +
        geom_abline(
            color="red", 
            linetype="dashed"
        )

    ggsave(
        paste0(i, ".pdf"),
        plot = p,
        width = 6,
        height = 6
    )
}