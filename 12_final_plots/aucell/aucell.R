BiocManager::install(c("doMC", "doRNG"))
# For the main example:
BiocManager::install(c("mixtools", "GEOquery", "SummarizedExperiment"))
# For the examples in the follow-up section of the tutorial:
BiocManager::install(c("DT", "plotly", "NMF", "d3heatmap", "shiny", "rbokeh"))


library(AUCell)
library(Seurat)
library(ggrastr)

# obj <- readRDS("/mnt/raid62/Lung_cancer_10x/final_plots/02_rds/seurat_obj.rds")

# expr <- read.csv("/mnt/raid62/Lung_cancer_10x/final_plots/02_rds/hx_expr.csv.gz", row.names = 1)
expr = readRDS("/mnt/raid62/Lung_cancer_10x/final_plots/02_rds/all_cell_expr.rds")

cells_rankings <- AUCell_buildRankings(as.matrix(expr))


saveRDS(cells_rankings, "/mnt/raid62/Lung_cancer_10x/final_plots/02_rds/02_aucell/cells_rankings.rds")

markers = read.xlsx("/mnt/raid62/Lung_cancer_10x/final_plots/20190701_gene_markers.xlsx", sheet = 3)
markers <- na.omit(markers[, 1:2])

geneSets = list()
markers[, 2] <- as.character(markers[, 2])
for(i in unique(markers[, 2])) {
    temp = markers[markers[,2] == i, ]

    geneSets[[i]] = temp[, 1]
}

# geneSets <- GSEABase::GeneSet(genes, setName="geneSet1") # alternative
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)

pdf(
    "/mnt/raid62/Lung_cancer_10x/final_plots/01_first_plots/aucell/geneSets_threshold.pdf", 
    width = 18, 
    height = 6 * ceiling(length(geneSets) / 3),
)
par(mfrow=c(ceiling(length(geneSets) / 3), 3))
set.seed(123)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
dev.off()


saveRDS(cells_AUC, "/mnt/raid62/Lung_cancer_10x/final_plots/02_rds/02_aucell/cells_AUC.rds")
saveRDS(cells_assignment, "/mnt/raid62/Lung_cancer_10x/final_plots/02_rds/02_aucell/cells_assignment.rds")

selectedThresholds <- getThresholdSelected(cells_assignment)



## make AUCell plots

for(i in names(selectedThresholds)) {
    
    pdf(
        paste(
            "/mnt/raid62/Lung_cancer_10x/final_plots/01_first_plots/aucell", 
            paste0(str_replace(i, "\\s", "_"), "_aucell.pdf"), 
            sep = "/"
        ), 
        width = 18, 
        height = 6
    )
    par(mfrow=c(1,3))
    AUCell_plotTSNE(
        tSNE=obj@dr$tsne@cell.embeddings, 
        exprMat=obj@data, 
        cellsAUC=cells_AUC[i,], 
        thresholds=selectedThresholds
    )
    dev.off()
}


cells_AUC <- readRDS("/mnt/raid62/Lung_cancer_10x/final_plots/02_rds/02_aucell/cells_AUC.rds")
cells_assignment <- readRDS("/mnt/raid62/Lung_cancer_10x/final_plots/02_rds/02_aucell/cells_assignment.rds")

selectedThresholds <- getThresholdSelected(cells_assignment)

tsne = read.csv("/mnt/raid62/Lung_cancer_10x/final_plots/02_rds/02_aucell/tsne.csv") 

for(i in names(selectedThresholds)) {
    print(i)
    auc <- as.data.frame(getAUC(cells_AUC)[i, ])
    colnames(auc) <- "AUC"
    
    data = as.data.frame(tsne)
    data$AUC <- auc[rownames(data), "AUC"]
    
    data$AUC[data$AUC < selectedThresholds[i]] <- 0
    
    # print(head(data))
    # data = data[order(data$AUC, decreasing = F), ]
    
    p <- ggplot(
        data, aes(x=tSNE_1, y=tSNE_2, color=AUC, alpha = 0.5)
    ) + geom_point_rast(size = 0.1) +
        theme_bw() +
        theme(
            plot.title = element_text(hjust = 0.5, face="bold", size = 20), 
            legend.position = "none",
            axis.text = element_text(size = 12),
            axis.title = element_text(size = 15)
        ) +
        scale_color_gradient(low = "lightgrey", high = "blue") +
        labs(title = paste0(i, "(AUC > ", round(selectedThresholds[i], 2), ")"))
    
    ggsave(
        filename = paste0("/mnt/raid62/Lung_cancer_10x/final_plots/01_first_plots/aucell/", i, "_tsne.pdf"),
        plot = p,
        width = 6,
        height = 6,
        dpi = 600
    )
}



# library(stringr)
# 
# 
# 
# colnames(expr) = sapply(colnames(expr), function(x) {
#     str_replace(str_replace(x, "^X", ""), "\\.", "-")
# })
# meta <- read.csv("/mnt/raid62/Lung_cancer_10x/final_plots/02_rds/meta.csv", row.names = 1)
# 
# obj <- CreateSeuratObject(
#     expr,
#     meta.data = meta[colnames(expr), ]
# )
# 
# 
# saveRDS(expr, "/mnt/raid62/Lung_cancer_10x/final_plots/02_rds/all_cell_expr.rds")
