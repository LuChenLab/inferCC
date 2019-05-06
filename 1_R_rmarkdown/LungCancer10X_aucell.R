## load expression data

library(here)
library(AUCell)
library(openxlsx)
library(GSEABase)

obj =  readRDS(here("00_data_ingest", "04_rds_generated", "seurat_obj_selected_patients_base.rds"))

## read markers
markers <- read.xlsx(here("00_data_ingest", "03_annotation_csv", "gene_markers.xlsx"), sheet = 3)


set.seed(1024)

expMatrix <- obj@raw.data[, sample(1:ncol(obj@raw.data), round(0.1 * ncol(obj@raw.data)))]

pdf(file = here("01_figures/aucell", "cells_rankings_All.pdf"))
cells_rankings <- AUCell_buildRankings(as.matrix(data), nCores=1, plotStats=TRUE)
dev.off()

saveRDS(cells_rankings, file=here("00_data_ingest", "04_rds_generated", "AUCell_cells_rankings_all.rds"))


cells_rankings <- readRDS(here("00_data_ingest", "04_rds_generated", "AUCell_cells_rankings.rds"))


#### generate gene set 
geneSets = c()
for (i in unique(markers$Cells)) {
    
    if (!is.na(i)) {
        print(i)
        geneSets <- c(geneSets, GeneSet(intersect(markers[markers$Cells == i, "Markers"], rownames(data)), setName=i))
    }
}

geneSets <- GeneSetCollection(geneSets)


#### calculate auc and make auc plots
pdf(file = here("01_figures/aucell", "stats_auc_all.pdf"))

par(mfrow=c(3,3)) 
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, assign=FALSE) 
dev.off()


saveRDS(cells_assignment, file=here("00_data_ingest", "04_rds_generated", "AUCell_cells_assignment_all.rds"))
saveRDS(cells_AUC, file=here("00_data_ingest", "04_rds_generated", "AUCell_cells_auc_al.rds"))
saveRDS(expMatrix, file=here("00_data_ingest", "04_rds_generated", "AUCell_expmatrix.rds"))

tsne <- readRDS(here("00_data_ingest", "04_rds_generated", "seurat_harmony_tsne.rds"))


selectedThresholds <- getThresholdSelected(cells_assignment)

selectedThresholds[2] <-  0.25


for (i in 1:nrow(cells_AUC)) {
    print(cells_AUC[i, ]@NAMES)
    png(here("01_figures/aucell/final", paste(gsub("[ /]", "_", cells_AUC[i, ]@NAMES), ".png", sep = "")), res = 600, width = 18, height = 6, units = "in")
    par(mfrow=c(1, 3))
    AUCell_plotTSNE(tSNE=tsne$Y, cellsAUC=cells_AUC[i, ], exprMat = data)
    dev.off()
}


#####################
marker <- read.table(here("00_data_ingest/03_annotation_csv", "Single_cell_markers.txt"), header = TRUE, sep = "\t")
temp <- marker[marker$tissueType == "Lung" & marker$speciesType != "Human", ]
marker <- marker[marker$tissueType == "Lung" & marker$speciesType == "Human", ]

temp <- temp[!temp$cellName %in% marker$cellName,]
marker <- rbind(marker, temp)


geneSets = c()

for(i in 1:nrow(marker)) {
    res = list()
    print(as.character(marker[i, "cellName"]))
    
    temp <- intersect(strsplit(as.character(marker[i, "geneSymbol"]), ", ")[[1]], cells_rankings@NAMES)
    
    if (length(temp) > 0) {
        geneSets <- c(geneSets, GeneSet(temp, setName=paste(marker[i, "cellName"], marker[i, "PMID"], sep = "_")))
    }
    
}


geneSets <- GeneSetCollection(geneSets)

# geneSets <- subsetGeneSets(geneSets, rownames(expMatrix)) 


#### calculate auc and make auc plots
pdf(file = here("01_figures/aucell", "stats_auc_cellmarker.pdf"))

par(mfrow=c(3,3)) 
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, assign=FALSE) 
dev.off()


saveRDS(cells_assignment, file=here("00_data_ingest", "04_rds_generated", "AUCell_cells_assignment_cell_marker.rds"))
saveRDS(cells_AUC, file=here("00_data_ingest", "04_rds_generated", "AUCell_cells_auc_cellmarker.rds"))

tsne <- readRDS(here("00_data_ingest", "04_rds_generated", "seurat_harmony_tsne.rds"))


for (i in 1:nrow(cells_AUC)) {
    print(cells_AUC[i, ]@NAMES)
    tiff(here("01_figures/aucell", paste("tsne_cellmarker_", gsub("[\\s/]", "_", cells_AUC[i, ]@NAMES), ".tiff", sep = "")), res = 600, compression = "lzw", width = 18, height = 6, units = "in")
    par(mfrow=c(1, 3))
    AUCell_plotTSNE(tSNE=tsne$Y, cellsAUC=cells_AUC[i, ], exprMat = expMatrix)
    dev.off()
}

gc()


#####
cells <- c()
genes <- c()
for (i in geneSets) {
    genes <- c(genes, i@geneIds)
    cells <- c(cells, rep(i@setName, length(i@geneIds)))
}

genes <- data.frame(cells=cells, markers=genes)

write.csv(genes, file = here("00_data_ingest", "03_annotation_csv", "cellmarkers_used_marker.csv"))


#####  测试所有表达量
raw_counts <- here("00_data_ingest", "00_raw_data", "raw_counts_clean.csv.gz")

data <- read.csv(gzfile(raw_counts), header=T)

data <- fread(paste("zcat", raw_counts))
data <- as.matrix(data)
row.names(data) <- data[, 1]

rownames(data) <- data[, 1]

pdf(file = here("01_figures/aucell", "cells_rankings_all.pdf"))
cells_rankings <- AUCell_buildRankings(expMatrix, nCores=8, plotStats=TRUE)
dev.off()



####
cells_rankings <- readRDS(here("00_data_ingest", "04_rds_generated", "AUCell_cells_rankings.rds"))
expMatrix <- readRDS(here("00_data_ingest", "04_rds_generated", "AUCell_expmatrix.rds"))
marker <- read.table(here("00_data_ingest/03_annotation_csv", "all_cell_markers_extract.txt"), header = FALSE, sep = "\t")



geneSets = c()
for(i in 1:nrow(marker)) {
    res = list()
    print(as.character(marker[i, 1]))
    
    temp <- intersect(strsplit(as.character(marker[i, 2]), ", ")[[1]], cells_rankings@NAMES)
    
    if (length(temp) > 0) {
        geneSets <- c(geneSets, GeneSet(temp, setName=as.character(marker[i, 1])))
    }
    
}

geneSets <- GeneSetCollection(geneSets)

pdf(file = here("01_figures/aucell", "stats_auc_allcellmarker.pdf"))

par(mfrow=c(3,3)) 
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, assign=FALSE) 
dev.off()


saveRDS(cells_assignment, file=here("00_data_ingest", "04_rds_generated", "AUCell_cells_assignment_all_cell_marker.rds"))
saveRDS(cells_AUC, file=here("00_data_ingest", "04_rds_generated", "AUCell_cells_auc_all_cellmarker.rds"))

tsne <- readRDS(here("00_data_ingest", "04_rds_generated", "seurat_harmony_tsne.rds"))


for (i in 1:nrow(cells_AUC)) {
    print(cells_AUC[i, ]@NAMES)
    tiff(here("01_figures/aucell", paste("allcellmarker_",  cells_AUC[i, ]@NAMES, ".tiff", sep = "")), res = 600, compression = "lzw", width = 18, height = 6, units = "in")
    par(mfrow=c(1, 3))
    AUCell_plotTSNE(tSNE=tsne$Y, cellsAUC=cells_AUC[i, ], exprMat = expMatrix)
    dev.off()
}