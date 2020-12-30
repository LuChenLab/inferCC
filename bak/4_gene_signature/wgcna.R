library(ggplot2)
library(dplyr)
library(reshape2)
library(WGCNA)
library(Seurat)
library(stringr)
library(openxlsx)

if(str_detect(packageVersion("Seurat"), "^3.\\d.\\d")) {
    detach("package:Seurat")
    devtools::install("Seurat", version="2.3.4")
    library("Seurat")
}

set.seed(1)
setwd("/mnt/raid62/Lung_cancer_10x/02_figures_each_cell/Basal")

obj <- readRDS("Basal_cells.rds")
r = gzfile("paga/normalized_counts.csv.gz")
data = read.csv(r, row.names = 1)

genes = read.xlsx("annotation_results_by_stage.xlsx", rowNames = T)
genes = genes[genes$p_val < 0.05 & genes$avg_logFC > 0,]


# set global variable
workingDir = getwd()

setwd(workingDir);
options(stringsAsFactors = FALSE);
enableWGCNAThreads()
